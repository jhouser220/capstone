"""Explicit AWS reruns with a shared spending guard and private receipt output."""
import argparse,json,os,shutil,threading,time,uuid,gzip
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import httpx
from .reproduce import run_one,install_transport,IntegrityError,load,digest
RATES={'google.gemma-4-31b':(.14,.40),'qwen.qwen3-32b-v1:0':(.15,.60)}

class Budget:
 def __init__(self,cap,out):
  self.cap=cap;self.out=Path(out);self.lock=threading.RLock();self.stop=threading.Event();self.inflight={};self.rows=[]
 def totals(self):return {'known_usage_usd':sum(x['cost'] for x in self.rows),'uncertain_billing_usd':sum(x['uncertain'] for x in self.rows),'physical_calls':len(self.rows),'inflight':len(self.inflight),'cap_usd':self.cap}
 def transport(self,actor):
  def call(request,state,step,used):
   token=os.environ.get('AWS_BEARER_TOKEN_BEDROCK')
   if not token:raise IntegrityError('AWS_BEARER_TOKEN_BEDROCK is required')
   body=json.loads(request.content)
   if body['model']!=actor or body['max_tokens']!=8192 or body['temperature']!=1.:raise IntegrityError('Serving configuration drift')
   if set(body)&{'top_p','top_k','min_p','reasoning_effort','chat_template_kwargs'}:raise IntegrityError('Unexpected sampling override')
   if shutil.disk_usage(self.out).free<10*1024**3:raise IntegrityError('Storage reserve reached')
   ri,ro=RATES[actor];ceiling=((262144 if actor.startswith('google.') else 32768)*ri+16384*ro)/1e6
   attempt=uuid.uuid4().hex
   with self.lock:
    t=self.totals()
    if self.stop.is_set() or t['known_usage_usd']+t['uncertain_billing_usd']+sum(self.inflight.values())+ceiling>self.cap:raise IntegrityError('Budget/access stop')
    self.inflight[attempt]=ceiling
   endpoint='https://bedrock-mantle.us-east-1.api.aws/openai/v1/chat/completions' if actor.startswith('google.') else 'https://bedrock-runtime.us-east-1.amazonaws.com/model/qwen.qwen3-32b-v1%3A0/invoke'
   wire=body if actor.startswith('google.') else {k:v for k,v in body.items() if k!='model'}
   row={'attempt_id':attempt,'actor':actor,'cost':0.,'uncertain':ceiling,'status_code':None}
   receipt={'request':body,'endpoint':endpoint,'wire_request':wire,'state':state,'step':step,'response':None,'status_code':None};start=time.perf_counter()
   try:
    with httpx.Client(timeout=300) as client:response=client.post(endpoint,headers={'Authorization':'Bearer '+token},json=wire)
    data=response.json();receipt.update(response=data,status_code=response.status_code);row['status_code']=response.status_code
    usage=data.get('usage',{})
    if 'prompt_tokens' in usage and 'completion_tokens' in usage:
     row.update(input_tokens=int(usage['prompt_tokens']),output_tokens=int(usage['completion_tokens']),cost=(usage['prompt_tokens']*ri+usage['completion_tokens']*ro)/1e6,uncertain=0.)
    if response.status_code!=200:
     self.stop.set();raise IntegrityError('Provider returned HTTP '+str(response.status_code))
    if data.get('model')!=actor:self.stop.set();raise IntegrityError('Provider actor mismatch')
    text=json.dumps(data.get('choices',[])).lower()
    if any(s in text for s in ['i cannot assist with','i cannot help with',"i can\'t assist with"]):self.stop.set();raise IntegrityError('Actor refusal')
    used.append(receipt)
    return httpx.Response(200,json=data)
   except BaseException:
    self.stop.set();raise
   finally:
    row['seconds']=time.perf_counter()-start
    raw=json.dumps(receipt)
    if token in raw:raise IntegrityError('Credential in provider evidence')
    p=self.out/'receipts'/f'{attempt}.json.gz';p.parent.mkdir(exist_ok=True)
    with gzip.open(p,'wt') as f:f.write(raw)
    with self.lock:
     self.inflight.pop(attempt);self.rows.append(row)
     with (self.out/'ATTEMPTS.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
     (self.out/'COSTS.json').write_text(json.dumps(self.totals(),indent=2))
  return call

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--bundle',required=True);ap.add_argument('--output',required=True);ap.add_argument('--cap-usd',type=float,required=True);ap.add_argument('--allow-paid',action='store_true');ap.add_argument('--workers',type=int,default=2);ap.add_argument('--limit',type=int);args=ap.parse_args()
 if not args.allow_paid:raise ValueError('Live execution requires --allow-paid')
 if not 0<args.cap_usd<=1:raise ValueError('This verification permits a maximum $1 incremental allocation')
 O=Path(args.output);O.mkdir(parents=True,exist_ok=False);B=Path(args.bundle)
 for rel,h in load(B/'SHA256.json').items():
  import hashlib
  if hashlib.sha256((B/rel).read_bytes()).hexdigest()!=h:raise IntegrityError('Asset hash mismatch')
 entries=[e for e in load(B/'manifest.json') if e['cell']['arm'] in ('incumbent','native_filled')]
 # Alternate arms/actors within each publication; every trajectory gets fresh requests.
 entries.sort(key=lambda e:(e['cell']['screen_id'],e['cell']['draw'],e['cell']['actor'],e['cell']['arm']))
 if args.limit:entries=entries[:args.limit]
 install_transport();budget=Budget(args.cap_usd,O);start=time.perf_counter();done=[]
 def work(e):
  try:return run_one(B,e,O,mode='fresh',live_transport=budget.transport(e['cell']['actor']))
  except BaseException as exc:return {'cell':e['cell'],'status':'failed','error':str(exc)}
 with ThreadPoolExecutor(max_workers=args.workers) as pool:
  for f in as_completed([pool.submit(work,e) for e in entries]):
   row=f.result();done.append(row)
   with (O/'ledger.jsonl').open('a') as out:out.write(json.dumps(row)+'\n')
   if len(done)%4==0 or row['status']=='failed':print(json.dumps({'done':len(done),'total':len(entries),'failed':sum(x['status']=='failed' for x in done),'seconds':time.perf_counter()-start,**budget.totals()}),flush=True)
 summary={'expected':len(entries),'passed':sum(x['status']=='passed' for x in done),'failed':sum(x['status']=='failed' for x in done),'seconds':time.perf_counter()-start,'costs':budget.totals(),'rows':done}
 (O/'SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)
 if summary['failed']:raise SystemExit(1)
if __name__=='__main__':main()
