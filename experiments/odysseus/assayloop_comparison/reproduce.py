"""Run complete campaigns from explicit frozen assets and checked LLM receipts.

Run as a module: python -m experiments.odysseus.assayloop_comparison.reproduce
Offline mode reconstructs decisions; recorded selections are comparison targets.
"""
import argparse,contextvars,gzip,hashlib,json,os,time,traceback
from pathlib import Path
import numpy as np
from .adapter import SourceState,run_campaign
from .assets import HistoricalBank

ACTIVE=contextvars.ContextVar('comparison_transport',default=None)
class IntegrityError(BaseException):pass

def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,allow_nan=True).encode()).hexdigest()
def load(path):
 p=Path(path)
 return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def observations(obs):return [{'candidate':o.candidate,'label':dict(o.label),'metadata':dict(o.metadata)} for o in obs]
def state(kw):return {'history':[{'round':x.step,'accepted':x.acquired_batch,'observations':observations(x.new_observations)} for x in kw['history']],'candidates':kw['candidates'],'context':kw['task_context']}

def install_transport():
 import logging
 logging.disable(logging.CRITICAL)
 import openai
 import assayloop.tracing as tracing
 if getattr(openai.OpenAI,'_comparison_bound',False):return
 real=openai.OpenAI
 def factory(**kw):
  client=ACTIVE.get()
  if client is None:raise IntegrityError('Missing explicitly bound transport')
  return real(**kw,http_client=client)
 factory._comparison_bound=True
 openai.OpenAI=factory
 # Safe evidence is written below; upstream tracing includes client configuration.
 tracing.log_completion_call=lambda **kw:None

class ReceiptActor:
 def __init__(self,actor,records,live_transport=None):
  from assayloop.acquisitions.llm_single_acq import LLMSingleAcquisition
  from assayloop.llm.client import LLMClientConfig
  self.actor=actor;self.records=records;self.live=live_transport;self.events=[];self.step=0
  self.original=LLMSingleAcquisition(llm=LLMClientConfig(provider='vllm',model=actor,
   base_url='https://bedrock-mantle.us-east-1.api.aws/openai/v1',api_key='externally-bound-transport',
   max_tokens=8192,temperature=1.,disable_thinking=False,timeout=300,max_retries=8,retry_on_empty=False),oversample_factor=1.)
 def reset(self):self.step=0;self.events=[];self.original.reset()
 def suggest(self,**kw):
  import httpx
  current=state(kw);rec=self.records[self.step] if self.live is None else None
  if rec and digest(current)!=digest(rec['state']):raise IntegrityError(f'Language history/context/candidate mismatch at step {self.step+1}')
  used=[]
  def handle(request):
   body=json.loads(request.content)
   if self.live is not None:
    return self.live(request,current,self.step+1,used)
   if len(used)>=len(rec['events']):raise IntegrityError('More requests than saved receipts')
   saved=rec['events'][len(used)]
   binding=saved.get('binding',{})
   expected=binding.get('request',saved.get('request'))
   if expected is None:raise IntegrityError('Unknown receipt request schema')
   if digest(body)!=digest(expected):raise IntegrityError('Exact request body mismatch')
   response=saved['response'];status=saved.get('status_code',200)
   if isinstance(response,str):response=json.loads(response)
   used.append({'request':body,'response':response,'status_code':status,'source_digest':digest(saved)})
   return httpx.Response(status,json=response)
  with httpx.Client(transport=httpx.MockTransport(handle)) as client:
   token=ACTIVE.set(client)
   try:proposal=self.original.suggest(**kw)
   finally:ACTIVE.reset(token)
  if rec and (len(used)!=len(rec['events']) or proposal!=rec['proposal']):raise IntegrityError('Receipt count or original parser output mismatch')
  self.events.append({'step':self.step+1,'state':current,'proposal':proposal,'events':used})
  self.step+=1
  return proposal

def metrics(batches,screen):
 library=set(screen.genes);truth={g for g,h in zip(screen.genes,screen.hits) if h};selected=[g for b in batches for g in b]
 assert len(set(selected))==len(selected)
 hits=len(set(selected)&truth);measured=sum(g in library for g in selected)
 return dict(unique_hits=hits,selected=len(selected),measured=measured,unavailable=len(selected)-measured,requested=1000,
  measured_precision=hits/measured if measured else 0.,ef=hits/((1000-len(selected)+measured)*len(truth)/len(library)))

def run_one(bundle,entry,out,mode='replay',live_transport=None):
 from assaybench.data.screens import ScreenRecord
 from assayloop.models.amortized_ranker import AmortizedRankerModel
 from assayloop.tasks.gene_batch import AssayBenchGeneBatchTask
 import torch
 torch.set_num_threads(1)
 start=time.perf_counter();B=Path(bundle);cell=entry['cell'];sid=cell['screen_id'];original=load(B/entry['file'])
 assert hashlib.sha256((B/entry['file']).read_bytes()).hexdigest()==entry['sha256']
 sc=ScreenRecord(**load(B/'screens'/f'{sid}.json'));desc=np.load(B/'screens'/f'{sid}.npy',allow_pickle=False)
 source=AmortizedRankerModel(checkpoint=B/'checkpoint',ckpt_file='model_last.pt',device='cpu',desc_emb_by_name={sid:desc})
 for p in source.net.parameters():p.requires_grad_(False)
 task=AssayBenchGeneBatchTask(sc,universe_genes=load(B/f"UNIVERSE_{entry['fold']}.json"),seed=0)
 arm={'source':'assayformer_only','scalar':'historical_only'}.get(cell['arm'],cell['arm'])
 bank=HistoricalBank.load(B/'banks'/f'{sid}.npz') if arm in ('incumbent','historical_only') else None
 actor=ReceiptActor(cell['actor'],original['language'],live_transport) if cell['actor']!='none' else None
 result,audit=run_campaign(arm=arm,task=task,source_model=source,screen_id=sid,draw=cell['draw'],source_state=SourceState(source,desc),bank=bank,actor_acquisition=actor,allow_live_calls=actor is not None)
 rows=[{'round':x.step,'accepted':x.acquired_batch,'observations':observations(x.new_observations)} for x in result.history]
 if mode=='replay':
  for actual,expected in zip(rows,original['rounds']):
   if digest(actual)!=digest(expected):raise IntegrityError(f"Decision/outcome mismatch: {cell}, round {actual['round']}")
  assert len(rows)==len(original['rounds'])==10
 m=metrics([x['accepted'] for x in rows],sc)
 if mode=='replay':
  assert m['unique_hits']==original['primary']['unique_hits'];assert abs(m['ef']-original['primary']['ef'])<1e-10
 summary=dict(cell=cell,cohort=entry['cohort'],publication_id=entry['publication_id'],mode=mode,status='passed',metrics=m,seconds=time.perf_counter()-start,
  verified_decisions=len(rows),verified_language_requests=sum(len(x['events']) for x in actor.events) if actor else 0)
 output=Path(out)/Path(entry['file']).name
 output.parent.mkdir(parents=True,exist_ok=True)
 with gzip.open(output,'wt') as f:json.dump(dict(summary=summary,rounds=rows,audit=audit,language=actor.events if actor else []),f)
 return summary

def worker(args):
 install_transport()
 try:return run_one(*args)
 except BaseException as e:return dict(cell=args[1]['cell'],status='failed',error=str(e),traceback=traceback.format_exc())

def main():
 from concurrent.futures import ProcessPoolExecutor,as_completed
 import multiprocessing
 ap=argparse.ArgumentParser();ap.add_argument('--bundle',required=True);ap.add_argument('--output',required=True);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--limit',type=int);ap.add_argument('--arms',nargs='+');ap.add_argument('--screen');args=ap.parse_args()
 B=Path(args.bundle);O=Path(args.output);O.mkdir(parents=True,exist_ok=True)
 if (O/'SUMMARY.json').exists():raise ValueError('Use a fresh output directory')
 for rel,h in load(B/'SHA256.json').items():
  if hashlib.sha256((B/rel).read_bytes()).hexdigest()!=h:raise IntegrityError('Asset hash mismatch: '+rel)
 entries=[x for x in load(B/'manifest.json') if (not args.arms or x['cell']['arm'] in args.arms) and (not args.screen or x['cell']['screen_id']==args.screen)]
 if args.limit:entries=entries[:args.limit]
 start=time.perf_counter();done=[]
 with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as pool:
  futures=[pool.submit(worker,(str(B),e,str(O))) for e in entries]
  for f in as_completed(futures):
   row=f.result();done.append(row)
   with (O/'ledger.jsonl').open('a') as out:out.write(json.dumps(row)+'\n')
   if row['status']!='passed':print(json.dumps(row),flush=True)
   elif len(done)%8==0:print(json.dumps({'completed':len(done),'total':len(entries),'seconds':time.perf_counter()-start}),flush=True)
 summary={'passed':sum(x['status']=='passed' for x in done),'failed':sum(x['status']!='passed' for x in done),'expected':len(entries),'seconds':time.perf_counter()-start,'fresh_provider_calls':0,'rows':done}
 (O/'SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)
 if summary['failed']:raise SystemExit(1)
if __name__=='__main__':main()
