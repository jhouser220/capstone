"""Fresh numerical reports and LLM probabilities. Calibration uses only old revealed labels."""
from study import *
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics import brier_score_loss
SMALL='google/gemini-2.5-flash-lite';BIG='google/gemini-2.5-flash'
def get_llm(data,menu,report,model,seed,role):
 candidates=[{'id':i,'gene':str(data['genes'][g])} for i,g in enumerate(menu)]
 prompt='Task: '+json.dumps(data['task'])+'\nEstimate the probability each candidate is a hit in this screen (a large measured effect, in either direction). These outcomes are unknown. Use prior biology and the report if supplied. Probabilities must reflect uncertainty; numerical predictions are not measurements.\nCandidates: '+json.dumps(candidates)+'\nReport: '+json.dumps(report)+'\nReturn JSON {"probabilities":{"0":0.15,"1":0.20,"2":0.05,...},"reason":"one short paragraph describing evidence, including uncertainty"}. Include every candidate ID as a key and a numeric probability from0 to1. Supply exactly '+str(len(menu))+' probabilities. No other genes.'
 cid=role+'_'+api.digest([prompt,model,seed])[:20];reused=(R/'calls'/('rich_'+cid)/'response.json').exists();response=api.request_llm(cid,prompt,model,seed,2000);obj=api.parse_result(response);repaired=False
 if isinstance(obj,dict) and 'probabilities' not in obj:
  import re
  content=response['choices'][0]['message']['content'];fixed=re.sub(r'(:\s*[0-9]+(?:\.[0-9]+)?)"(?=\s*[,}])',r'\1',content)
  try:obj=json.loads(fixed);repaired=fixed!=content
  except ValueError:pass
 mapping=obj.get('probabilities',{}) if isinstance(obj,dict) else {};arr=[mapping.get(str(i)) for i in range(len(menu))] if isinstance(mapping,dict) and set(mapping)=={str(i) for i in range(len(menu))} else []
 string_count=sum(isinstance(v,str) for v in arr)
 try:arr=[float(v) if isinstance(v,str) else v for v in arr]
 except ValueError:arr=[]
 valid=len(arr)==len(menu) and all(type(v) in [int,float] and np.isfinite(v) and 0<=v<=1 for v in arr)
 if not valid:raise RuntimeError('Invalid probabilities; response preserved, no silent repair: '+cid)
 return np.array(arr),dict(call_id='rich_'+cid,prompt_sha256=api.digest(prompt),reason=obj.get('reason'),syntax_repaired=repaired,numeric_strings_converted=string_count,reused_identical_response=reused)

def report_for(data,h,menu,method):
 pred=fit_predict(data,h,menu,method,diagnostics=method=='xgb_hit');xp=fit_predict(data,h,menu,'xgb_hit',diagnostics=True);kn=fit_predict(data,h,menu,'knn');quant=fit_predict(data,h,menu,'xgb_effect');hist=np.array(h['ids']);sim=data['x'][menu]@data['x'][hist].T;nearest=np.argmax(sim,axis=1)
 cases=[]
 for j,i in enumerate(menu):
  obs=[]
  # Closest success and non-hit, both measured, without choosing on candidate outcomes.
  for label in [1,0]:
   pos=np.flatnonzero(h['y']==label)
   if len(pos):
    k=int(pos[np.argmax(sim[j,pos])]);obs.append(dict(gene=str(data['genes'][hist[k]]),hit=int(h['y'][k]),effect=round(float(h['effect'][k]),5),similarity=round(float(sim[j,k]),4)))
  card=dict(id=j,numerical_probability=round(float(pred['p'][j]),5),knn_probability=round(float(kn['p'][j]),5),predicted_absolute_effect_quantiles=[round(float(x),5) for x in [quant['interval'][j,0],quant['effect'][j],quant['interval'][j,1]]],observed_neighbours=obs)
  if 'contributions'in xp:
   card['xgb_probability']=round(float(xp['p'][j]),5)
   cc=xp['contributions'][j,:-1];top=np.argsort(-np.abs(cc))[:2];card['feature_contributions_log_odds']=[dict(component=int(k+1),value=round(float(cc[k]),4)) for k in top]
  cases.append(card)
 return dict(training_count=len(hist),hits=int(h['y'].sum()),warning='Quantile intervals are uncalibrated. PCA component contributions explain the predictor, not biological causes. Nearest examples are selected, not random samples.',candidates=cases),pred

def experiment_llm(name,datasets,seeds,method,model=SMALL,steps=1):
 folder=R/name;folder.mkdir(exist_ok=True);write(folder/'protocol.json',dict(time=api.stamp(),datasets=datasets,seeds=seeds,numerical=method,model=model,steps=steps,menu='union numerical top48 and geometry top48, deduplicated, shuffled; same across one-step arms',arms=['numeric','llm_prior','llm_report','fixed_blend','calibrated_blend'],weight='wNumeric=exp(-10*BrierNumeric)/(exp(-10*BrierNumeric)+exp(-10*BrierLLM)); calibrated on32 masked initial observations',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
 rows=[]
 for ds in datasets:
  data,oracle=dataset(ds)
  for seed in seeds:
   rng=np.random.default_rng(seed);ini=initial(data,seed);h=reveal(ini,oracle)
   # Old observed calibration set is not supplied to prediction models or LLM until scores saved.
   tr=ini[:-32];val=ini[-32:];report,vpred=report_for(data,reveal(tr,oracle),val,method)
   vl,info=get_llm(data,val,report,model,seed,'calibration');write(folder/f'{ds}_{seed}_calibration_predictions.json',dict(ids=val,numerical=vpred['p'].tolist(),llm=vl.tolist(),call=info))
   bn=float(np.mean((vpred['p']-oracle['y'][val])**2));bl=float(np.mean((vl-oracle['y'][val])**2));w=float(np.exp(-10*bn)/(np.exp(-10*bn)+np.exp(-10*bl)))
   q=np.array([i for i in range(len(data['genes'])) if i not in set(ini)]);pp=fit_predict(data,h,q,method);top=choose(pp,q,48);prior=sorted(q,key=lambda i:-data['central'][i])[:48];menu=[int(i) for i in dict.fromkeys(top+prior)];rng.shuffle(menu)
   report,pred=report_for(data,h,menu,method);write(folder/f'{ds}_{seed}_visible.json',dict(history=[dict(gene=str(data['genes'][g]),hit=int(oracle['y'][g]),effect=float(oracle['effect'][g])) for g in ini],menu=[str(data['genes'][g]) for g in menu],report=report,numerical_weight=w,validation_brier=dict(numeric=bn,llm=bl)))
   lp,ci=get_llm(data,menu,'No historical observations supplied.',model,seed,'prior');lr,cr=get_llm(data,menu,report,model,seed,'report')
   scores={'numeric':pred['p'],'llm_prior':lp,'llm_report':lr,'fixed_blend':.5*pred['p']+.5*lr,'calibrated_blend':w*pred['p']+(1-w)*lr}
   decisions={m:choose({'p':v},menu) for m,v in scores.items()};write(folder/f'{ds}_{seed}_predictions_decisions.json',dict(scores={m:v.tolist() for m,v in scores.items()},menu=menu,decisions=decisions,calls=[ci,cr]))
   base=set(decisions['numeric'])
   for m,sel in decisions.items():
    added=set(sel)-base;removed=base-set(sel);row=dict(dataset=ds,seed=seed,method=m,hits=int(oracle['y'][sel].sum()),tests=32,hit_rate=float(oracle['y'][sel].mean()),added_hits=int(oracle['y'][list(added)].sum()),removed_hits=int(oracle['y'][list(removed)].sum()),changes=len(added),brier=float(np.mean((scores[m]-oracle['y'][menu])**2)),numerical_weight=w);rows.append(row)
   print(name,ds,seed,[(r['method'],r['hits']) for r in rows[-5:]],flush=True)
 write(folder/'results.json',rows);aggregate(folder,rows)
if __name__=='__main__':experiment_llm(**json.loads(Path(sys.argv[1]).read_text()))
