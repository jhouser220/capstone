"""Counterfactual one-step replay: change history, keep available actions fixed."""
from run import *
def main():
 folder=ROOT/'phase5_matched_pool';folder.mkdir(exist_ok=True);write(folder/'protocol.json',dict(time=old.api.stamp(),datasets=['IL2','Carnevale22_Adenosine','Sanchez21_down'],seeds=[8101,8102,8103],history_size=256,batch=32,query_pool='exclude UNION of both histories in both arms',methods=['prior','gp_iso','full:effect_median'],steps=1,interpretation='adaptive counterfactual replay, not an online campaign or independent confirmation',source='LLMNN prior/feedback separation, projection to matched history information',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()));rows=[]
 for ds in ['IL2','Carnevale22_Adenosine','Sanchez21_down']:
  data,o=dataset(ds)
  for seed in [8101,8102,8103]:
   histories=dict(mixed=old.base.initial(data,seed),random=np.random.default_rng(seed).choice(len(data['genes']),256,replace=False).tolist());excluded=set(histories['mixed'])|set(histories['random']);q=np.array([i for i in range(len(data['genes'])) if i not in excluded])
   for mode,ini in histories.items():
    h=old.base.reveal(ini,o)
    for method in ['prior','gp_iso','full:effect_median']:
     pred=old.predict(data,h,q,method);sel=old.base.choose(pred,q);stem=f'{ds}_{seed}_{mode}_{method.replace(":","-")}'
     np.savez_compressed(folder/(stem+'.npz'),history=ini,excluded=sorted(excluded),query=q,score=pred['p'],selected=sel)
     row=dict(dataset=ds,seed=seed,history_mode=mode,method=method,selected=sel,hits=int(o['y'][sel].sum()),tests=32,initial_hits=int(h['y'].sum()),query_n=len(q));rows.append(row);print(ds,seed,mode,method,row['hits'],flush=True)
 write(folder/'results.json',rows)
if __name__=='__main__':main()
