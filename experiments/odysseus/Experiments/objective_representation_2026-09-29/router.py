from study2 import *
EXPERTS=['gp_iso','full:effect_median']
def validate(data,oracle,ini):
 train=ini[:192];valid=ini[192:];h=base.reveal(train,oracle);pred={m:predict(data,h,valid,m) for m in EXPERTS}
 # These labels are already observed campaign history; NEVER future candidate labels.
 out={m:dict(top16_hits=int(oracle['y'][base.choose(pred[m],valid,16)].sum()),mae=float(np.mean(np.abs(pred[m]['effect']-np.abs(oracle['effect'][valid])))),spearman=float(spearmanr(pred[m]['effect'],np.abs(oracle['effect'][valid])).statistic),scores=pred[m]['p'].tolist()) for m in EXPERTS}
 choices={'cv_hits':max(EXPERTS,key=lambda m:out[m]['top16_hits']),'cv_mae':min(EXPERTS,key=lambda m:out[m]['mae'])}
 return dict(train=train,valid=valid,metrics=out,choices=choices)
def run(name,datasets,seeds,steps=3):
 folder=R/name;folder.mkdir();write(folder/'protocol.json',dict(time=api.stamp(),datasets=datasets,seeds=seeds,steps=steps,validation='first192 initial history train,last64 validation; top16hits or MAE; tieGP; freeze expert before new selections',experts=EXPERTS));rows=[]
 for ds in datasets:
  data,oracle=base.dataset(ds)
  for seed in seeds:
   ini=base.initial(data,seed);v=validate(data,oracle,ini);write(folder/f'{ds}_{seed}_validation.json',v);runs={}
   # Run each actual distinct selected policy once. Routers re-use deterministic expert trajectory explicitly.
   for method in EXPERTS+['prior']:
    ids=list(ini);cum=0;curve=[];rec=[]
    for step in range(steps):
     h=base.reveal(ids,oracle);q=np.array([i for i in range(len(data['genes'])) if i not in set(ids)]);pred=predict(data,h,q,method);sel=base.choose(pred,q)
     np.savez_compressed(folder/f'{ds}_{seed}_{method.replace(":","-")}_{step}_predictions.npz',history=ids,query=q,score=pred['p'],selected=sel)
     nh=int(oracle['y'][sel].sum());cum+=nh;curve.append(cum);rec.append(dict(step=step,hits=nh,selected=sel));ids+=sel
    row=dict(dataset=ds,seed=seed,method=method,hits=cum,tests=32*steps,hit_rate=cum/(32*steps),curve=curve,discovery_auc=float(np.mean(curve)),rounds=rec);runs[method]=row;rows.append(row);print(name,ds,seed,method,cum,flush=True)
   for policy,method in v['choices'].items():rows.append(dict(runs[method],method=policy,selected_expert=method,shared_deterministic_trajectory=True))
 write(folder/'results.json',rows);base.aggregate(folder,rows)
if __name__=='__main__':run(**json.loads(Path(sys.argv[1]).read_text()))
