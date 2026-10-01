from pathlib import Path
import sys
R=Path(__file__).resolve().parent;OLD=R.parent/'rich_numeric_feedback_2026-09-29';sys.path.insert(0,str(OLD))
import study as base
from study import np,json,time,hashlib,api,P
from xgboost import XGBClassifier,XGBRegressor,XGBRanker
from scipy.stats import spearmanr
write=api.write_json

def predict(data,h,q,method):
 if method in ['prior','gp_iso','knn']:return base.fit_predict(data,h,q,method)
 feature,obj=method.split(':');x=data['f'] if feature=='pca' else data['x'];ids=h['ids'];y=h['y'];mag=np.abs(h['effect']);kw=dict(n_estimators=100,max_depth=3,learning_rate=.05,reg_lambda=5,tree_method='hist',n_jobs=1,random_state=0,subsample=1,colsample_bytree=1)
 if obj=='hit_logistic':
  model=XGBClassifier(**kw,objective='binary:logistic');model.fit(x[ids],y);s=model.predict_proba(x[q])[:,1]
 elif obj=='hit_rank':
  model=XGBRanker(**kw,objective='rank:ndcg',lambdarank_pair_method='topk',lambdarank_num_pair_per_sample=32);model.fit(x[ids],y,group=[len(ids)]);s=model.predict(x[q])
 else:
  target=y if obj=='hit_squared' else mag
  model=XGBRegressor(**kw,objective='reg:quantileerror' if obj=='effect_median' else 'reg:squarederror',**({'quantile_alpha':.5} if obj=='effect_median' else {}));model.fit(x[ids],target);s=model.predict(x[q])
 return {'p':s,'effect':s}

def experiment(name,datasets,seeds,methods,steps=1,modes=['updated'],initial_mode='mixed'):
 folder=R/name;folder.mkdir();cfg=dict(name=name,datasets=datasets,seeds=seeds,methods=methods,steps=steps,modes=modes,initial_mode=initial_mode,initial=256,batch=32,pool='all feature-covered untested; unselected retained',time=api.stamp(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());write(folder/'protocol.json',cfg);(folder/'code_snapshot.py').write_text(Path(__file__).read_text());rows=[]
 for ds in datasets:
  data,oracle=base.dataset(ds)
  for seed in seeds:
   ini=base.initial(data,seed) if initial_mode=='mixed' else np.random.default_rng(seed).choice(len(data['genes']),256,replace=False).tolist()
   for method in methods:
    for mode in modes:
     ids=list(ini);cum=0;curve=[];rec=[];t0=time.monotonic()
     for step in range(steps):
      hist=ini if mode=='frozen' else ids;h=base.reveal(hist,oracle)
      if mode=='shuffle_new' and len(hist)>256:
       perm=np.random.default_rng(seed+step).permutation(len(hist)-256)+256;h['y'][256:]=h['y'][perm];h['effect'][256:]=h['effect'][perm]
      q=np.array([i for i in range(len(data['genes'])) if i not in set(ids)]);pred=predict(data,h,q,method);sel=base.choose(pred,q)
      assert len(set(sel))==32 and not set(sel)&set(ids)
      stem=f'{ds}_{seed}_{method.replace(":","-")}_{mode}_{step}'
      np.savez_compressed(folder/(stem+'_predictions.npz'),history=hist,history_y=h['y'],history_effect=h['effect'],query=q,score=pred['p'],selected=sel)
      nh=int(oracle['y'][sel].sum());cum+=nh;curve.append(cum);diag={'rank_effect_spearman':float(spearmanr(pred['p'],np.abs(oracle['effect'][q])).statistic)}
      if method.endswith(('hit_logistic','hit_squared')) or method in ['gp_iso','knn']:
       pp=np.clip(pred['p'],1e-6,1-1e-6);diag['brier']=float(np.mean((pp-oracle['y'][q])**2));diag['log_loss']=float(-np.mean(oracle['y'][q]*np.log(pp)+(1-oracle['y'][q])*np.log(1-pp)))
      rec.append(dict(step=step,hits=nh,selected=sel,diagnostics=diag));ids+=sel
     row=dict(dataset=ds,seed=seed,method=method,mode=mode,initial_mode=initial_mode,hits=cum,tests=32*steps,hit_rate=cum/(32*steps),curve=curve,discovery_auc=float(np.mean(curve)),rounds=rec,runtime_seconds=time.monotonic()-t0);rows.append(row);write(folder/f'{ds}_{seed}_{method.replace(":","-")}_{mode}.json',row);print(name,ds,seed,method,mode,cum,flush=True)
 write(folder/'results.json',rows);summary={}
 for m in methods:
  for mode in modes:
   rr=[r for r in rows if r['method']==m and r['mode']==mode];summary[m+'/'+mode]=dict(hits=sum(r['hits'] for r in rr),tests=sum(r['tests'] for r in rr),cases=[r['hits'] for r in rr],by_dataset={ds:sum(r['hits'] for r in rr if r['dataset']==ds) for ds in datasets})
 write(folder/'summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':experiment(**json.loads(Path(sys.argv[1]).read_text()))
