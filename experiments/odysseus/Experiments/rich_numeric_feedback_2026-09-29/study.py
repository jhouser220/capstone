"""New numerical/LLM pilots. Oracle outcomes separated from the fitting interface."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']: os.environ[k]='1'
import json,csv,hashlib,sys,time,warnings
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm,spearmanr,rankdata
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel,RBF,WhiteKernel
from sklearn.metrics import brier_score_loss,log_loss
from xgboost import XGBClassifier,XGBRegressor
import api
R=Path(__file__).resolve().parent;P=R.parents[1];DATA=api.DATA
FEATURE=Path(os.environ.get('BDA_FEATURE_FILE', P/'data/achilles_normalized.npz')).expanduser().resolve()
if not FEATURE.is_file():
 raise FileNotFoundError('Build features with prepare_features.py or set BDA_FEATURE_FILE. See experiments/odysseus/README.md.')
if os.environ.get('BDA_FEATURE_SHA256') and hashlib.sha256(FEATURE.read_bytes()).hexdigest()!=os.environ['BDA_FEATURE_SHA256']:
 raise ValueError('Feature checksum mismatch')
write=api.write_json
CACHE={}
def dataset(ds):
 if ds in CACHE:return CACHE[ds]
 rows=list(csv.DictReader((DATA/f'ground_truth_{ds}.csv').open()));scores={r['Gene']:float(r['Score']) for r in rows}
 z=np.load(FEATURE,allow_pickle=False); idx=[i for i,g in enumerate(z['genes']) if g in scores];genes=z['genes'][idx];x=z['values'][idx].astype(float)
 cache=R/(ds+'_geometry.npz')
 if cache.exists(): geom=np.load(cache); f=geom['f'];central=geom['central']
 else:
  pca=PCA(n_components=16,random_state=4201,svd_solver='randomized');f=pca.fit_transform(x);f/=f.std(0);central=x@x.mean(0);np.savez_compressed(cache,f=f,central=central,components=pca.components_,variance=pca.explained_variance_ratio_,genes=genes)
 hits=set(map(str,np.load(DATA/f'topmovers_{ds}.npy',allow_pickle=True)))
 oracle={'y':np.array([int(g in hits) for g in genes]),'effect':np.array([scores[g] for g in genes])}
 obj=(dict(ds=ds,genes=genes,x=x,f=f,central=central,task=json.loads((DATA/'task_prompts'/f'{ds}.json').read_text())),oracle);CACHE[ds]=obj;return obj

def initial(data,seed):
 rng=np.random.default_rng(seed);prior=np.argsort(-data['central'],kind='stable');ini=np.r_[prior[:64],rng.choice(prior[64:],192,replace=False)];rng.shuffle(ini);return ini.tolist()
def reveal(ids,oracle):return {'ids':list(ids),'y':oracle['y'][ids].copy(),'effect':oracle['effect'][ids].copy()}
def optimize(obj,initial_theta,bounds):
 r=minimize(obj,initial_theta,jac=True,bounds=bounds,method='L-BFGS-B',options={'maxiter':35});return r.x,r.fun

def fit_predict(data,h,q,method,diagnostics=False):
 """No oracle argument: receives ONLY revealed outcomes, public features and query IDs."""
 ids=np.array(h['ids']);q=np.array(q);x=data['f'];y=np.array(h['y']);eff=np.abs(h['effect']);extras={};model=None
 if method=='prior':p=rankdata(data['central'][q])/len(q);pred=p
 elif method=='knn':
  sim=data['x'][q]@data['x'][ids].T;nn=np.argsort(-sim,axis=1,kind='stable')[:,:20];p=(y[nn].sum(1)+1)/22;pred=eff[nn].mean(1)
 elif method.startswith('xgb'):
  common=dict(n_estimators=100,max_depth=3,learning_rate=.05,subsample=1,colsample_bytree=1,n_jobs=1,random_state=0,tree_method='hist',reg_lambda=5)
  if method=='xgb_hit':model=XGBClassifier(**common,objective='binary:logistic');model.fit(x[ids],y);p=model.predict_proba(x[q])[:,1];pred=p
  else:
   model=XGBRegressor(**common,objective='reg:quantileerror',quantile_alpha=np.array([.1,.5,.9]));model.fit(x[ids],eff);quant=np.sort(model.predict(x[q]),axis=1);pred=quant[:,1];p=pred;extras['interval']=quant[:,[0,2]]
  if diagnostics and method=='xgb_hit':
   import xgboost as xgb
   contrib=model.get_booster().predict(xgb.DMatrix(x[q]),pred_contribs=True);extras['contributions']=contrib
 elif method.startswith('gp'):
  # ARD changes diagonal metric over fixed standardized PCA components. No encoder training.
  ard=method=='gp_ard';kernel=ConstantKernel(1,(.05,10))*RBF(np.ones(16)*4 if ard else 4,(.3,50))+WhiteKernel(.5,(.05,3))
  model=GaussianProcessRegressor(kernel=kernel,normalize_y=True,optimizer=optimize,random_state=0,alpha=1e-5)
  with warnings.catch_warnings():warnings.simplefilter('ignore');model.fit(x[ids],eff)
  pred,std=model.predict(x[q],return_std=True)
  # Screen threshold estimated ONLY from revealed labels/effects; imperfect for asymmetric thresholds.
  threshold=(np.max(eff[y==0])+np.min(eff[y==1]))/2 if np.any(y==1) and np.any(y==0) else np.quantile(eff,.95)
  p=norm.sf((threshold-pred)/np.maximum(std,1e-6));extras.update(std=std,interval=np.stack([pred-1.28155*std,pred+1.28155*std],1),kernel=str(model.kernel_),threshold=float(threshold))
 else:raise ValueError(method)
 if method.startswith('gp'):extras['_model']=model
 return dict(p=np.array(p),effect=np.array(pred),**extras)

def choose(pred,q,n=32):return np.array(q)[np.argsort(-pred['p'],kind='stable')[:n]].tolist()
def diagnostic(pred,q,oracle,method):
 y=oracle['y'][q];e=np.abs(oracle['effect'][q]);out={'effect_spearman':float(spearmanr(pred['effect'],e).statistic)}
 if method not in ['prior','xgb_effect']:
  pp=np.clip(pred['p'],1e-5,1-1e-5);out.update(brier=float(np.mean((pp-y)**2)),log_loss=float(log_loss(y,pp,labels=[0,1])))
 if 'interval'in pred:out.update(coverage80=float(np.mean((e>=pred['interval'][:,0])&(e<=pred['interval'][:,1]))),interval_width=float(np.mean(pred['interval'][:,1]-pred['interval'][:,0])))
 return out

def experiment(name,datasets,seeds,methods,steps=1):
 folder=R/name;folder.mkdir(exist_ok=False);write(folder/'protocol.json',dict(time=api.stamp(),datasets=datasets,seeds=seeds,methods=methods,steps=steps,initial=256,batch=32,pool='all feature-covered untested genes; unselected retained',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
 results=[]
 for ds in datasets:
  data,oracle=dataset(ds)
  for seed in seeds:
   ini=initial(data,seed)
   for method in methods:
    ids=list(ini);cum=0;curve=[];rec=[];start=time.monotonic()
    for step in range(steps):
     q=np.array([i for i in range(len(data['genes'])) if i not in set(ids)]);h=reveal(ids,oracle)
     pred=fit_predict(data,h,q,method);selected=choose(pred,q)
     np.savez_compressed(folder/f'{ds}_{seed}_{method}_{step}_predictions.npz',query=q,p=pred['p'],effect=pred['effect'],history=np.array(ids),selected=selected,**{k:v for k,v in pred.items() if k in ['std','interval']})
     assert len(set(selected))==32 and not set(selected)&set(ids)
     nh=int(oracle['y'][selected].sum());cum+=nh;curve.append(cum);rec.append(dict(step=step,hits=nh,selected=selected,diagnostics=diagnostic(pred,q,oracle,method),kernel=pred.get('kernel'),threshold=pred.get('threshold')));ids+=selected
    row=dict(dataset=ds,seed=seed,method=method,hits=cum,tests=32*steps,hit_rate=cum/(32*steps),discovery_auc=float(np.mean(curve)),curve=curve,rounds=rec,runtime=time.monotonic()-start);results.append(row);write(folder/f'{ds}_{seed}_{method}.json',row);print(name,ds,seed,method,cum,flush=True)
 write(folder/'results.json',results);aggregate(folder,results);return results

def aggregate(folder,results):
 summary={}
 for m in sorted({r['method'] for r in results}):
  rr=[r for r in results if r['method']==m];summary[m]=dict(hits=sum(r['hits'] for r in rr),tests=sum(r['tests'] for r in rr),by_dataset={ds:sum(r['hits'] for r in rr if r['dataset']==ds) for ds in sorted({r['dataset'] for r in rr})},cases=[r['hits'] for r in rr])
 write(folder/'summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':
 cfg=json.loads(Path(sys.argv[1]).read_text());experiment(**cfg)
