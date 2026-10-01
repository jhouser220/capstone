from study import *
from scipy.linalg import solve_triangular

def acquire(data,h,q,pred,arm):
 if arm=='gp_iso':return choose(pred,q),{}
 gp=pred['_model'];x=data['f'][q];xt=data['f'][h['ids']];kernel=gp.kernel_.k1;noise=float(gp.kernel_.k2.noise_level);v=solve_triangular(gp.L_,kernel(xt,x),lower=True);var=np.maximum(kernel.diag(x)-np.sum(v*v,axis=0),1e-9);original=var.copy();factors=[];chosen=np.argsort(-pred['p'],kind='stable')[:24].tolist();diagnostics={'noise':noise,'initial_latent_variance_mean':float(var.mean())}
 if arm=='gp_uncertainty8':
  order=np.argsort(-var,kind='stable');chosen+= [int(i) for i in order if i not in chosen][:8]
 else:
  # Exact posterior Gaussian information increment about latent candidate function;
  # noisy observations; no knowledge of selected future outcomes needed.
  for j in range(32):
   if j>=24:
    scores=.5*np.log1p(np.maximum(var,0)/noise);scores[chosen]=-np.inf;chosen.append(int(np.argmax(scores)))
   i=chosen[j];cov=kernel(x,x[i:i+1]).ravel()-v.T@v[:,i]
   if factors:
    ff=np.array(factors).T;cov-=ff@ff[i,:]
   col=cov/np.sqrt(max(var[i],0)+noise);factors.append(col);var=np.maximum(var-col*col,0)
  diagnostics['remaining_latent_variance_mean']=float(var.mean())
 diagnostics['exploration_ids']=[int(q[i]) for i in chosen[24:]]
 return [int(q[i]) for i in chosen],diagnostics

def run(name,datasets,seeds,steps=3):
 folder=R/name;folder.mkdir(exist_ok=True);write(folder/'protocol.json',dict(time=api.stamp(),datasets=datasets,seeds=seeds,steps=steps,arms=['gp_iso','gp_uncertainty8','gp_eig8'],batch=32,exploration=8,source='BATCHIE projection: Gaussian conditional information, not original drug tensor model',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()));rows=[]
 for ds in datasets:
  data,oracle=dataset(ds)
  for seed in seeds:
   for arm in ['gp_iso','gp_uncertainty8','gp_eig8']:
    ids=initial(data,seed);curve=[];rec=[];cum=0
    for step in range(steps):
     q=np.array([i for i in range(len(data['genes'])) if i not in set(ids)]);h=reveal(ids,oracle);pred=fit_predict(data,h,q,'gp_iso');sel,diag=acquire(data,h,q,pred,arm)
     assert len(sel)==len(set(sel))==32 and not set(sel)&set(ids)
     np.savez_compressed(folder/f'{ds}_{seed}_{arm}_{step}_predictions.npz',history=ids,query=q,p=pred['p'],effect=pred['effect'],std=pred['std'],selected=sel)
     hit=int(oracle['y'][sel].sum());cum+=hit;curve.append(cum);rec.append(dict(step=step,hits=hit,selected=sel,acquisition=diag));ids+=sel
    row=dict(dataset=ds,seed=seed,method=arm,hits=cum,tests=32*steps,hit_rate=cum/(32*steps),curve=curve,discovery_auc=float(np.mean(curve)),rounds=rec);rows.append(row);write(folder/f'{ds}_{seed}_{arm}.json',row);print(name,ds,seed,arm,cum,flush=True)
 write(folder/'results.json',rows);aggregate(folder,rows)
if __name__=='__main__':run(**json.loads(Path(sys.argv[1]).read_text()))
