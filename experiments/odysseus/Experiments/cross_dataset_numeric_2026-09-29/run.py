"""Fresh cross-screen campaigns. No paid APIs. Targets are evaluator-only except revealed history."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,csv,json,time,hashlib,warnings
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'objective_representation_2026-09-29'))
import study2 as old
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
write=old.write
DATA=old.base.DATA
FEATURE=old.base.FEATURE
CACHE={}
def dataset(ds):
 if ds in CACHE:return CACHE[ds]
 rows=list(csv.DictReader((DATA/f'ground_truth_{ds}.csv').open()));cols=('0','1') if ds=='Steinhart_crispra_GD2_D22' else ('Gene','Score')
 from collections import Counter
 counts=Counter(r[cols[0]] for r in rows);duplicates={g for g,n in counts.items() if n>1}
 assert not duplicates or ds=='Scharenberg22'
 scores={r[cols[0]]:float(r[cols[1]]) for r in rows if r[cols[0]] not in duplicates}
 z=np.load(FEATURE,allow_pickle=False);keep=[i for i,g in enumerate(z['genes']) if g in scores];genes=z['genes'][keep];x=z['values'][keep].astype(float);assert len(genes)>352;assert np.isfinite(x).all()
 cache=ROOT/(ds+'_features.npz')
 if cache.exists():q=np.load(cache);assert np.array_equal(q['genes'],genes);f=q['f'];central=q['central']
 else:
  p=PCA(n_components=16,random_state=4201,svd_solver='randomized');f=p.fit_transform(x);f/=np.maximum(f.std(0),1e-12);central=x@x.mean(0);np.savez_compressed(cache,genes=genes,f=f,central=central)
 hits=set(map(str,np.load(DATA/f'topmovers_{ds}.npy',allow_pickle=True)));o=dict(y=np.array([int(g in hits) for g in genes]),effect=np.array([scores[g] for g in genes]));assert np.isfinite(o['effect']).all();data=dict(genes=genes,x=x,f=f,central=central)
 CACHE[ds]=(data,o);return data,o

def run(cfg):
 folder=ROOT/cfg['name'];folder.mkdir(exist_ok=True)
 if (folder/'protocol.json').exists():
  previous=json.loads((folder/'protocol.json').read_text());(folder/('protocol_'+previous['code_sha256'][:10]+'.json')).write_text(json.dumps(previous,indent=2));(folder/('code_'+previous['code_sha256'][:10]+'.py')).write_text((folder/'code_snapshot.py').read_text())
 n_initial=cfg.get('initial',256)
 write(folder/'protocol.json',dict(cfg,time=old.api.stamp(),initial=n_initial,batch=32,steps=3,feature_sha256=hashlib.sha256(FEATURE.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),data_hashes={str(p.relative_to(DATA)):hashlib.sha256(p.read_bytes()).hexdigest() for ds in cfg['datasets'] for p in [DATA/f'ground_truth_{ds}.csv',DATA/f'topmovers_{ds}.npy']},api_cost=0,duplicate_policy='Exclude all ambiguous duplicate gene IDs; six in Scharenberg; same pool in all arms'));(folder/'code_snapshot.py').write_text(Path(__file__).read_text());rows=[]
 for ds in cfg['datasets']:
  data,o=dataset(ds)
  for seed in cfg['seeds']:
   assert cfg['initial_mode']!='mixed' or n_initial==256
   ini=old.base.initial(data,seed) if cfg['initial_mode']=='mixed' else np.random.default_rng(seed).choice(len(data['genes']),n_initial,replace=False).tolist()
   for method in cfg['methods']:
    path=folder/f'{ds}_{seed}_{method.replace(":","-")}.json'
    if path.exists():rows.append(json.loads(path.read_text()));continue
    start=time.monotonic();ids=list(ini);curve=[];steps=[];cum=0
    for step in range(3):
     q=np.array([i for i in range(len(data['genes'])) if i not in set(ids)]);hist=ini if cfg.get('frozen',False) else ids;h=old.base.reveal(hist,o)
     if method=='random':pred=dict(p=np.random.default_rng(seed+10000+step).random(len(q)))
     elif method=='full:hit_logistic' and len(set(h['y']))==1:pred=dict(p=np.repeat(float(h['y'].mean()),len(q)))
     elif method=='gp_mean':
      pred=old.predict(data,h,q,'gp_iso');pred['p']=pred['effect'].copy()
     else:pred=old.predict(data,h,q,method)
     assert np.isfinite(pred['p']).all();sel=old.base.choose(pred,q);assert len(sel)==len(set(sel))==32 and not set(sel)&set(ids)
     np.savez_compressed(folder/f'{ds}_{seed}_{method.replace(":","-")}_{step}.npz',query=q,score=pred['p'],history=hist,selected=sel)
     hits=int(o['y'][sel].sum());cum+=hits;curve.append(cum)
     record=dict(step=step,selected=sel,hits=hits,remaining_hits=int(o['y'][q].sum()),remaining_n=len(q),initial_hits=int(o['y'][ini].sum()),spearman_abs_effect=float(spearmanr(pred['p'],np.abs(o['effect'][q])).statistic))
     if method in ['gp_iso','full:hit_logistic']:
      pp=np.clip(pred['p'],1e-6,1-1e-6);record['brier']=float(np.mean((pp-o['y'][q])**2))
     steps.append(record);ids+=sel
    r=dict(dataset=ds,seed=seed,method=method,initial_mode=cfg['initial_mode'],frozen=cfg.get('frozen',False),initial=ini,hits=cum,tests=96,hit_rate=cum/96,curve=curve,discovery_auc=float(np.mean(curve)),initial_hits=int(o['y'][ini].sum()),remaining_initial_hits=int(o['y'].sum()-o['y'][ini].sum()),incremental_recall=cum/max(1,int(o['y'].sum()-o['y'][ini].sum())),rounds=steps,runtime_seconds=time.monotonic()-start);write(path,r);rows.append(r);print(cfg['name'],ds,seed,method,cum,flush=True)
    write(folder/'results_so_far.json',rows)
 write(folder/'results.json',rows);print('DONE',cfg['name'],len(rows),flush=True)

if __name__=='__main__':
 for arg in sys.argv[1:]:run(json.loads(Path(arg).read_text()))
