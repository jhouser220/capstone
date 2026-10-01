from study import *
from sklearn.feature_selection import mutual_info_regression

def run():
 out=[]
 for ds in ['IFNG','IL2']:
  data,oracle=dataset(ds);ids=initial(data,6301);tr=ids[:192];va=ids[192:];h=reveal(tr,oracle);x=data['f'];y=oracle['y'][tr];eff=np.abs(oracle['effect'][tr]);model=XGBClassifier(n_estimators=100,max_depth=3,learning_rate=.05,n_jobs=1,random_state=0,reg_lambda=5,tree_method='hist');model.fit(x[tr],y);p=model.predict_proba(x[va])[:,1];base=float(np.mean((p-oracle['y'][va])**2));rng=np.random.default_rng(6301)
  perm=[]
  for group in np.array_split(np.arange(16),4):
   diffs=[]
   for rep in range(20):
    xx=x[va].copy();xx[:,group]=xx[rng.permutation(len(va))][:,group];pp=model.predict_proba(xx)[:,1];diffs.append(float(np.mean((pp-oracle['y'][va])**2))-base)
   perm.append(dict(components=(group+1).tolist(),brier_increase_mean=float(np.mean(diffs)),sd=float(np.std(diffs))))
  mi=mutual_info_regression(x[tr],eff,random_state=6301);corr=[float(spearmanr(x[tr,k],eff).statistic) for k in range(16)]
  boot=[]
  for rep in range(30):
   ix=rng.choice(len(tr),len(tr),replace=True);boot.append([float(spearmanr(x[np.array(tr)[ix],k],eff[ix]).statistic) for k in range(16)])
  out.append(dict(dataset=ds,training_ids=tr,validation_ids=va,validation_brier=base,group_permutation=perm,mutual_information=mi.tolist(),spearman=corr,spearman_bootstrap_10_90=np.quantile(boot,[.1,.9],axis=0).tolist(),interpretation='Observed-history diagnostic only; coordinate associations and model attributions are not causal mechanisms. Neither these diagnostics nor held-out candidate scores tune completed comparisons.'))
 write(R/'OBSERVED_HISTORY_TOOL_DIAGNOSTICS.json',out)
if __name__=='__main__':run()
