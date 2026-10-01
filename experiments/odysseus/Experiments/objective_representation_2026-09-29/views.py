from study2 import *
from router import validate
import paid_api
from scipy.stats import rankdata

def build(ds,seed):
 data,oracle=base.dataset(ds);ini=base.initial(data,seed);h=base.reveal(ini,oracle);q=np.array([i for i in range(len(data['genes'])) if i not in set(ini)]);gp=predict(data,h,q,'gp_iso');tree=predict(data,h,q,'full:effect_median');menu=[int(i) for i in dict.fromkeys(base.choose(gp,q,48)+base.choose(tree,q,48)+sorted(q,key=lambda i:-data['central'][i])[:16])];np.random.default_rng(seed).shuffle(menu);loc={g:i for i,g in enumerate(q)};rg=rankdata(gp['p'])/len(q);rt=rankdata(tree['p'])/len(q);sim=data['x'][menu]@data['x'][ini].T;cards=[]
 for j,g in enumerate(menu):
  i=loc[g];card=dict(id=j,gp_probability=round(float(gp['p'][i]),5),gp_mean_absolute_effect=round(float(gp['effect'][i]),5),tree_median_absolute_effect=round(float(tree['p'][i]),5),gp_rank_percentile=round(float(rg[i]),5),tree_rank_percentile=round(float(rt[i]),5));examples=[]
  for hit in [1,0]:
   candidates=np.flatnonzero(h['y']==hit)
   if len(candidates):
    k=int(candidates[np.argmax(sim[j,candidates])]);examples.append(dict(gene=str(data['genes'][ini[k]]),hit=int(h['y'][k]),effect=round(float(h['effect'][k]),5),similarity=round(float(sim[j,k]),4)))
  card['examples']=examples;cards.append(card)
 val=validate(data,oracle,ini);reliability={m:{k:v for k,v in met.items() if k!='scores'} for m,met in val['metrics'].items()}
 return dict(dataset=ds,seed=seed,task=data['task'],initial=ini,menu=menu,genes=[str(data['genes'][g]) for g in menu],cards=cards,reliability=reliability,gp_selected=base.choose(gp,q),tree_selected=base.choose(tree,q)),oracle

def run(name):
 folder=R/name;folder.mkdir(exist_ok=True);write(folder/'protocol.json',dict(time=api.stamp(),datasets=['IFNG','IL2'],lite_seeds=[7501,7502],astra_seeds=[7501],arms=['none','scores','examples','shuffled_scores(lite only)'],batch=32,menu='GPtop48+fullmedianTop48+geometrytop16, shuffled, shared',sources='MESA/JitMem view selection projection; LLMNN genuine/no-report/shuffled controls; LGBO numerical/LLM inspiration',new_batch_cap=2.5));rows=[]
 for ds in ['IFNG','IL2']:
  for seed in [7501,7502]:
   s,oracle=build(ds,seed);write(folder/f'{ds}_{seed}_visible.json',s)
   for model in (['google/gemini-2.5-flash-lite','openai/gpt-6-astra'] if seed==7501 else ['google/gemini-2.5-flash-lite']):
    for arm in (['none','scores','examples','shuffled_scores'] if model.startswith('google') else ['none','scores','examples']):
     report='No measured history or numerical advice supplied.'
     if arm!='none':
      cards=[dict(c) for c in s['cards']]
      if arm!='examples':cards=[{k:v for k,v in c.items() if k!='examples'} for c in cards]
      if arm=='shuffled_scores':
       perm=np.random.default_rng(seed+19).permutation(len(cards));cards=[dict(cards[int(perm[i])],id=i) for i in range(len(cards))]
      report=dict(validation='192 train /64 already-observed validation; top16_hits higher is better; MAE lower is better; small biased sample may not generalize',validation_metrics=s['reliability'],warning='GP probabilities are uncalibrated; tree magnitudes and rank percentiles are NOT hit probabilities. Neighbour similarity is not biological causation. Positive and negative signed effects can both be hits.',candidates=cards)
     prompt='Task:'+json.dumps(s['task'])+'\nSelect exactly32 candidate IDs to maximize the number of hits (large measured effects in either direction). Outcomes for these candidates are unknown. Use your biological prior and evidence if supplied.\nCandidates:'+json.dumps([dict(id=i,gene=g) for i,g in enumerate(s['genes'])],separators=(',',':'))+'\nReport:'+json.dumps(report,separators=(',',':'))+'\nReturn JSON {"selected_ids":[32 distinct integer candidate IDs],"preferred_model":"gp/tree/neither","reason":"one brief evidence-based sentence"}. Select only IDs in the menu; do not select historical examples.'
     cid=f'{ds}_{seed}_{arm}_'+paid_api.digest([model,prompt])[:16];response=paid_api.request_llm(cid,prompt,model,seed,4000 if model.startswith('openai') else 1000);obj=paid_api.parse_result(response);ids=obj.get('selected_ids',[]) if isinstance(obj,dict) else [];valid=len(ids)==32 and len(set(ids))==32 and all(type(i)==int and 0<=i<len(s['menu']) for i in ids)
     write(folder/f'{ds}_{seed}_{model.split("/")[-1]}_{arm}_decision.json',dict(raw_output=obj,valid=valid,call_id='obj_'+cid,selected_ids=ids,prompt_sha256=paid_api.digest(prompt)))
     if not valid:
      original=obj;repair_prompt=prompt+'\nFORMAT CORRECTION: your previous response was '+json.dumps(obj)+'. Return exactly 32 distinct candidate IDs. The validation top16 metric is NOT the requested batch size. This is a formatting correction only; no new experimental outcomes are provided.'
      response=paid_api.request_llm(cid+'_count_repair',repair_prompt,model,seed+100000,4000 if model.startswith('openai') else 1000);obj=paid_api.parse_result(response);ids=obj.get('selected_ids',[]) if isinstance(obj,dict) else [];valid=len(ids)==32 and len(set(ids))==32 and all(type(i)==int and 0<=i<len(s['menu']) for i in ids)
      write(folder/f'{ds}_{seed}_{model.split("/")[-1]}_{arm}_repair.json',dict(original=original,repaired=obj,valid=valid,call_id='obj_'+cid+'_count_repair'))
      if not valid:raise RuntimeError('One repair failed; stop without imputing genes; '+cid)
     repaired=(folder/f'{ds}_{seed}_{model.split("/")[-1]}_{arm}_repair.json').exists()
     selected=[s['menu'][i] for i in ids];draft=set(s['gp_selected']);added=list(set(selected)-draft);removed=list(draft-set(selected));r=dict(dataset=ds,seed=seed,model=model,arm=arm,hits=int(oracle['y'][selected].sum()),tests=32,hit_rate=float(oracle['y'][selected].mean()),added_hits=int(oracle['y'][added].sum()),removed_hits=int(oracle['y'][removed].sum()),changes=len(added),output_repaired=repaired,preferred_model=obj.get('preferred_model'),selected=selected);rows.append(r);write(folder/'results_so_far.json',rows);print(ds,seed,model,arm,r['hits'],flush=True)
   # Numerical baselines saved before outcomes and evaluated once per state.
   for arm in ['gp','tree']:
    sel=s[arm+'_selected'];rows.append(dict(dataset=ds,seed=seed,model='numerical',arm=arm,hits=int(oracle['y'][sel].sum()),tests=32,selected=sel))
 write(folder/'results.json',rows);summary={}
 for model in sorted({r['model'] for r in rows}):
  summary[model]={}
  for arm in sorted({r['arm'] for r in rows if r['model']==model}):
   rr=[r for r in rows if r['model']==model and r['arm']==arm];summary[model][arm]=dict(hits=sum(r['hits'] for r in rr),tests=sum(r['tests'] for r in rr),cases=[r['hits'] for r in rr])
 write(folder/'summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':run('r5_feedback_views')
