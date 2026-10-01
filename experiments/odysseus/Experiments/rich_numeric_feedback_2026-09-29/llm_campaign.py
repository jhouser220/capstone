from llm_study import *

def run(name,datasets,seeds,model=SMALL,steps=3):
 folder=R/name;folder.mkdir(exist_ok=True);arms=['prior','numeric','llm_prior','llm_report','online_blend'];write(folder/'protocol.json',dict(time=api.stamp(),datasets=datasets,seeds=seeds,model=model,steps=steps,arms=arms,batch=32,initial=256,menu='union GPtop48 and geometrytop48; all unselected genes remain in underlying pool',blend='initial numeric weight .5; subsequent softmax(-10*mean previous selected-gene squared errors) versus report LLM; no future labels',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()));rows=[]
 for ds in datasets:
  data,oracle=dataset(ds)
  for seed in seeds:
   for arm in arms:
    ids=initial(data,seed);curve=[];cum=0;rec=[];errn=[];errl=[]
    for step in range(steps):
     h=reveal(ids,oracle);q=np.array([i for i in range(len(data['genes'])) if i not in set(ids)]);pr=fit_predict(data,h,q,'gp_iso');numeric=choose(pr,q);info={};w=.5;lp=None;report=None
     if arm=='prior':sel=sorted(q,key=lambda i:-data['central'][i])[:32];sel=[int(i) for i in sel];menu=q.tolist();scores=None
     elif arm=='numeric':sel=numeric;menu=q.tolist();scores=pr['p']
     else:
      menu=[int(i) for i in dict.fromkeys(choose(pr,q,48)+sorted(q,key=lambda i:-data['central'][i])[:48])];np.random.default_rng(seed+step).shuffle(menu)
      report,pred=report_for(data,h,menu,'gp_iso')
      if arm=='llm_prior':report='No historical observations supplied.'
      lp,info=get_llm(data,menu,report,model,seed+step, 'campaign_prior' if arm=='llm_prior' else 'campaign_report')
      if arm=='online_blend':
       if errn:
        w=float(np.exp(-10*np.mean(errn))/(np.exp(-10*np.mean(errn))+np.exp(-10*np.mean(errl))))
       scores=w*pred['p']+(1-w)*lp
      else:scores=lp
      sel=choose({'p':scores},menu)
     assert len(set(sel))==32 and not set(sel)&set(ids)
     # Persist complete agent-visible report and predictions before scoring the intervention.
     write(folder/f'{ds}_{seed}_{arm}_{step}_decision.json',dict(history_ids=ids,menu=menu,selected=sel,numeric_draft=numeric,scores=None if scores is None else scores.tolist(),report=report,numerical_weight=w,call=info))
     if arm=='online_blend':
      loc=[menu.index(i) for i in sel];errn.extend(((pred['p'][loc]-oracle['y'][sel])**2).tolist());errl.extend(((lp[loc]-oracle['y'][sel])**2).tolist())
     nh=int(oracle['y'][sel].sum());cum+=nh;curve.append(cum);added=list(set(sel)-set(numeric));removed=list(set(numeric)-set(sel));rec.append(dict(step=step,hits=nh,selected=sel,added_hits=int(oracle['y'][added].sum()),removed_hits=int(oracle['y'][removed].sum()),changes=len(added),weight=w,call=info));ids+=sel
    row=dict(dataset=ds,seed=seed,method=arm,model=model,hits=cum,tests=32*steps,hit_rate=cum/(32*steps),curve=curve,discovery_auc=float(np.mean(curve)),rounds=rec);rows.append(row);write(folder/f'{ds}_{seed}_{arm}.json',row);print(name,ds,seed,arm,cum,flush=True)
 write(folder/'results.json',rows);aggregate(folder,rows)
if __name__=='__main__':run(**json.loads(Path(sys.argv[1]).read_text()))
