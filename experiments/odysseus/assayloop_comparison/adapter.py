"""Adapters for the authors' AssayBench SequentialLoop and AssayLoop classes.

Inputs are explicit source model, source state callback, bank and acquisition.
The module itself performs no inference requests, loading, or credential access.
"""
import numpy as np
from .numerics import map_state_precision, posterior_scores, rank_fusion

ARMS=('native','native_filled','incumbent','ordinary_full','assayformer_only','historical_only')

class SourceState:
    """Frozen 10-D AssayFormer state using a supplied cached description vector."""
    def __init__(self, source, description):
        self.source=source;self.description=np.asarray(description,dtype=np.float32)
        if getattr(source,'encoder_type','')=='modernbert_text':
            raise ValueError("This comparison requires the frozen numeric-context encoder")
        self.W=source.net.gene_emb.weight.detach().cpu().numpy().copy()
        self.b=source.net.gene_bias.detach().cpu().numpy().copy()
        if self.W.shape[1]!=10:raise ValueError("Expected the frozen 10-D source")
    def ids(self, genes):
        out=np.asarray([self.source.vocab.to_idx(str(g)) for g in genes],dtype=int)
        if len(set(out.tolist()))!=len(out):raise ValueError("Gene identities collide in source vocabulary")
        if len(out) and np.any(out==0):raise ValueError("Candidate absent from source vocabulary")
        return out
    def __call__(self, observations):
        import torch
        net=self.source.net;dev=self.source.device
        with torch.no_grad():
            ids=torch.tensor(self.ids([o.candidate for o in observations])[None],device=dev)
            hits=torch.tensor([[int(bool(o.label.get('hit'))) for o in observations]],dtype=torch.long,device=dev)
            desc=torch.tensor(self.description[None],device=dev)
            pad=torch.zeros_like(ids,dtype=torch.bool)
            return net.out_proj(net._encode_hidden(desc,ids,hits,pad))[0].cpu().numpy()

class Recipient:
    def __init__(self, arm, source, state, bank, screen_id, draw):
        self.arm=arm;self.source=source;self.state=state;self.bank=bank
        self.screen_id=str(screen_id);self.draw=draw;self.audit=[];self.selector=None
        if arm in ('incumbent','historical_only'):
            if bank is None:raise ValueError("Historical bank required")
            self.bank_index={g:i for i,g in enumerate(bank.genes)}
            # Bank names are checked against source IDs at every used identity.
    def name(self):return 'assayloop_comparison_'+self.arm
    def reset(self):self.audit=[];self.selector=None;self.source.reset()
    def predict(self, observations, candidates, task_context=None):
        from assaybench.core.types import ModelPrediction
        step=len(self.audit)+1
        if self.arm in ('native','native_filled','assayformer_only'):
            result=self.source.predict(observations,candidates,task_context)
            self.audit.append({'step':step,'recipient':'source'})
            return result
        z=self.state(observations);ids=self.state.ids(candidates)
        selected=[o for o in observations if o.metadata.get('in_library',False)]
        observed=self.state.ids([o.candidate for o in selected])
        labels=np.asarray([float(bool(o.label['hit'])) for o in selected])
        W,b=self.state.W,self.state.b
        mode,_=map_state_precision(z,np.eye(10),W[observed],b[observed],labels)
        # Proposal ranking preserves the E159 float32 source dot product.
        proposal_base=W[ids]@z+b[ids];post=W[ids]@mode+b[ids]
        self.selector=(list(candidates),proposal_base,post)
        if self.arm=='ordinary_full':
            # Numerical source score uses the source's original torch scoring.
            source_prediction=self.source.predict(observations,candidates,task_context)
            scores=rank_fusion([source_prediction.scores[g] for g in candidates],post)
            detail={'mode':mode.tolist(),'measured_labels':len(labels)}
        else:
            all_genes=list(candidates)+[o.candidate for o in selected]
            if any(str(g) not in self.bank_index for g in all_genes):raise ValueError("Bank vocabulary mismatch")
            ci=np.asarray([self.bank_index[str(g)] for g in candidates],int)
            oi=np.asarray([self.bank_index[str(o.candidate)] for o in selected],int)
            design=np.concatenate((W[ids].astype(float),.25*self.bank.residual[ci]),axis=1)
            obs_design=np.concatenate((W[observed].astype(float),.25*self.bank.residual[oi]),axis=1)
            prior=np.r_[z.astype(float),self.bank.mean]
            precision=np.diag(np.r_[np.ones(10),np.full(32,1/self.bank.scalar_c)])
            scores,detail=posterior_scores(prior,precision,obs_design,b[observed],labels,
                design,b[ids],self.bank.availability[ci],screen_id=self.screen_id,draw=self.draw,step=step)
        self.audit.append({'step':step,'recipient':self.arm,**detail})
        return ModelPrediction(scores=dict(zip(candidates,map(float,scores))))

class Acquisition:
    def __init__(self, arm, recipient, actor):
        from assayloop.acquisitions.greedy_from_model import GreedyFromModel
        self.arm=arm;self.recipient=recipient;self.actor=actor
        self.greedy=GreedyFromModel(seed=0,epsilon=0)
        self.fill=GreedyFromModel(seed=0,epsilon=0)
        self.audit=[]
    def name(self):return self.arm
    def reset(self):
        self.greedy.reset();self.fill.reset();self.audit=[]
        if self.actor is not None:self.actor.reset()
    def last_trace(self):return self.audit[-1] if self.audit else {}
    def suggest(self, **kw):
        if kw['batch_size']!=100:raise ValueError("Frozen comparison requires100 slots per decision")
        step=len(kw['history'])+1;eligible=list(kw['candidates'])
        if len(set(eligible))!=len(eligible):raise ValueError("Duplicate eligible candidates")
        if self.arm in ('assayformer_only','historical_only') or step>3:
            picked=self.greedy.suggest(**kw);proposal=[]
        else:
            proposal=list(self.actor.suggest(**kw))
            if len(proposal)>100 or len(set(proposal))!=len(proposal) or not set(proposal)<=set(eligible):
                raise ValueError("Actor returned invalid/duplicate/unbound candidates")
            if self.arm=='native':picked=proposal
            elif self.arm=='native_filled':
                used=set(proposal);remaining=[g for g in eligible if g not in used]
                fill=self.fill.suggest(**{**kw,'candidates':remaining,'batch_size':100-len(proposal)}) if len(proposal)<100 else []
                picked=proposal+fill
            else:
                candidates,base,post=self.recipient.selector
                if candidates!=eligible:raise ValueError("Selector state is stale")
                pos={g:i for i,g in enumerate(candidates)};pi=np.asarray([pos[g] for g in proposal],int)
                ranked=list(proposal)
                if len(proposal)>25:
                    values=rank_fusion(base[pi],post[pi])
                    ranked=[proposal[i] for i in np.argsort(-values,kind='stable')]
                numeric=self.greedy.suggest(**kw);kept=ranked[:25];used=set(kept)
                picked=(kept+[g for g in numeric if g not in used])[:100]
        if len(set(picked))!=len(picked) or not set(picked)<=set(eligible):raise ValueError("Invalid selected batch")
        if self.arm!='native' and len(picked)!=min(100,len(eligible)):raise ValueError("Budget completion failed")
        self.audit.append({'step':step,'arm':self.arm,'requested':100,'selected':list(picked),
                           'proposal_count':len(proposal)})
        return picked

def run_campaign(*, arm, task, source_model, screen_id, draw=0, source_state=None,
                 bank=None, actor_acquisition=None, allow_live_calls=False):
    """Execute the full authors' loop with explicit adapters and fresh task/actor.

    Caller owns provider budget, request/response capture and dataset admission.
    This opt-in is additional to that caller's normal funding safeguards.
    """
    from assaybench.core.loop import SequentialLoop
    if list(task.initial_observations()):raise ValueError('Frozen comparison requires an empty initial history')
    if arm not in ARMS:raise ValueError(f"Unknown arm: {arm}")
    if draw not in (0,1):raise ValueError("Only repetitions0/1 are pinned")
    if arm not in ('assayformer_only','historical_only') and (actor_acquisition is None or not allow_live_calls):
        raise ValueError("Provide a genuine actor and explicitly enable its calls")
    if arm in ('incumbent','ordinary_full','historical_only') and source_state is None:
        raise ValueError("Explicit source-state adapter required")
    model=Recipient(arm,source_model,source_state,bank,screen_id,draw)
    acquisition=Acquisition(arm,model,actor_acquisition)
    result=SequentialLoop(task,model,acquisition,[],batch_size=100,max_shortfall_frac=1.).run(10,verbose=False)
    return result,{'model':model.audit,'acquisition':acquisition.audit}
