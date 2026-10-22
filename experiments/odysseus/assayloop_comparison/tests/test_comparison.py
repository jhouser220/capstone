import csv,io,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from assayloop_comparison.numerics import map_state_precision,posterior_scores,rank_fusion
from assayloop_comparison.assets import HistoricalBank,build_bank
from assayloop_comparison.audit_results import recount
D=Path(__file__).resolve().parents[1]

class NumericalTests(unittest.TestCase):
    def test_exact_frozen_recipient_fixture(self):
        a=json.loads((D/'tests/numerical_reference.json').read_text())
        W=np.asarray(a['W']);V=np.asarray(a['V']);b=np.asarray(a['bias']);oi=a['observed_indices'];ci=a['candidate_indices']
        X=np.c_[W,.25*V];prior=np.r_[a['z'],a['u0']];precision=np.diag(np.r_[np.ones(10),np.full(32,.5)])
        score,detail=posterior_scores(prior,precision,X[oi],b[oi],np.asarray(a['labels']),X[ci],b[ci],np.asarray(a['q'])[ci],screen_id='fixture',draw=0,step=1)
        np.testing.assert_allclose(score,a['expected_scores'],rtol=0,atol=1e-12)
        np.testing.assert_allclose(detail['mode'],a['expected_mode'],rtol=0,atol=1e-12)
    def test_empty_history_returns_prior(self):
        z=np.array([.2,-.4]);P=np.diag([2.,3.])
        mu,L=map_state_precision(z,P,np.empty((0,2)),np.empty(0),np.empty(0))
        np.testing.assert_array_equal(mu,z);np.testing.assert_allclose(L@L.T,P)
    def test_exact_ties_have_equal_fusion_scores(self):
        np.testing.assert_array_equal(rank_fusion(np.array([1.,1.,0.]),np.array([2.,2.,0.])),[-1.5,-1.5,-3])
    def test_bank_rejects_query_donor(self):
        with self.assertRaises(ValueError):HistoricalBank(['g'],np.zeros((1,32)),np.zeros(32),np.ones(1),1,'p',('p',))
    def test_query_outcomes_do_not_affect_constructed_bank(self):
        rng=np.random.default_rng(4);M=(rng.random((40,45))>.2).astype(float);Y=(rng.random((40,45))>.7).astype(float)
        kwargs=dict(genes=[str(i) for i in range(45)],source_weights=rng.normal(size=(45,10)),donor_mask=M,donor_hits=Y,donor_publications=['query']+['donor']*39,donor_descriptions=rng.normal(size=(40,6)),query_description=np.ones(6),query_publication='query')
        a=build_bank(**kwargs);Y2=Y.copy();Y2[0]=1-Y2[0];kwargs['donor_hits']=Y2;b=build_bank(**kwargs)
        np.testing.assert_allclose(a.residual,b.residual,atol=1e-12)
        np.testing.assert_array_equal(a.mean,b.mean);self.assertEqual(a.donor_publications,('donor',))
    def test_historical_results_and_paired_budgets(self):
        r=recount();self.assertEqual(r['campaigns'],342);self.assertEqual((r['publication_wins'],r['publication_losses']),(11,7))
        expected={'original_test3':11.75,'additional_test6':26.33333333333334,'validation6':15.20833333333334,'training3':16.416666666666664}
        for k,v in expected.items():self.assertAlmostEqual(r['incumbent_minus'][k]['native_filled'],v)
    def test_duplicate_records_fail(self):
        text=(D/'recorded_results/campaign_metrics.csv').read_text().splitlines()
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bad.csv';p.write_text('\n'.join(text+[text[1]])+'\n')
            with self.assertRaises(ValueError):recount(p)

try:
    from assaybench.core.types import Observation,ModelPrediction
    from assayloop.acquisitions.greedy_from_model import GreedyFromModel
    UPSTREAM=True
except ImportError: UPSTREAM=False

@unittest.skipUnless(UPSTREAM,'Install pinned upstream AssayLoop/AssayBench for integration checks')
class LoopTests(unittest.TestCase):
    def make(self):
        from assayloop_comparison.adapter import run_campaign
        class Task:
            def __init__(self):self.seen={}
            def task_id(self):return 'synthetic_fixture'
            def context(self):return {}
            def initial_observations(self):return []
            def candidates(self):return [str(i) for i in range(1100) if str(i) not in self.seen]
            def reveal(self,genes):
                assert len(set(genes))==len(genes) and not set(genes)&set(self.seen)
                rows=[Observation(g,{'hit':int(g)%7==0},{'in_library':int(g)%5!=0}) for g in genes]
                self.seen.update({o.candidate:o for o in rows});return rows
        class Source:
            def reset(self):pass
            def predict(self,observations,candidates,task_context=None):return ModelPrediction(scores={g:-float(g) for g in candidates})
        class Actor:
            def reset(self):self.calls=0
            def suggest(self,**kw):
                self.calls+=1
                assert len(kw['history'])==self.calls-1
                return kw['candidates'][-40:]
        class State:
            W=np.random.default_rng(8).normal(size=(1100,10))*.01;b=np.zeros(1100)
            def ids(self,genes):return np.asarray([int(g) for g in genes],int)
            def __call__(self,obs):return np.zeros(10)
        bank=HistoricalBank([str(i) for i in range(1100)],np.random.default_rng(9).normal(size=(1100,32))*.01,np.zeros(32),np.full(1100,.8),1,'query',('donor',))
        return run_campaign,Task,Source,Actor,State,bank
    def test_full_ten_decision_arms(self):
        run,Task,Source,Actor,State,bank=self.make()
        for arm in ('native','native_filled','incumbent','ordinary_full','assayformer_only','historical_only'):
            with self.subTest(arm=arm):
                task=Task();actor=Actor()
                result,audit=run(arm=arm,task=task,source_model=Source(),source_state=State(),bank=bank,actor_acquisition=actor,screen_id='synthetic',allow_live_calls=True)
                self.assertEqual(len(result.history),10)
                expected=820 if arm=='native' else 1000
                self.assertEqual(len(task.seen),expected)
                self.assertEqual(actor.calls,0 if arm in ('assayformer_only','historical_only') else 3)
                self.assertTrue(all(r.metrics.get('batch_size',0)==len(r.acquired_batch) for r in result.history))
    def test_transport_setup_is_idempotent(self):
        import openai
        from assayloop_comparison.reproduce import install_transport
        install_transport(); first=openai.OpenAI
        install_transport(); self.assertIs(openai.OpenAI,first)

    def test_live_execution_requires_explicit_opt_in(self):
        run,Task,Source,Actor,State,bank=self.make()
        with self.assertRaises(ValueError):run(arm='native',task=Task(),source_model=Source(),actor_acquisition=Actor(),screen_id='synthetic')
    def test_native_fill_preserves_all_proposals(self):
        run,Task,Source,Actor,State,bank=self.make();task=Task();proposed=task.candidates()[-40:]
        result,_=run(arm='native_filled',task=task,source_model=Source(),actor_acquisition=Actor(),screen_id='synthetic',allow_live_calls=True)
        self.assertEqual(result.history[0].acquired_batch[:40],proposed)

if __name__=='__main__':unittest.main()
