import json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from assayloop_comparison.audit_external import recount,D
class ExternalResultsTests(unittest.TestCase):
 def test_complete_paired_panel(self):
  r=recount();self.assertEqual((r['campaigns'],r['assays'],len(r['publications'])),(88,11,7))
  expected={'Schukken':(47.75,50),'Abt–Wang':(322,365.5),'Wilson':(97.25,107.75),'Nguyen':(13.5,15),'Shan':(76.4166666667,86.8333333333),'Herranz-Ors':(657.0833333333,712),'Samelson':(113,155)}
  for p in r['publications']:
   a,b=expected[p['source']];self.assertAlmostEqual(p['original_assayloop'],a);self.assertAlmostEqual(p['assayloop_hsa'],b)
 def test_missing_or_duplicate_campaign_rejected(self):
  rows=json.loads((D/'external_campaigns.json').read_text())
  for value in (rows[:-1],rows+[rows[0]]):
   with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)/'rows.json';p.write_text(json.dumps(value))
    with self.assertRaises(ValueError):recount(p)
 def test_reward_mutation_rejected(self):
  rows=json.loads((D/'external_campaigns.json').read_text());rows[0]['hits']+=1
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'rows.json';p.write_text(json.dumps(rows))
   with self.assertRaises(ValueError):recount(p)
try:
 from assayloop_comparison.feedback import MeasuredOnlyAcquisition,LLMSingleAcquisition
 AVAILABLE=True
except ImportError:AVAILABLE=False
@unittest.skipUnless(AVAILABLE,'Pinned upstream required')
class FeedbackTests(unittest.TestCase):
 def test_unknown_outcomes_removed_without_mutating_numerical_history(self):
  measured=SimpleNamespace(candidate='measured',metadata={'in_library':True},label={'hit':True})
  missing=SimpleNamespace(candidate='missing',metadata={'in_library':False},label={'hit':False})
  history=[SimpleNamespace(step=1,new_observations=[measured,missing])]
  actor=object.__new__(MeasuredOnlyAcquisition)
  def check_parent(self,rows,batch,context):
   selftest.assertEqual(rows[0].new_observations,[measured]);selftest.assertNotIn('total_hits',context)
   return 'system','history of measured',1
  selftest=self
  with patch.object(LLMSingleAcquisition,'_build_prompt',check_parent):
   _,prompt,n=actor._build_prompt(history,100,{'total_hits':99,'name':'screen'})
  self.assertEqual(history[0].new_observations,[measured,missing]);self.assertIn('Round 1: missing',prompt);self.assertIn('no outcome is known',prompt)
if __name__=='__main__':unittest.main()
