"""Recount compact external evidence without dependencies, network or credentials."""
import collections,json,statistics
from pathlib import Path
D=Path(__file__).resolve().parent/'recorded_results'
def recount(path=D/'external_campaigns.json'):
 rows=json.loads(Path(path).read_text());seen=set();groups=collections.defaultdict(list)
 for r in rows:
  key=(r['screen'],r['method'],r['actor'],r['draw'])
  if key in seen:raise ValueError('Duplicate campaign')
  seen.add(key)
  if len(r['round_hits'])!=10 or any(not isinstance(n,int) or n<0 or n>100 for n in r['round_hits']):raise ValueError('Invalid round hit counts')
  if sum(r['round_hits'])!=r['hits'] or r['accepted']!=1000 or r['measured']+r['unknown']!=1000 or r['hits']>r['measured']:raise ValueError('Invalid hit/budget accounting')
  cumulative=[0]
  for x in r['round_hits']:cumulative.append(cumulative[-1]+x)
  auc=sum((a+b)/2 for a,b in zip(cumulative,cumulative[1:]))/10
  if abs(auc-r['requested_auc'])>1e-8:raise ValueError('Curve area mismatch')
  groups[r['source'],r['screen'],r['method']].append(r)
 screens=[]
 for src,sid in sorted({(r['source'],r['screen']) for r in rows}):
  a,b=groups[src,sid,'Native'],groups[src,sid,'HSA']
  if len(a)!=4 or len(b)!=4 or {(x['actor'],x['draw']) for x in a}!={(x['actor'],x['draw']) for x in b}:raise ValueError('Incomplete or unmatched four-path assay')
  screens.append(dict(source=src,screen=sid,original_assayloop=statistics.mean(x['hits'] for x in a),assayloop_hsa=statistics.mean(x['hits'] for x in b)))
 pubs=[]
 for src in sorted({r['source'] for r in screens}):
  group=[r for r in screens if r['source']==src];n=statistics.mean(r['original_assayloop'] for r in group);h=statistics.mean(r['assayloop_hsa'] for r in group)
  pubs.append(dict(source=src,assays=len(group),original_assayloop=n,assayloop_hsa=h,delta=h-n,relative_gain_percent=100*(h-n)/n))
 return dict(campaigns=len(rows),assays=len(screens),publications=pubs)
if __name__=='__main__':print(json.dumps(recount(),indent=2))
