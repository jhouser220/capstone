"""Offline recount from public campaign-level historical records; no model calls."""
import csv,json,math
from collections import defaultdict
from pathlib import Path

def recount(path=None):
    path=Path(path) if path else Path(__file__).parent/'recorded_results/campaign_metrics.csv'
    with path.open() as stream: rows=list(csv.DictReader(stream))
    seen=set();groups=defaultdict(list)
    for r in rows:
        key=tuple(r[k] for k in ('screen_id','method','actor','draw'))
        if key in seen:raise ValueError(f'Duplicate campaign: {key}')
        seen.add(key)
        if int(r['decision'])!=10 or int(r['requested'])!=1000:raise ValueError('Incomplete campaign')
        if r['method'] in ('incumbent','native_filled','ordinary_full') and int(r['selected'])!=1000:
            raise ValueError('Unequal selected budget')
        if not 0<=int(r['hits'])<=int(r['measured'])<=int(r['selected'])<=1000:raise ValueError('Invalid counts')
        groups[(r['cohort'],r['method'],r['publication_id'])].append(r)
    means={};effects={};wins=losses=0
    for cohort in sorted({r['cohort'] for r in rows}):
        means[cohort]={}
        for method in sorted({r['method'] for r in rows}):
            pubs=[sum(float(r['hits']) for r in rr)/len(rr) for (c,m,p),rr in groups.items() if c==cohort and m==method]
            if pubs:means[cohort][method]=sum(pubs)/len(pubs)
        effects[cohort]={m:means[cohort]['incumbent']-means[cohort][m]
                         for m in ('original_native','native_filled','ordinary_full','AssayFormer_only')}
        pubs={r['publication_id'] for r in rows if r['cohort']==cohort}
        for pub in pubs:
            a=groups[cohort,'incumbent',pub];b=groups[cohort,'native_filled',pub]
            cells=lambda rr:{(r['screen_id'],r['actor'],r['draw']) for r in rr}
            if cells(a)!=cells(b) or len(a)!=4:raise ValueError('Missing matched actor/repetition cell')
            delta=sum(int(r['hits']) for r in a)/4-sum(int(r['hits']) for r in b)/4
            wins+=delta>0;losses+=delta<0
    return {'campaigns':len(rows),'means':means,'incumbent_minus':effects,
            'publication_wins':wins,'publication_losses':losses}

if __name__=='__main__':print(json.dumps(recount(),indent=2))
