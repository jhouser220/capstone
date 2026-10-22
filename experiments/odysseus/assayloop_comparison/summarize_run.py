"""Publication-balanced summaries and optional historical CSV equivalence."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def summarize(summary,reference=None):
    rows=summary['rows']
    if summary['failed'] or len(rows)!=summary['expected']:
        raise ValueError('A complete successful matrix is required')
    groups=defaultdict(list);seen=set()
    names={'native':'original_native','source':'AssayFormer_only','scalar':'cold_incumbent'}
    ref={}
    if reference:
        for r in csv.DictReader(Path(reference).open()):
            ref[(r['screen_id'],r['method'],r['actor'],int(r['draw']))]=r
    for r in rows:
        c=r['cell'];key=(c['screen_id'],names.get(c['arm'],c['arm']),c['actor'],c['draw'])
        if key in seen:raise ValueError('Duplicate cell')
        seen.add(key)
        if reference:
            old=ref[key];m=r['metrics']
            for new,field in [('unique_hits','hits'),('selected','selected'),('measured','measured'),('unavailable','unavailable')]:
                if m[new]!=int(old[field]):raise ValueError(f'Historical difference: {key} {field}')
            if abs(m['ef']-float(old['official_ef']))>1e-10:raise ValueError('Historical EF difference')
        groups[(r['cohort'],c['arm'],r['publication_id'])].append(r)
    output={}
    for cohort in sorted({r['cohort'] for r in rows}):
        methods={}
        for arm in sorted({r['cell']['arm'] for r in rows}):
            pubs=[v for (co,a,p),v in groups.items() if co==cohort and a==arm]
            if pubs:
                methods[arm]={metric:sum(sum(r['metrics'][metric] for r in v)/len(v) for v in pubs)/len(pubs) for metric in ('unique_hits','selected','measured','measured_precision','ef')}
        output[cohort]={'means':methods,'incumbent_minus_native_filled':methods['incumbent']['unique_hits']-methods['native_filled']['unique_hits']}
    wins=losses=ties=0
    for cohort,arm,pub in groups:
        if arm!='incumbent':continue
        a=groups[cohort,arm,pub];b=groups[cohort,'native_filled',pub]
        cells=lambda rr:{(r['cell']['actor'],r['cell']['draw']) for r in rr}
        if cells(a)!=cells(b) or len(a)!=4:raise ValueError('Incomplete matched publication')
        d=sum(r['metrics']['unique_hits'] for r in a)/4-sum(r['metrics']['unique_hits'] for r in b)/4
        wins+=d>0;losses+=d<0;ties+=d==0
    if reference and set(ref)!=seen:raise ValueError('Historical matrix coverage differs')
    return {'campaigns':len(rows),'cohorts':output,'publication_wins':wins,'publication_losses':losses,'publication_ties':ties,'historical_csv_verified':bool(reference)}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--summary',required=True);p.add_argument('--reference');a=p.parse_args()
    print(json.dumps(summarize(json.loads(Path(a.summary).read_text()),a.reference),indent=2))
