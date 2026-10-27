"""Offline action replay of saved external campaigns using the packaged method.

Requires locally trusted research snapshot and released AssayFormer checkpoint.
No experiment runner, credentials, network client or provider call is used.
Saved language proposals are replayed only against their exact preceding history.
This does not sample a fresh LLM or certify the prompt-generation path.
"""
import argparse, gzip, hashlib, json, sys
from pathlib import Path
import numpy as np
import torch
from threadpoolctl import threadpool_limits
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from assayloop_comparison.assets import HistoricalBank
from assayloop_comparison.adapter import SourceState, run_campaign
from assaybench.data.screens import ScreenRecord
from assayloop.tasks.gene_batch import AssayBenchGeneBatchTask
from assayloop.models.amortized_ranker import AmortizedRankerModel


def load(p):
    with (gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)) as f:
        return json.load(f)


def digest(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def observations(rows):
    return [dict(candidate=o.candidate, label=o.label, metadata=o.metadata) for o in rows]


class ExternalTask(AssayBenchGeneBatchTask):
    def context(self):
        result = super().context()
        result.pop('total_hits', None)
        return result


class RecordedActor:
    def __init__(self, run):
        self.run = run
    def reset(self):
        self.calls = 0
    def suggest(self, **kw):
        assert self.calls < 3
        expected = self.run['rounds'][:self.calls]
        assert len(kw['history']) == len(expected)
        for actual, saved in zip(kw['history'], expected):
            assert list(actual.acquired_batch) == saved['accepted'], 'Changed selection history'
            assert json.dumps(observations(actual.new_observations), sort_keys=True) == json.dumps(saved['observations'], sort_keys=True), 'Changed feedback'
        assert 'total_hits' not in kw.get('task_context', {})
        proposal = self.run['language_records'][self.calls]['proposal']
        self.calls += 1
        assert set(proposal) <= set(kw['candidates'])
        return proposal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--research-root', type=Path, required=True)
    ap.add_argument('--checkpoint', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--all-paths', action='store_true', help='Default: one actor/draw per assay and arm')
    args = ap.parse_args()
    root = args.research_root.resolve()
    E = root / 'Experiments'
    package = Path(__file__).resolve().parent
    rows = load(package / 'recorded_results/external_campaigns.json')
    source_meta = load(package / 'recorded_results/external_sources.json')
    pool = load(E / 'native_subset_training_2026-10-17/UNIVERSE_validation.json')
    archive_path = E / 'full_native_benchmark_2026-10-17/train_1349_archive.npz'
    factor_path = E / 'native_frozen_readout_2026-10-17/artifacts/pub_subspaces_140.pt'
    # These are the project-owned, previously loaded fitted numerical artifacts.
    # Never substitute untrusted pickle/checkpoint files.
    factors = torch.load(factor_path, map_location='cpu', weights_only=False)
    with np.load(archive_path, allow_pickle=False) as archive:
        pubs = archive['pubs'].astype(str)
        genes = archive['u_genes'].astype(str).tolist()
        M = archive['M_u'].astype(float)
    pub_masks = np.asarray([M[pubs == p].mean(0) for p in sorted(set(pubs))])
    q = (1 + pub_masks.sum(0)) / (2 + len(pub_masks))
    global_keys = sorted(set(factors) - set(pubs))
    factor = factors[global_keys[0]]
    assert all(np.array_equal(factor['V_res32'], factors[k]['V_res32']) for k in global_keys)
    scales = load(E / 'boundary_calibration_e209_2026-10-25/SCALAR_PRIORS.json')
    c = scales[global_keys[0]]
    chosen, seen = [], set()
    for row in rows:
        key = (row['screen'], row['method'])
        if not args.all_paths and key in seen:
            continue
        seen.add(key)
        chosen.append(row)
    outcomes = []
    for row in chosen:
        meta = source_meta[row['screen']]
        base = E / meta['source_directory']
        target_path = base / ('TARGET.json' if (base / 'TARGET.json').exists() else 'TARGETS.json')
        targets = load(target_path)
        if 'screen' in targets:
            target, description = targets, base / 'DESCRIPTION.npy'
        else:
            key = next(k for k,v in targets.items() if v['screen_id'] == row['screen'])
            target, description = targets[key], base / ('DESCRIPTION_' + key + '.npy')
        pub = str(target['publication_id'])
        assert pub not in set(pubs)
        sc = ScreenRecord(**target['screen'])
        desc = np.load(description, allow_pickle=False)
        source = AmortizedRankerModel(checkpoint=args.checkpoint, ckpt_file='model_last.pt',
            device='cpu', desc_emb_by_name={row['screen']: desc})
        state = SourceState(source, desc)
        vi = state.ids(genes)
        residual = np.asarray(factor['V_res32'])[vi]
        norm_desc = desc.astype(float) / np.linalg.norm(desc.astype(float))
        mean = .1 * (norm_desc @ np.asarray(factor['W_text32'], dtype=float))
        bank = HistoricalBank(genes, residual, mean, q, c, pub, tuple(sorted(set(pubs))))
        path = root / row['artifact']
        assert digest(path) == row['sha256'], 'Saved campaign changed'
        saved = load(path)
        actor = RecordedActor(saved)
        task = ExternalTask(sc, universe_genes=pool, seed=0)
        arm = 'incumbent' if row['method'] == 'HSA' else 'native_filled'
        result, audit = run_campaign(arm=arm, task=task, source_model=source,
            screen_id=row['screen'], draw=row['draw'], source_state=state, bank=bank,
            actor_acquisition=actor, allow_live_calls=True)
        for i,(actual, expected) in enumerate(zip(result.history, saved['rounds'])):
            assert list(actual.acquired_batch) == expected['accepted'], (row['screen'], arm, i+1, 'action mismatch')
            assert json.dumps(observations(actual.new_observations), sort_keys=True) == json.dumps(expected['observations'], sort_keys=True)
            if arm == 'incumbent':
                np.testing.assert_allclose(audit['model'][i]['mode'], saved['model_audit'][i]['posterior_mode'], atol=2e-5, rtol=0)
        assert len(result.history) == 10 and actor.calls == 3
        hits = sum(o.metadata['in_library'] and o.label['hit'] for r in result.history for o in r.new_observations)
        assert hits == row['hits']
        outcomes.append(dict(source=row['source'], screen=row['screen'], method=row['method'],
            actor=row['actor'], draw=row['draw'], decisions=10, hits=int(hits), all_selected_ids_exact=True,
            saved_campaign_sha256=row['sha256']))
        print(json.dumps(outcomes[-1]), flush=True)
    report = dict(status='passed', kind='exact-history saved-proposal full-loop replay',
        fresh_llm_calls=0, campaigns=len(outcomes), decisions=10*len(outcomes), rows=outcomes,
        factor_sha256=digest(factor_path), archive_sha256=digest(archive_path),
        checkpoint_sha256=digest(args.checkpoint/'model_last.pt'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    torch.set_num_threads(1)
    with threadpool_limits(limits=1):
        main()
