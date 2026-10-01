"""No API calls: check outcome isolation, valid choices, reports and paid-call gate."""
import os
import sys
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'Experiments/cross_dataset_numeric_2026-09-29'))
import run
import llm_study
import paid_api


def main():
    data, oracle = run.dataset('IFNG')
    # The numerical interface accepts revealed measurements, not the oracle.
    initial = run.old.base.initial(data, 8101)
    remaining = run.np.array([i for i in range(len(data['genes'])) if i not in set(initial)])
    query = remaining[:96]
    altered = {k: v.copy() for k, v in oracle.items()}
    altered['y'][remaining] = 1 - altered['y'][remaining]
    altered['effect'][remaining] += 100
    history = run.old.base.reveal(initial, oracle)
    changed_history = run.old.base.reveal(initial, altered)
    for method in ['prior', 'gp_iso', 'full:effect_median', 'full:hit_logistic']:
        a = run.old.predict(data, history, query, method)
        b = run.old.predict(data, changed_history, query, method)
        assert run.np.array_equal(a['p'], b['p']), method
        selected = run.old.base.choose(a, query)
        assert len(selected) == len(set(selected)) == 32
        assert not set(selected) & set(initial)
        assert set(selected) <= set(query)
    report, _ = llm_study.report_for(data, history, query[:4], 'gp_iso')
    observed = {str(data['genes'][i]) for i in initial}
    assert all(e['gene'] in observed for c in report['candidates'] for e in c['observed_neighbours'])
    for api in [run.old.api, paid_api]:
        with patch.dict(os.environ, {'BDA_ENABLE_PAID_CALLS': '0'}):
            with patch.object(api.urllib.request, 'urlopen', side_effect=AssertionError('Unexpected network call')):
                try:
                    api.request_llm('packaging_gate_check', 'test')
                except RuntimeError as exc:
                    assert 'Paid calls disabled' in str(exc)
                else:
                    raise AssertionError('Missing paid-call gate')
    print('PASS: four learners, hidden-outcome mutation, valid selections, observed-only examples, both paid-call gates. No API calls.')


if __name__ == '__main__':
    main()
