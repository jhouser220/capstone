"""Small, budgeted retrospective BDA-style pilot. No paid work on import."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import sqlite3
import time
import urllib.request
from datetime import datetime, timezone

import numpy as np

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
DATA = Path(os.environ.get('BDA_DATA_DIR', PROJECT.parents[1] / 'datasets')).expanduser().resolve()
SMALL = 'qwen/qwen3-8b'
LARGER = 'qwen/qwen3-14b'
PILOT_CAP = 1.0
BATCH = 32
SYSTEM = '''You are selecting gene perturbations for a retrospective biological discovery task.
Use only the supplied task, available gene names, observed experiments, and your prior knowledge.
Unobserved outcomes are unknown. Treat critiques as suggestions, not ground truth.
Use the requested output type: either exact gene symbols or integer candidate IDs. Select only from the current candidate menu. Historical measurements are evidence, never available actions. Never present a prediction as a measured outcome.
Return the requested JSON object only. /no_think'''


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False))
    tmp.replace(path)


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


class Budget:
    """Atomic reservations count unresolved requests conservatively at their maximum."""
    def __init__(self, path=None, cap=PILOT_CAP):
        ledger_root = Path(os.environ.get('BDA_BUDGET_ROOT', PROJECT / 'Experiments'))
        path = Path(path) if path is not None else ledger_root / ROOT.name / 'budget.sqlite'
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path, self.cap = path, cap
        with sqlite3.connect(path) as db:
            db.execute('CREATE TABLE IF NOT EXISTS calls (id TEXT PRIMARY KEY, status TEXT, reserve REAL, actual REAL, meta TEXT)')

    def reserve(self, call_id, amount, meta):
        with sqlite3.connect(self.path, timeout=60) as db:
            db.execute('BEGIN IMMEDIATE')
            spent = db.execute('SELECT COALESCE(SUM(COALESCE(actual,reserve)),0) FROM calls').fetchone()[0]
            if db.execute('SELECT 1 FROM calls WHERE id=?', (call_id,)).fetchone():
                raise RuntimeError('Existing unresolved or completed call; inspect before retry: ' + call_id)
            other = 0.0
            for ledger in Path(os.environ.get('BDA_BUDGET_ROOT', PROJECT / 'Experiments')).glob('*/budget.sqlite'):
                if ledger.resolve() == self.path.resolve(): continue
                with sqlite3.connect('file:'+str(ledger)+'?mode=ro', uri=True) as prior:
                    other += prior.execute('SELECT COALESCE(SUM(COALESCE(actual,reserve)),0) FROM calls').fetchone()[0]
            if other + spent + amount >= 40.0:
                raise RuntimeError('Project working allocation exhausted; preserve reserve below $50')
            if spent + amount > self.cap:
                raise RuntimeError(f'Pilot cap would be exceeded: accounted={spent:.6f}, next={amount:.6f}')
            db.execute('INSERT INTO calls VALUES (?,?,?,?,?)', (call_id, 'pending', amount, None, json.dumps(meta)))

    def settle(self, call_id, cost, status):
        with sqlite3.connect(self.path) as db:
            db.execute('UPDATE calls SET actual=?,status=? WHERE id=?', (cost, status, call_id))

    def summary(self):
        with sqlite3.connect(self.path) as db:
            rows = db.execute('SELECT id,status,reserve,actual,meta FROM calls').fetchall()
        return {'cap_usd': self.cap, 'reported_cost_usd': sum(x[3] or 0 for x in rows),
                'accounted_usd': sum(x[3] if x[3] is not None else x[2] for x in rows),
                'calls': len(rows), 'unresolved': [x[0] for x in rows if x[3] is None]}


def api_get(endpoint):
    req = urllib.request.Request('https://openrouter.ai/api/v1/' + endpoint,
                                 headers={'Authorization': 'Bearer ' + os.environ['OPENROUTER_API_KEY']})
    with urllib.request.urlopen(req, timeout=45) as response:
        return json.load(response)


def request_llm(call_id, prompt, model=SMALL, seed=101, max_tokens=1024):
    call_id = 'rich_'  + call_id
    folder = ROOT / 'calls' / call_id
    response_path = folder / 'response.json'
    if response_path.exists():
        saved=json.loads((folder/'request.json').read_text())
        if saved['model']!=model or saved['messages'][-1]['content']!=prompt or saved.get('seed')!=seed or saved.get('max_tokens')!=max_tokens or saved['messages'][0]['content']!=SYSTEM:
            raise RuntimeError('Cached request differs from current configuration: '+call_id)
        result = json.loads(response_path.read_text())
        if 'choices' in result:
            return result
        raise RuntimeError('Previous failed response needs inspection: ' + call_id)
    if os.environ.get('BDA_ENABLE_PAID_CALLS') != '1' or not os.environ.get('BDA_BUDGET_ROOT'):
        raise RuntimeError('Paid calls disabled. Set BDA_ENABLE_PAID_CALLS=1 and BDA_BUDGET_ROOT to the reconciled shared experiment-ledger directory before opting in.')
    payload = {'model': model, 'messages': [{'role': 'system', 'content': SYSTEM}, {'role': 'user', 'content': prompt}],
               'temperature': 0.5, 'seed': seed, 'max_tokens': max_tokens,
               'response_format': {'type': 'json_object'},
               'provider': {'only': ['google-ai-studio' if model.startswith('google/') else 'alibaba'], 'allow_fallbacks': False, 'require_parameters': True,
                            'max_price': {'prompt': 0.5, 'completion': 3.0, 'request': 0}},
               'usage': {'include': True}}
    if model == SMALL or model.startswith('google/'): payload['reasoning'] = {'enabled': False}
    # UTF-8 byte count + generous framing margin bounds byte-level tokenizers.
    # Reserve at provider price ceilings, double output allowance, plus 10% margin.
    upper_input = len(json.dumps(payload, ensure_ascii=False).encode()) + 4096
    reservation = 1.10 * (upper_input * 0.5e-6 + 2 * max_tokens * 3e-6)
    budget = Budget()
    budget.reserve(call_id, reservation, {'model': model, 'prompt_sha256': digest(prompt), 'time': stamp()})
    write_json(folder / 'request.json', payload)
    req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions',
                                 data=json.dumps(payload).encode(), method='POST',
                                 headers={'Authorization': 'Bearer ' + os.environ['OPENROUTER_API_KEY'],
                                          'Content-Type': 'application/json', 'X-Title': 'UCB Capstone Feedback Pilot'})
    started = time.monotonic()
    try:
        with urllib.request.urlopen(req, timeout=180) as response:
            result = json.load(response)
        write_json(response_path, result)
        usage = result.get('usage', {})
        cost = usage.get('cost')
        if cost is None and result.get('id'):
            details = api_get('generation?id=' + result['id'])
            write_json(folder / 'generation.json', details)
            cost = details.get('data', {}).get('total_cost')
        if cost is None:
            budget.settle(call_id, None, 'unknown_cost')
            raise RuntimeError('No cost returned; reserved upper bound retained')
        budget.settle(call_id, float(cost), 'completed')
        write_json(folder / 'accounting.json', {'time': stamp(), 'latency_seconds': time.monotonic()-started,
                                               'reserved_usd': reservation, 'usage': usage})
        if float(cost) > reservation:
            raise RuntimeError('Cost exceeded conservative reservation; stop and audit')
        if 'choices' not in result:
            raise RuntimeError('API response has no choices')
        print(json.dumps({'call': call_id, 'cost': cost, 'input': usage.get('prompt_tokens'),
                          'output': usage.get('completion_tokens')}), flush=True)
        return result
    except Exception as exc:
        # Never dump headers, environment or credentials. Unknown bill remains reserved.
        write_json(folder / 'error.json', {'type': type(exc).__name__, 'http_status': getattr(exc, 'code', None),
                                          'time': stamp(), 'note': 'No automatic retry; inspect request state.'})
        raise


def parse_result(result):
    content=(result['choices'][0]['message'].get('content') or '').strip()
    try:
        value=json.loads(content)
    except (ValueError,TypeError):
        value={}
        for i,char in enumerate(content):
            if char in '{[':
                try:value=json.JSONDecoder().raw_decode(content[i:])[0];break
                except ValueError:continue
    return value if isinstance(value,(dict,list)) else {}
