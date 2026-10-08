"""
candidate_pool.py

Optional "overgenerate, then filter by confidence" strategy for BioDiscoveryAgent.

Instead of asking the model for exactly --num_genes genes per round, ask for a
larger pool of candidates (--candidate_pool, e.g. 350), each with a 1-5
confidence score, then test only the --num_genes highest-confidence candidates.

Enabled with --candidate_pool N (N > 0). With N = 0 (the default) the original
code path in tools.py is used unchanged.

Gene validity rules are deliberately identical to the original code: a name is
valid only if it exactly matches a gene in the dataset's measured-gene list and
has not been tested already. That keeps an A/B comparison against the original
method about the new strategy only.

Example usage:

python research_assistant.py --task perturb-genes-brief \
    --model gemini:gemini-3.5-flash-lite \
    --reasoning_effort low \
    --candidate_pool 350 --candidate_min 300 --prompt_tries 8 \
    --run_name test --data_name IFNG --steps 1 --num_genes 128 \
    --log_dir gemini_pool
"""

import os
import re

import numpy as np
import pandas as pd

MIN_PROGRESS = 5     # a retry adding fewer new valid candidates than this is a "stall"
MAX_STALLS = 2       # consecutive stalls before finalizing early
DEFAULT_SCORE = 3    # used when the model omits or garbles a confidence score
MAX_INVALID_SHOWN = 100


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def normalize_completion(completion):
    """Strip visible scratchpad blocks and turn markdown 'Solution' headings
    (e.g. '## 3. Solution') into the plain 'Solution:' label."""
    completion = re.sub(r'<(thought|thinking|think)>[\s\S]*?</\1>', '', completion)
    # unclosed block (truncated output): drop everything from the opening tag on
    completion = re.sub(r'<(thought|thinking|think)>[\s\S]*$', '', completion)
    completion = re.sub(r'(?im)^\s*#{1,6}\s*\d*\.?\s*Solution\s*$', 'Solution:', completion)
    return completion


def extract_solution_text(completion):
    """Text after the LAST 'Solution:' label, or None if there isn't one.
    (Last rather than first, so an earlier mention of the word can't swallow
    the preceding prose into the gene list.)"""
    completion = normalize_completion(completion)
    if "Solution:" not in completion:
        return None
    return completion.rsplit("Solution:", 1)[1]


_NUMBERING = re.compile(r'^\s*\d+\s*[\.\)]\s*')
_SYMBOL = re.compile(r'^[\s\-\*\u2022]*[\[\(]?\s*([A-Za-z0-9][A-Za-z0-9_\-\.]*)')
_SCORE = re.compile(r'(?<![\w.])([1-5])(?!\d|\.\d)')


def _split_items(text):
    """Split on newlines, and on commas that are not inside parentheses."""
    items, cur, level = [], [], 0
    for ch in text:
        if ch == '\n':
            items.append(''.join(cur))
            cur, level = [], 0
        elif ch == '(':
            level += 1
            cur.append(ch)
        elif ch == ')':
            level = max(0, level - 1)
            cur.append(ch)
        elif ch == ',' and level == 0:
            items.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    items.append(''.join(cur))
    return [i.strip() for i in items if i.strip()]


def parse_scored_genes(solution_text):
    """Parse '1. GENE (5), 2. GENE2 (3), ...' into [(gene, score_or_None), ...]
    in listed order. Tolerates 'GENE: 5', 'GENE (5/5)', 'GENE (confidence: 4)',
    bullets, numbering, and newline- or comma-separated lists."""
    out = []
    for item in _split_items(solution_text):
        item = _NUMBERING.sub('', item, count=1)
        m = _SYMBOL.match(item)
        if not m:
            continue
        gene = m.group(1).rstrip('.-_')
        if not gene:
            continue
        sm = _SCORE.search(item[m.end():])
        out.append((gene, int(sm.group(1)) if sm else None))
    return out


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------

def select_top(candidates, k):
    """candidates: [(gene, score)] in listed order. Highest score first; the
    model's listed order breaks ties (stable)."""
    order = sorted(range(len(candidates)), key=lambda i: (-candidates[i][1], i))
    return [candidates[i][0] for i in order[:k]]


def _topup_prompt(candidates, invalid, pool_target):
    n_more = max(1, pool_target - len(candidates))
    s = ("\n You have so far proposed {} of the {} candidate genes requested for "
         "this round (each with a confidence score). These genes were: \n"
         .format(len(candidates), pool_target))
    s += str([g for g, _ in candidates])
    s += ("\n Please add {} more NEW candidate genes, each with a 1-5 confidence "
          "score, using the same response format as before. Do not repeat any "
          "gene listed above.".format(n_more))
    if invalid:
        s += (" The following names were not recognized as valid gene symbols "
              "for this screen, so do not use them: \n "
              + str(invalid[:MAX_INVALID_SHOWN]))
    return s


# --------------------------------------------------------------------------
# One round
# --------------------------------------------------------------------------

def run_candidate_pool_round(prompt, log_file, curr_step, gene_sampled,
                             measured_genes, all_hit_genes, log_dir, args,
                             complete_fn):
    """Run one round in candidate-pool mode and save sampled_genes_{step+1}.npy.

    complete_fn is tools.complete_text (passed in so this module has no heavy
    imports). Returns the list of genes selected for testing this round.
    """
    num_genes = args.num_genes
    pool_target = args.candidate_pool
    pool_min = max(num_genes, min(args.candidate_min, pool_target))

    measured_set = set(measured_genes)
    tested_set = set(gene_sampled)
    hit_set = set(all_hit_genes)

    call_kwargs = {"max_tokens_to_sample": args.pool_max_tokens}
    if getattr(args, "reasoning_effort", None) is not None:
        call_kwargs["reasoning_effort"] = args.reasoning_effort

    candidates = []          # [(gene, score)] in first-seen order
    seen = set()
    invalid, invalid_seen = [], set()
    n_defaulted = 0
    stalls = 0
    prompt_try = prompt

    for itr in range(args.prompt_tries):
        completion = complete_fn(prompt_try, model=args.model, log_file=log_file,
                                 **call_kwargs)
        solution = extract_solution_text(completion)
        if solution is None:
            print(itr, 'Invalid output')
            prompt_try = prompt + (_topup_prompt(candidates, invalid, pool_target)
                                   if candidates else '')
            continue

        pairs = parse_scored_genes(solution)
        n_new = n_dropped = 0
        for gene, score in pairs:
            if gene in measured_set:
                if gene in tested_set or gene in seen:
                    continue
                seen.add(gene)
                if score is None:
                    score = DEFAULT_SCORE
                    n_defaulted += 1
                candidates.append((gene, score))
                n_new += 1
            else:
                n_dropped += 1
                if gene not in invalid_seen:
                    invalid_seen.add(gene)
                    invalid.append(gene)

        print('Dropped genes:', n_dropped)
        print('New genes predicted:', n_new)
        print("[pool step {}] try {}: parsed {} items, pool now {} "
              "(need >= {}, ask {})".format(curr_step, itr, len(pairs),
                                            len(candidates), pool_min, pool_target))

        if len(candidates) >= pool_min:
            break
        stalls = stalls + 1 if n_new < MIN_PROGRESS else 0
        if stalls >= MAX_STALLS and len(candidates) >= num_genes:
            print("[pool step {}] diminishing returns, finalizing with {} "
                  "candidates".format(curr_step, len(candidates)))
            break
        prompt_try = prompt + _topup_prompt(candidates, invalid, pool_target)

    out_path = os.path.join(log_dir, 'sampled_genes_{}.npy'.format(curr_step + 1))

    if not candidates:
        print("WARNING [pool step {}]: no valid candidates after {} tries; "
              "nothing tested this round.".format(curr_step, args.prompt_tries))
        np.save(out_path, gene_sampled)
        return []

    selected = select_top(candidates, num_genes)
    if len(candidates) < num_genes:
        print("WARNING [pool step {}]: only {} valid candidates (< num_genes={}); "
              "testing all of them.".format(curr_step, len(candidates), num_genes))
    elif len(candidates) < pool_min:
        print("NOTE [pool step {}]: pool of {} is below the requested minimum {}; "
              "selecting the top {} anyway.".format(curr_step, len(candidates),
                                                    pool_min, num_genes))
    np.save(out_path, list(gene_sampled) + selected)

    # ---- analysis-only outputs (never shown to the model) ----
    sel_set = set(selected)
    df = pd.DataFrame([{"gene": g, "confidence": s, "listed_position": i,
                        "selected": g in sel_set, "is_hit": g in hit_set}
                       for i, (g, s) in enumerate(candidates)])
    df.to_csv(os.path.join(log_dir, 'step_{}_candidates.csv'.format(curr_step)),
              index=False)

    print("[pool step {}] pool={} selected={} invalid names seen={} "
          "scores defaulted={}".format(curr_step, len(candidates), len(selected),
                                       len(invalid), n_defaulted))
    summary = df.groupby("confidence").agg(n=("gene", "size"),
                                           selected=("selected", "sum"),
                                           hits=("is_hit", "sum"))
    print(summary.sort_index(ascending=False).to_string())
    print("hits among selected: {}/{} | hits left in unselected pool: {}/{}".format(
        int(df[df.selected].is_hit.sum()), len(selected),
        int(df[~df.selected].is_hit.sum()), len(candidates) - len(selected)))
    return selected
