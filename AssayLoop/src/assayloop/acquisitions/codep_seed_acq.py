"""Co-dependency seed acquisition: AssayFormer fused with DepMap co-essentiality.

Motivation. AssayFormer scores a gene as ``û·v_g + b_g`` where ``v_g`` is a
10-number BPMF embedding that is the SAME for every screen. Offline screening
(``scripts/diagnostics/screen_profiles.py``) found that a gene's *profile* of
knockout effects across DepMap cell lines carries hit-relevant information that
``span(v_g)`` does not: appended as a feature, its margin over a width-matched
noise block was +0.35 EF [+0.11, +0.61] on validation and +1.36 [+0.69, +2.06] on
test, both bootstrapped over screen families, while a "gene is in DepMap at all"
confound control came out at zero. A single *scalar* summary of the same data
(mean essentiality) was null, so the profile structure is what matters.

The feature. Each round, take the genes revealed as hits so far, average their
row-centred unit-normalised essentiality profiles into a query vector ``q``, and
score every untested gene by ``corr(profile(g), q)`` -- genes that co-vary across
cell lines with what already worked in this screen. This is pathway propagation
with a data-derived graph instead of a curated one, and it is conditioned on the
screen, which a fixed gene embedding cannot be.

Two deliberate choices, both to keep the baseline nested so the worst case is
parity rather than regression:

* **Round 1 has no hits**, so there is no query vector and the acquisition is
  exactly ``greedy`` on AssayFormer's scores. The same holds for any screen that
  never reveals a hit.
* **Genes absent from DepMap are made NEUTRAL, not unranked.** Under RRF, simply
  omitting them is not neutral: possessing a co-dependency rank is itself worth up
  to ``weight/RRF_K``, so the 3.4k genes missing from the atlas would be
  systematically outranked by covered ones — the coverage confound re-entering
  through the fusion rule. They instead receive the MEAN co-dependency
  contribution of the covered genes, so coverage alone conveys no advantage either
  way. (Under ``zsum`` a standardised zero is already the covered mean, so that
  path was correct as written.)

``weight=0`` reproduces ``greedy`` exactly under either fusion rule.

Second scorer (optional): screen-kNN. ``v_g`` is a 10-number compression of 1349
training screens, and §4.9 of the handoff found only ~3.5% of a screen's hit
structure survives it. Given the raw training hit matrix, this scores a gene by how
often it is a hit in the training screens whose hit sets most resemble what THIS
screen has revealed. Offline the two scorers are near-perfectly additive -- their
combined oracle margin (+0.260 validation, +0.733 test) is about the sum of the
parts -- so they carry largely independent information. Using training labels at
inference is not leakage: they are training data and the folds are disjoint.

Hierarchical shrinkage (optional). The query direction is estimated from the ~20-30
hits revealed so far, in 808-dimensional space. Three attempts to extract more
structure from those few profiles all lost to the plain difference of means, which
says the binding constraint is sample size. The 1349 training screens have COMPLETE
labels, so each one's own optimal direction can be computed offline from hundreds of
hits (``scripts/diagnostics/build_dir_bank.py``). Given ``dir_bank``, the directions
of the most similar training screens -- the same similarity the screen-kNN scorer
already computes -- are averaged into a prior ``p`` and the live estimate is pulled
toward it:

    q = normalise((1 - alpha) * q_own + alpha * p)

``shrink_alpha=0`` is the default and reproduces the unshrunk query bit for bit, so
the baseline stays nested. ``alpha=1`` discards the live screen's own direction
entirely and is the informative extreme of the sweep, not a recommendation. Note
this borrows a DIRECTION, not a gene list: unlike the screen-kNN scorer it can
promote a gene that was never a hit in any training screen, provided it co-varies
with the genes that were.

Fusion. ``rrf`` (default) converts both to ranks and sums ``w/(k+rank)``, which
needs no calibration between an unbounded bilinear score and a correlation in
[-1, 1]. ``zsum`` standardises both and adds them, which is sharper when the
weight is right and worse when it is not. Offline fitting of that weight on the
first 300 labels was measured to be too noisy to trust (a pure confound control
scored +0.58 EF under the same procedure), which is why a fixed weight and a
scale-free rule are preferred here over anything learned online.
"""

from __future__ import annotations

import logging
import re
import threading
from pathlib import Path
from typing import Any

import numpy as np

from assaybench.core.acquisition import AcquisitionFunction
from assaybench.core.types import ModelPrediction, StepRecord

log = logging.getLogger("assayloop.acquisitions.codep_seed")

RRF_K = 60.0

# Shared across screens in a process: the paper Screen-kNN state is ~230 MB and
# identical for every screen in a sweep, so one copy, not one per acquisition. The
# lock matters under --parallel: without it every worker starts its own build on the
# first round and the peak is one copy PER THREAD, not one per process.
_PKNN_CACHE: dict = {}
_PKNN_LOCK = threading.Lock()
# The universe is supposed to be one fixed pool per sweep, so this dict is supposed to
# hold ONE entry. It held one per screen once (see `_pknn_scores`) and 218 x 230 MB
# took the process out with the OOM killer. Cap it: a broken invariant now costs a
# rebuild and a warning instead of the machine.
_PKNN_CACHE_MAX = 2

# Same reasoning for the GenePT table: 139 MB, identical for every screen, and
# one acquisition is built PER SCREEN -- without this a 218-screen sweep reads it
# 218 times.
_GPT_CACHE: dict = {}
_GPT_LOCK = threading.Lock()


class CodepSeedAcquisition(AcquisitionFunction):
    """Greedy on AssayFormer's scores, fused with DepMap co-essentiality.

    Args:
        profile_npz: ``.npz`` written by ``scripts/diagnostics/depmap_to_npz.py``
            holding ``values`` (cell lines x genes), ``rows``, ``cols``.
        weight: this scorer's weight in the fusion. ``0`` reproduces ``greedy``.
        fuse: ``"rrf"`` (scale-free, default) or ``"zsum"``.
        min_seeds: hits required before the feature switches on. Below this the
            query vector is dominated by one or two genes.
        seed: tie-breaking RNG seed, matching ``GreedyFromModel``.
        dir_bank: ``.npz`` from ``scripts/diagnostics/build_dir_bank.py``: one
            offline-optimal profile direction per training screen.
        shrink_alpha: weight on the borrowed prior direction, in [0, 1]. ``0``
            (default) leaves the live estimate untouched.
        shrink_k: training screens averaged into the prior.
        shrink_drop_top: nearest neighbours DISCARDED before averaging, to test
            whether any gain depends on near-duplicate screens.
        shrink_min_hits: a training screen needs this many atlas-covered hits for
            its direction to enter the bank at all.
        shrink_shuffle: >= 0 permutes the bank with this seed, breaking the
            screen-to-direction correspondence while keeping the distribution.
            A CONTROL ARM: it isolates "borrowed from SIMILAR screens" from
            "regularised toward a generic direction".
        pknn_gate: hold the paper's Screen-kNN silent until ``min_seeds`` confirmed
            hits exist, which is how every scorer of ours behaves. ``True``
            (default) makes its arms NESTED in the baseline, so round 1 is plain
            greedy and any difference is attributable to the estimator. ``False``
            is the ungated behaviour the pre-2026-10-28 ``sweep-pkv-*`` arms ran:
            with no observations their score is the pure global hit prior, which
            reorders round 1 toward famous genes before anything is known. Keep
            ``False`` to reproduce those arms and to price the round-1 prior alone.
            See §4.56(e).
    """

    def __init__(
        self,
        *,
        profile_npz: str | Path,
        weight: float = 1.0,
        fuse: str = "rrf",
        min_seeds: int = 2,
        seed: int = 0,
        train_hits: str | Path | None = None,
        knn_weight: float = 0.0,
        knn_k: int = 50,
        query: str = "diffmeans",
        min_nonhits: int = 10,
        max_labels: int = 600,
        profile_npz2: str | Path | None = None,
        weight2: float = 0.0,
        dir_bank: str | Path | None = None,
        shrink_alpha: float = 0.0,
        shrink_k: int = 25,
        shrink_drop_top: int = 0,
        shrink_min_hits: int = 50,
        shrink_shuffle: int = -1,
        ppi_npz: str | Path | None = None,
        ppi_weight: float = 0.0,
        ppi_degree_residual: bool = False,
        gpt_npz: str | Path | None = None,
        gpt_weight: float = 0.0,
        gpt_k: int = 10,
        pknn_weight: float = 0.0,
        pknn_train: str = "train",
        pknn_tau: float = 0.1,
        pknn_kappa: float = 100.0,
        pknn_sim: str = "agreement",
        pknn_k: int = 0,
        pknn_shrink_m: float = 0.0,
        pknn_gate: bool = True,
        drop_pcs: int = 0,
        amp_power: float = 0.0,
        gene_locus: str | Path | None = None,
        prox_level: str = "arm",
        oracle_seeds: list[str] | None = None,
        inject_genes: list[str] | None = None,
        inject_per_round: int = 0,
        inject_rounds: int = 3,
        af_weight: float = 1.0,
        hsa_weight: float = 0.0,
        hsa_scorer: Any = None,
    ):
        if fuse not in ("rrf", "zsum"):
            raise ValueError(f"fuse must be 'rrf' or 'zsum', got {fuse!r}")
        self.weight = float(weight)
        self.weight2 = float(weight2)
        # AssayFormer's own rank enters the ensemble at `af_weight`. 1.0 is the
        # historical behaviour and is bit-for-bit what every existing arm ran.
        # 0.0 drops it, which is only sensible when another term already carries
        # AssayFormer internally -- see `hsa_scorer`.
        self.af_weight = float(af_weight)
        # A seventh fusion term, supplied as a callable rather than loaded from a
        # file, because its state (a 42-D historical bank, a frozen source-state
        # encoder) belongs to the caller and must not become a dependency of this
        # package. `hsa_scorer(observations, candidates, task_context)` returns a
        # per-candidate score, as a dict keyed by candidate or an aligned array.
        # It is called EXACTLY ONCE PER ROUND and only when `hsa_weight != 0`,
        # because a stateful scorer may count rounds (the reference implementation
        # raises outside steps 1..10).
        self.hsa_weight = float(hsa_weight)
        self._hsa_scorer = hsa_scorer
        if self.hsa_weight != 0.0 and hsa_scorer is None:
            raise ValueError("hsa_weight != 0 requires an hsa_scorer")
        self.fuse = fuse
        self.min_seeds = int(min_seeds)
        self.seed = int(seed)
        self.knn_weight = float(knn_weight)
        self.knn_k = int(knn_k)
        self.query = str(query)
        self.min_nonhits = int(min_nonhits)
        self.max_labels = int(max_labels)
        self.shrink_alpha = float(shrink_alpha)
        self.shrink_k = int(shrink_k)
        self.shrink_drop_top = int(shrink_drop_top)
        self.shrink_min_hits = int(shrink_min_hits)
        self.shrink_shuffle = int(shrink_shuffle)
        if self.query not in ("diffmeans", "centroid"):
            raise ValueError(f"query must be 'diffmeans' or 'centroid', got {query!r}")
        if not 0.0 <= self.shrink_alpha <= 1.0:
            raise ValueError(f"shrink_alpha must be in [0, 1], got {shrink_alpha!r}")
        self._last_trace: dict[str, Any] = {}
        self._pknn_universe = None
        self._pknn_at: dict[str, int] = {}

        # The training hit matrix drives BOTH the screen-kNN scorer and the
        # neighbour selection shrinkage borrows from, so load it for either.
        self._K_csc = self._K_csr = None
        self._K_col: dict[str, int] = {}
        self._K_size = None
        self._K_screens: list[str] = []
        want_train = self.knn_weight != 0.0 or self.shrink_alpha > 0.0
        if train_hits is not None and want_train:
            from scipy.sparse import load_npz
            pre = Path(train_hits)
            K = load_npz(pre.with_suffix(".mat.npz"))
            meta = np.load(pre.with_suffix(".npz"), allow_pickle=True)
            self._K_csc, self._K_csr = K.tocsc(), K.tocsr()
            self._K_col = {str(g): j for j, g in enumerate(meta["gene_names"])}
            self._K_size = np.sqrt(np.maximum(
                np.asarray(K.sum(axis=1)).ravel(), 1.0))
            self._K_screens = [str(s) for s in meta["screen_names"]]
            log.info("codep_seed: screen-kNN over %d training screens x %d genes",
                     K.shape[0], K.shape[1])

        # Direction bank, aligned to the rows of K by screen NAME rather than by
        # position: both come from `load_screens(target_set="train")` and so agree
        # today, but a silent misalignment here would look like a weak result
        # rather than a bug.
        self._bank = None
        if dir_bank is not None and self.shrink_alpha > 0.0:
            if self._K_csc is None:
                raise ValueError("shrink_alpha > 0 needs train_hits for the "
                                 "neighbour similarity as well as dir_bank")
            z = np.load(Path(dir_bank), allow_pickle=True)
            D = np.asarray(z["dirs"], dtype=np.float64)
            names = [str(s) for s in z["screen_names"]]
            nh = np.asarray(z["n_hits"])
            at = {n: i for i, n in enumerate(names)}
            B = np.zeros((len(self._K_screens), D.shape[1]))
            for i, nm in enumerate(self._K_screens):
                j = at.get(nm)
                # an under-determined direction is dropped, not down-weighted:
                # its row stays zero and _prior_dir skips it
                if j is not None and nh[j] >= self.shrink_min_hits:
                    B[i] = D[j]
            self._bank = B
            usable = int((np.linalg.norm(B, axis=1) > 0).sum())
            if self.shrink_shuffle >= 0:
                # CONFOUND CONTROL. Directions share a large common component
                # (pairwise cosine ~ +0.43), so shrinking toward an ARBITRARY
                # training screen may help on its own. Permuting which direction
                # belongs to which screen keeps the bank's marginal distribution
                # and destroys only the neighbour correspondence.
                perm = np.random.default_rng(self.shrink_shuffle).permutation(len(B))
                self._bank = B[perm]
                log.warning("codep_seed: dir_bank SHUFFLED (seed %d) -- this is a "
                            "control arm, not a method", self.shrink_shuffle)
            log.info("codep_seed: direction bank %d usable of %d training screens, "
                     "alpha=%.2f k=%d drop_top=%d",
                     usable, len(B), self.shrink_alpha, self.shrink_k,
                     self.shrink_drop_top)

        # ORACLE CONTROL (§4.49), never a deployable setting. Genes known from the
        # answer key to be hits, supplied to the scorers as if already confirmed, to
        # bound "would a better nominator -- an LLM, say -- help by handing the kernel
        # seeds AssayFormer cannot reach?" (§4.48). They shape the query DIRECTION and
        # are REMOVED FROM THE CANDIDATE POOL, so they cannot be picked for credit;
        # otherwise this would measure "we gave it free hits", which is trivially
        # positive and says nothing.
        self.oracle_seeds = [str(g) for g in (oracle_seeds or [])]
        # INJECTION (§4.50), a strictly better-posed version of §4.49's oracle. Instead
        # of handing the scorers seeds from outside the campaign -- which forced those
        # seeds OUT of the candidate pool and cost the arm any credit for finding them
        # (§4.49's one confound) -- this substitutes them INTO the early batches as real
        # picks. Budget stays exactly 1000, the protocol stays the paper's 10x100, and
        # nothing is confiscated. It models what a perfect nominator actually does:
        # occupy a few early slots with genes the model would never choose.
        #
        # The injected hits credit the score, so this is only interpretable as a 2x2:
        # run it with the kernel ON and OFF and difference the differences. The free-hit
        # credit is identical in both arms and cancels exactly.
        self.inject_genes = [str(g) for g in (inject_genes or [])]
        self.inject_per_round = int(inject_per_round)
        self.inject_rounds = int(inject_rounds)
        self._injected: set[str] = set()
        self.drop_pcs = int(drop_pcs)
        self.amp_power = float(amp_power)
        self.ppi_weight = float(ppi_weight)
        self.ppi_degree_residual = bool(ppi_degree_residual)
        self._P, self._row, self._amp = self._load_profiles(
            profile_npz, drop_pcs=self.drop_pcs, want_amp=True)
        log.info("codep_seed: %d genes x %d dims from %s",
                 self._P.shape[0], self._P.shape[1], Path(profile_npz).name)

        # CHROMOSOMAL-PROXIMITY CORRECTION. Neighbouring genes show inflated
        # co-essentiality because a CRISPR knockout's fitness effect carries the local
        # copy-number state, so two genes in one amplified region look co-dependent for
        # a reason that is about the region, not about either gene. The co-essentiality
        # literature excludes proximal pairs before calling a network; we did not, and
        # this makes it an option. Default off: `gene_locus=None` leaves every score
        # bit-identical.
        self.prox_level = str(prox_level)
        if self.prox_level not in ("band", "arm", "chrom"):
            raise ValueError(f"prox_level must be band|arm|chrom, got {prox_level!r}")
        self._locus = None
        if gene_locus is not None:
            self._locus = self._load_locus(gene_locus, self._row, self.prox_level)
            known = int((self._locus >= 0).sum())
            log.info("codep_seed: proximity exclusion ON at level %r -- %d of %d "
                     "profile genes have a locus (%d distinct)", self.prox_level,
                     known, len(self._locus), len(set(self._locus[self._locus >= 0])))

        # Optional SECOND profile source, scored exactly like the first and fused as
        # its own term. Built for gene-gene co-hit (§4.21): a gene's profile is its
        # hit vector across the 1349 training screens, so two genes are close if they
        # tend to be hits in the same screens -- the gene-gene complement to
        # screen-kNN's screen-screen voting. Any atlas in the shared layout works.
        self._P2 = self._row2 = None
        if profile_npz2 is not None and self.weight2 != 0.0:
            self._P2, self._row2 = self._load_profiles(profile_npz2)[:2]
            log.info("codep_seed: second source %d genes x %d dims from %s, weight %.3g",
                     self._P2.shape[0], self._P2.shape[1],
                     Path(profile_npz2).name, self.weight2)

        # Optional SPARSE interaction kernel (BioPlex AP-MS, §4.28). Physical complex
        # co-membership is a different estimand from co-essentiality: two genes can
        # share a complex with uncorrelated dependency profiles, which is why §4.28
        # found it the second-best kernel. It is kept sparse because the dense form is
        # ~2.5 GB AND because it does not need to be dense -- see _ppi_scores.
        self._A = self._A_row = None
        if ppi_npz is not None and self.ppi_weight != 0.0:
            from scipy.sparse import load_npz
            meta = np.load(Path(ppi_npz), allow_pickle=True)
            self._A = load_npz(Path(ppi_npz).with_suffix(".mat.npz")).tocsr()
            self._A_row = {str(g): i for i, g in enumerate(meta["cols"])}
            self._A_norm = np.asarray(meta["norm"], dtype=np.float64)
            self._A_mean = (np.asarray(meta["row_sum"], dtype=np.float64)
                            / float(meta["n_dims"]))
            self._A_deg = np.diff(self._A.indptr).astype(np.float64)
            log.info("codep_seed: ppi %d genes x %d dims, %d nnz, weight %.3g",
                     self._A.shape[0], self._A.shape[1], self._A.nnz, self.ppi_weight)

        # THE PAPER'S OWN kNN BASELINE, as a fused ranker rather than a replacement
        # model. Standalone it scores EF 1.93 against AssayFormer's 4.83 on paper_test,
        # but it picks an almost disjoint gene set -- 0.040 Jaccard against
        # AssayFormer's 1000 picks, where our method sits at 0.743 -- and 66% of the
        # hits it finds are ones AssayFormer never reaches. Weak and orthogonal is the
        # profile that can still earn a place in a fusion, which is why this exists.
        # Embeddings come from `genept_to_npz.py` with rows ALREADY unit-normalised.
        self.gpt_weight = float(gpt_weight)
        self.gpt_k = int(gpt_k)
        self._G = None
        self._G_row: dict[str, int] = {}
        if gpt_npz is not None and self.gpt_weight != 0.0:
            key = str(Path(gpt_npz).resolve())
            cached = _GPT_CACHE.get(key)
            if cached is None:
                with _GPT_LOCK:
                    cached = _GPT_CACHE.get(key)
                    if cached is None:
                        z = np.load(Path(gpt_npz), allow_pickle=True)
                        cached = (np.asarray(z["values"], dtype=np.float32),
                                  {str(g): i for i, g in enumerate(z["cols"])})
                        _GPT_CACHE[key] = cached
            self._G, self._G_row = cached
            log.info("codep_seed: genept %d genes x %d dims, k=%d, weight %.3g",
                     self._G.shape[0], self._G.shape[1], self.gpt_k, self.gpt_weight)

        # THE PAPER'S OWN Screen-kNN BASELINE (`models/screen_knn.py`, the workbook's
        # `Screen-kNN` row at mean e_t 3.302), fused as a ranker rather than run as a
        # replacement model. Its estimator is NOT ours: distance is squared error over
        # the genes two screens SHARE -- so agreement on misses counts -- exponentially
        # weighted, then shrunk toward a global hit prior by a factor that grows with
        # the number of labels revealed. Ours is hit-overlap only, size-corrected, hard
        # top-50, no prior. Having both lets the project say which estimator the gain
        # belongs to.
        self.pknn_weight = float(pknn_weight)
        self.pknn_train = str(pknn_train)
        self.pknn_tau = float(pknn_tau)
        self.pknn_kappa = float(pknn_kappa)
        # Which distance feeds exp(-dist/tau). `agreement` is THEIRS: squared error over
        # every gene two screens share, so agreeing on a MISS counts -- and because hits
        # are ~5% sparse that is nearly free, which is why the median disagreement rate
        # is 0.33 and four fifths of the vote comes from outside the top 50.
        # Restricting the same quantity to the UNION OF HITS is exactly 1 - Jaccard:
        # over that union a gene is either hit by both (error 0) or by one (error 1), so
        # mean error = (|U| - |inter|)/|U|. `cosine` swaps the union normaliser for
        # sqrt(|H|.|H_s|), which is the similarity our own scorer uses. Holding the
        # exp/tau weighting and the lambda shrinkage fixed and changing only this
        # isolates "which similarity" from "how it is aggregated".
        # Restrict the vote to the k most similar screens, 0 = no truncation (theirs).
        # The decomposition in §4.55 put the ENTIRE difference between the two
        # estimators here: the prior is worth +0.004, switching the similarity to
        # hits-only is worth -0.057, and restricting the neighbourhood is worth +0.204.
        # This exists to test whether that effect transfers to THEIR similarity, which
        # measured slightly better than ours.
        self.pknn_k = int(pknn_k)
        # GATE THE PAPER'S SCORER THE WAY OURS ARE GATED. Default True, which CHANGED
        # the behaviour of this parameter's arms -- deliberately, and §4.56(e) is the
        # reason. Every scorer of ours returns None below `min_seeds` confirmed hits, so
        # round 1 is plain greedy and `weight=0` reproduces AssayFormer pick for pick.
        # `_pknn_scores` had no such check: with zero observations every training screen
        # sits at distance 0, the weights go uniform, `lam = n/(n+kappa) = 0`, and the
        # score is the PURE GLOBAL HIT PRIOR. So it reordered round 1 toward famous genes
        # before the campaign knew anything, its arms were NOT nested in the baseline
        # (round 1 differed on 218 of 218 screens even on matched CPU), and the round-1
        # invariant could not validate them.
        #
        # This is faithful to `ScreenKNNModel` -- a standalone model with no data SHOULD
        # predict the prior -- so it was never a reimplementation bug. It is a
        # FUSION-HARNESS asymmetry, and the harness is ours to choose. Comparing
        # estimators requires both terms silent in the same rounds.
        #
        # `pknn_gate=False` restores the ungated behaviour exactly and is what the
        # pre-2026-10-28 `sweep-pkv-*` arms ran; keep it for reproducing them and as the
        # control that prices the round-1 prior on its own.
        self.pknn_gate = bool(pknn_gate)
        # PER-GENE shrinkage, in units of "equivalent neighbours". Their lambda shrinks
        # every gene by the same amount; the uncertainty is actually per gene, in how
        # many neighbours measured it at all. Their estimator scores a gene measured by
        # ONE neighbour that hit it at exactly 1.0 -- the same as a gene all fifty
        # measured and all fifty hit -- and a ranking punishes that variance. Ours
        # shrinks it the other way, treating "untested by my neighbours" as evidence
        # against. A pseudo-count toward the prior sits between the two and is neither.
        self.pknn_shrink_m = float(pknn_shrink_m)
        self.pknn_sim = str(pknn_sim)
        if self.pknn_sim not in ("agreement", "jaccard", "cosine"):
            raise ValueError("pknn_sim must be agreement|jaccard|cosine, got "
                             f"{pknn_sim!r}")

    @staticmethod
    def _load_profiles(path: str | Path, drop_pcs: int = 0,
                       want_amp: bool = False):
        """Row-centred, unit-normalised gene profiles, so a dot product is a
        correlation. Centring per GENE is what makes this work for a co-hit matrix
        too: it subtracts the gene's marginal hit rate across training screens,
        leaving co-occurrence structure rather than "this gene is often a hit",
        which §10 records as useless on its own.

        ``drop_pcs`` projects out the top principal components of the profile matrix
        BEFORE normalising (§4.44). The 808 columns are not 808 independent
        measurements: a few large factors -- proliferation rate, screen quality,
        lineage -- are shared by every gene, so two genes both loading on "this line
        grows fast" get a high cosine for a reason that is not about either of them.
        Difference-of-means already cancels this in the DIRECTION (which is part of
        why it beats the centroid), but nothing cancels it in the SPACE, and those are
        two separate corrections. Removing shared components is standard in the
        co-essentiality literature. ``drop_pcs=0`` leaves the matrix bit-identical.

        ``want_amp`` also returns each gene's pre-normalisation amplitude, median-
        scaled. Unit-normalising destroys amplitude, so a gene with a nearly FLAT
        profile contributes a pure noise direction that can score arbitrarily high by
        chance -- and §4.25 measured the hard screens' hits at amplitude 0.138 against
        0.137 for random genes, i.e. exactly that population. ``amp_power`` in the
        scorer reweights by it.
        """
        z = np.load(Path(path), allow_pickle=True)
        vals = np.asarray(z["values"], dtype=np.float32)        # (dims, genes)
        cols = [str(c) for c in z["cols"]]
        # posinf/neginf -> 0, not numpy's default 3.4e38: a single infinity in a
        # gene's row makes its norm astronomical and its correlation with every
        # query ~0. DepMap has none, but the Replogle Perturb-seq atlas has ~6.8k
        # (zero-variance genes in its z-scoring), and the failure is silent -- it
        # looks like the feature having no signal.
        P = np.nan_to_num(vals.T.astype(np.float64),
                          posinf=0.0, neginf=0.0)               # (genes, dims)
        P -= P.mean(axis=1, keepdims=True)                      # centre per gene
        if drop_pcs > 0:
            # top right-singular directions of the gene x dim matrix, via the small
            # dim x dim gram matrix -- exact, and 808x808 rather than an SVD of
            # 18119x808
            C = P.T @ P
            w, Vv = np.linalg.eigh(C)
            Vk = Vv[:, np.argsort(-w)[:drop_pcs]]
            P -= (P @ Vk) @ Vk.T
        n = np.linalg.norm(P, axis=1, keepdims=True)
        out = P / np.where(n > 0, n, 1.0), {g: i for i, g in enumerate(cols)}
        if not want_amp:
            return out
        a = n.ravel()
        med = float(np.median(a[a > 0])) or 1.0
        return out[0], out[1], a / med

    def name(self) -> str:
        base = f"codep_seed_{self.query}_w{self.weight:g}_{self.fuse}"
        if self.weight2 != 0 and self._P2 is not None:
            base = f"{base}_p2{self.weight2:g}"
        if self.knn_weight != 0:
            base = f"{base}_knn{self.knn_weight:g}"
        if self.hsa_weight != 0:
            base = f"{base}_hsa{self.hsa_weight:g}"
        if self.af_weight != 1.0:
            base = f"{base}_af{self.af_weight:g}"
        if self.inject_per_round and self.inject_genes:
            base = f"{base}_INJ{self.inject_per_round}x{self.inject_rounds}"
        if self.oracle_seeds:
            base = f"{base}_ORC{len(self.oracle_seeds)}"
        if self.drop_pcs:
            base = f"{base}_pc{self.drop_pcs}"
        if self.amp_power:
            base = f"{base}_amp{self.amp_power:g}"
        if self.ppi_weight != 0 and self._A is not None:
            base = f"{base}_ppi{self.ppi_weight:g}"
            if self.ppi_degree_residual:
                base = f"{base}dr"
        if self.shrink_alpha > 0:
            base = f"{base}_shr{self.shrink_alpha:g}k{self.shrink_k}"
            if self.shrink_drop_top:
                base = f"{base}d{self.shrink_drop_top}"
            if self.shrink_shuffle >= 0:
                base = f"{base}_SHUF"
        return base

    def reset(self) -> None:
        self._last_trace = {}
        self._injected = set()

    def last_trace(self) -> dict[str, Any]:
        return self._last_trace

    # -- internals ---------------------------------------------------------

    def _hit_symbols(self, history: list[StepRecord]) -> list[str]:
        # oracle seeds first, so they survive any later truncation to max_labels
        out = list(self.oracle_seeds)
        for step in history:
            for obs in step.new_observations or []:
                lbl = obs.label
                hit = bool(lbl.get("hit")) if isinstance(lbl, dict) else bool(lbl)
                if hit:
                    out.append(str(obs.candidate))
        return out

    @staticmethod
    def _load_locus(path, row: dict[str, int], level: str) -> np.ndarray:
        """Map each profile row to an integer locus id; -1 where unknown.

        Reads HGNC's `location` column (cytobands such as ``19q13.43``). ``band``
        keys on the whole cytoband, ``arm`` on chromosome plus p/q -- the scale a
        copy-number event actually spans -- and ``chrom`` on the chromosome.
        """
        import csv
        pat = re.compile(r"^([0-9]{1,2}|X|Y|MT)([pq])?([0-9.]*)")
        keys: dict[str, str] = {}
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                sym, loc = r.get("symbol"), (r.get("location") or "").strip()
                if not sym or not loc:
                    continue
                m = pat.match(loc.split()[0])
                if not m or not m.group(2):
                    continue                      # unplaced, mitochondrial, no arm
                chrom, arm, band = m.group(1), m.group(2), m.group(3)
                keys[sym] = {"chrom": chrom, "arm": chrom + arm,
                             "band": chrom + arm + band}[level]
        ids: dict[str, int] = {}
        out = np.full(len(row), -1, dtype=np.int32)
        for sym, i in row.items():
            k = keys.get(sym)
            if k is not None:
                out[i] = ids.setdefault(k, len(ids))
        return out

    def _seed_query(self, history: list[StepRecord], P=None,
                    row=None) -> tuple[np.ndarray | None, int]:
        """Unit query direction from the labels revealed so far.

        ``centroid`` uses only the hits. ``diffmeans`` subtracts the non-hit mean,
        which cancels whatever a typical profile in AssayFormer's picked region looks
        like and leaves only what separates hits from non-hits. Both are closed form:
        online weight FITTING was measured to be far too noisy in this regime (a pure
        confound control scored +0.58 EF under it), so nothing here is fitted.
        """
        # defaults keep the single-source call sites unchanged
        P = self._P if P is None else P
        row = self._row if row is None else row
        H, M, n_hits = self._seed_sets(history, row)
        if H is None:
            return None, n_hits
        q = P[H].mean(axis=0)
        if M:
            q = q - P[M].mean(axis=0)
        nrm = np.linalg.norm(q)
        return (q / nrm if nrm > 0 else None), n_hits

    def _seed_sets(self, history: list[StepRecord], row) -> tuple:
        """``(hit rows, non-hit rows USED, total hits)`` in the profile matrix.

        The single definition of which labels feed the query, so the plain and the
        proximity-corrected scorers cannot drift apart. ``None`` for the rows when
        there are too few hits to point a direction. The non-hit list is empty under
        ``centroid``, and under ``diffmeans`` when too few non-hits are PRESENT IN THE
        ATLAS for a stable mean -- that counts covered genes, not labelled ones.
        """
        hit_rows, miss_rows = [], []
        for step in history:
            for obs in step.new_observations or []:
                lbl = obs.label
                hit = bool(lbl.get("hit")) if isinstance(lbl, dict) else bool(lbl)
                i = row.get(str(obs.candidate))
                if i is None:
                    continue
                (hit_rows if hit else miss_rows).append(i)
        for g in self.oracle_seeds:
            i = row.get(g)
            if i is not None:
                hit_rows.append(i)
        n_hits = len(hit_rows)
        if n_hits < self.min_seeds:
            return None, None, n_hits
        use_miss = (self.query == "diffmeans"
                    and len(miss_rows) >= self.min_nonhits)
        return (hit_rows[:self.max_labels],
                miss_rows[:self.max_labels] if use_miss else [],
                n_hits)

    def _codep_prox_scores(self, history: list[StepRecord], candidates: list):
        """Co-dependency score with chromosomally PROXIMAL seeds excluded, per gene.

        Exact, not approximate. The mean over seeds is a sum, so dropping the seeds
        that share a candidate's locus is the global seed sum minus that locus's seed
        sum -- one matrix-vector product per distinct locus among the candidates
        rather than one per gene, so the whole thing costs about what the uncorrected
        single product costs.

        A gene whose locus holds EVERY seed scores nan and is treated as uncovered,
        which the fusion already handles by giving it the mean covered contribution.
        """
        H, M, n_hits = self._seed_sets(history, self._row)
        if H is None:
            return None, n_hits, 0
        P, loc = self._P, self._locus
        S_H, n_H = P[H].sum(axis=0), len(H)
        S_M, n_M = (P[M].sum(axis=0), len(M)) if M else (None, 0)
        hit_by: dict[int, list[int]] = {}
        miss_by: dict[int, list[int]] = {}
        for r in H:
            L = int(loc[r])
            if L >= 0:
                hit_by.setdefault(L, []).append(r)
        for r in M:
            L = int(loc[r])
            if L >= 0:
                miss_by.setdefault(L, []).append(r)
        sumH = {L: P[rs].sum(axis=0) for L, rs in hit_by.items()}
        sumM = {L: P[rs].sum(axis=0) for L, rs in miss_by.items()}

        rows = np.array([self._row.get(str(c), -1) for c in candidates])
        out = np.full(len(candidates), np.nan)
        # -2 marks a gene absent from the atlas; -1 a gene with no known locus, which
        # is never excluded from anything
        L_of = np.where(rows >= 0, loc[np.clip(rows, 0, None)], -2)
        n_corrected = 0
        for L in np.unique(L_of):
            L = int(L)
            if L == -2:
                continue
            sel = np.where(L_of == L)[0]
            rr = rows[sel]
            sH, dH = S_H, n_H
            if L >= 0 and L in hit_by:
                sH, dH = S_H - sumH[L], n_H - len(hit_by[L])
                n_corrected += len(sel)
            if dH <= 0:
                continue
            v = (P[rr] @ sH) / dH
            if n_M:
                sM, dM = S_M, n_M
                if L >= 0 and L in miss_by:
                    sM, dM = S_M - sumM[L], n_M - len(miss_by[L])
                if dM > 0:
                    v = v - (P[rr] @ sM) / dM
            out[sel] = v
        return out, n_hits, n_corrected

    def _neighbor_sim(self, hit_symbols: list[str]) -> np.ndarray | None:
        """Similarity of every training screen to this one: how many of the hits
        revealed here are also hits there, corrected for the training screen's own
        size so a 7000-hit screen is not everyone's nearest neighbour. Shared by the
        screen-kNN scorer and by shrinkage, which must agree on "similar"."""
        if self._K_csc is None:
            return None
        kidx = [self._K_col[g] for g in hit_symbols if g in self._K_col]
        if len(kidx) < self.min_seeds:
            return None
        overlap = np.asarray(self._K_csc[:, kidx].sum(axis=1)).ravel()
        return overlap / self._K_size

    def _apply_injection(self, order: list[int], candidates: list,
                         history: list[StepRecord], k: int) -> tuple[list[int], int]:
        """Give the m lowest-ranked slots of this batch to genes the ranker would not
        have chosen (§4.50).

        Displacing the WEAKEST picks is what a mixed-batch acquisition does
        (`--mix-start`, §4) and is the conservative choice: displacing the strongest
        would overstate the displacement cost. The batch stays exactly k picks, so EF's
        denominator is untouched and the 2x2 cancels.

        MUST be applied on EVERY return path, including the weight-0 nested baseline --
        that path is precisely the kernel-OFF arm of the 2x2, and skipping it there
        would compare an injected arm against a non-injected one.
        """
        picked = order[:k]
        if not (self.inject_per_round and self.inject_genes):
            return picked, 0
        if len(history) >= self.inject_rounds:
            return picked, 0
        present = {str(c): i for i, c in enumerate(candidates)}
        avail = [g for g in self.inject_genes
                 if g not in self._injected and g in present]
        m = min(self.inject_per_round, len(avail), k)
        if not m:
            return picked, 0
        chosen = avail[:m]
        self._injected.update(chosen)
        return order[:k - m] + [present[g] for g in chosen], m

    def _ppi_scores(self, history: list[StepRecord],
                    candidates: list) -> np.ndarray | None:
        """Interaction-kernel score per candidate, computed WITHOUT densifying.

        The dense form of this atlas is ~2.5 GB, and §4.28 rejected the pair partly on
        that cost. It is unnecessary. The kernel wants
        ``cos(q, centred row of g)``, and for the centred rows

          * the NUMERATOR loses its centring term exactly. ``q`` is a combination of
            centred rows, each of which sums to zero, so ``sum(q) == 0`` and
            ``(A[g] - m_g * 1) . q == A[g] . q - m_g * (1 . q) == A[g] . q``.
          * the DENOMINATOR is precomputed by `bioplex_to_npz.py` as
            ``sqrt(ssq_g - sum_g**2 / d)``.

        So every candidate is scored by one sparse matvec. This is the exact same
        quantity the dense path would produce, not an approximation of it.
        """
        if self._A is None:
            return None
        hit_rows, miss_rows = [], []
        for step in history:
            for obs in step.new_observations or []:
                lbl = obs.label
                hit = bool(lbl.get("hit")) if isinstance(lbl, dict) else bool(lbl)
                i = self._A_row.get(str(obs.candidate))
                if i is None:
                    continue
                (hit_rows if hit else miss_rows).append(i)
        for g in self.oracle_seeds:
            i = self._A_row.get(g)
            if i is not None:
                hit_rows.append(i)
        if len(hit_rows) < self.min_seeds:
            return None

        def mean_centred(rows: list[int]) -> np.ndarray:
            """Mean of the centred, unit-normalised rows, as a dense vector. The
            ``-m_i`` offsets collect into a single scalar times the all-ones vector,
            so nothing sparse is expanded."""
            r = np.asarray(rows[:self.max_labels], dtype=np.int64)
            inv = 1.0 / np.where(self._A_norm[r] > 0, self._A_norm[r], 1.0)
            v = np.asarray(self._A[r].multiply(inv[:, None]).sum(axis=0)).ravel()
            return v / len(r) - float((self._A_mean[r] * inv).sum() / len(r))

        q = mean_centred(hit_rows)
        if self.query == "diffmeans" and len(miss_rows) >= self.min_nonhits:
            q = q - mean_centred(miss_rows)
        nrm = np.linalg.norm(q)
        if nrm <= 0:
            return None
        q = q / nrm

        have = np.array([self._A_row.get(str(c), -1) for c in candidates])
        ok = have >= 0
        out = np.full(len(candidates), np.nan)
        if not ok.any():
            return None
        rows = have[ok]
        num = np.asarray(self._A[rows] @ q).ravel()
        out[ok] = num / np.where(self._A_norm[rows] > 0, self._A_norm[rows], 1.0)
        if self.ppi_degree_residual:
            # CONTROL for §4.28's partial gate-1 failure: interaction degree
            # correlates +0.203 with cross-screen hit frequency, so part of this
            # kernel's value may be the §4.11 fame redundancy rather than new
            # information. Residualising on log degree removes the hub component
            # from the ranking; if the gain survives, it was not fame.
            x = np.log1p(self._A_deg[rows])
            X = np.column_stack([np.ones(len(x)), x])
            beta, *_ = np.linalg.lstsq(X, out[ok], rcond=None)
            out[ok] = out[ok] - X @ beta
        return out

    def _prior_dir(self, sim: np.ndarray | None) -> tuple[np.ndarray | None, int]:
        """Similarity-weighted mean of the nearest training screens' own optimal
        directions, as a unit vector.

        ``shrink_drop_top`` discards the most similar screens before averaging. The
        evaluation folds are disjoint from the training fold, but near-duplicate
        screens recur across this benchmark (HANDOFF §6 finds hit-set Jaccard up to
        1.00 WITHIN a curated set), so a neighbour can be effectively the same assay
        and the borrowed direction then approaches the live screen's own answer key.
        Dropping the top few is how that is tested rather than assumed away."""
        if self._bank is None or sim is None:
            return None, 0
        usable = np.linalg.norm(self._bank, axis=1) > 0
        s = np.where(usable, sim, -np.inf)
        order = np.argsort(-s)[self.shrink_drop_top:self.shrink_drop_top + self.shrink_k]
        top = [int(i) for i in order if np.isfinite(s[i]) and s[i] > 0]
        if not top:
            return None, 0
        w = sim[top]
        p = self._bank[top].T @ w / w.sum()
        nrm = np.linalg.norm(p)
        return (p / nrm if nrm > 0 else None), len(top)

    def _pknn_state(self, genes: tuple):
        """Cached (hit matrix, measured matrix, prior) over `genes`, built once per
        PROCESS rather than per screen: it is ~230 MB and identical for every screen
        in a sweep."""
        key = (self.pknn_train, genes)
        st = _PKNN_CACHE.get(key)
        if st is not None:
            return st
        with _PKNN_LOCK:
            st = _PKNN_CACHE.get(key)          # another worker may have won the race
            if st is None:
                from ..scripts.plot_lopo_results import _build_screen_knn_matrices
                from ..tasks import load_screens
                train = load_screens(target_set=self.pknn_train)
                built = _build_screen_knn_matrices(train, list(genes))
                names = [s.dataset_name for s in train]
                if len(_PKNN_CACHE) >= _PKNN_CACHE_MAX:
                    log.warning(
                        "pknn state cache at %d entries: the gene universe is varying "
                        "between calls when it should be fixed for the whole sweep. "
                        "Evicting; expect a rebuild per screen and check that "
                        "_pknn_universe is reconstructed from candidates + history.",
                        len(_PKNN_CACHE))
                    _PKNN_CACHE.pop(next(iter(_PKNN_CACHE)))
                _PKNN_CACHE[key] = st = (*built, names)
                log.info("codep_seed: paper Screen-kNN state %s over %d training "
                         "screens", built[0].shape, len(names))
        return st

    def _pknn_scores(self, history: list[StepRecord], candidates: list,
                     task_context: dict | None = None) -> np.ndarray | None:
        """The paper's Screen-kNN, vectorised. Reproduces `ScreenKNNModel` exactly;
        its own pure-Python loop is 21k candidates x 1349 screens per round, which is
        hours per sweep, and `plot_lopo_results` records this form as ~800x faster
        with a diff of 0."""
        if self._pknn_universe is None:
            # THE UNIVERSE MUST NOT DEPEND ON WHEN THIS SCORER FIRST RUNS.
            # `candidates` is the UNTESTED pool and shrinks by `batch_size` every
            # round, so pinning it at the first call was only ever safe because the
            # first call was round 1, where it equals the full pool for every screen.
            # `pknn_gate` delays that call to round 2-4, and the acquired genes differ
            # per screen -- which gave every screen its own cache key and its own
            # ~230 MB state (OOM after ~32 of 218 screens at 7.4 GB), AND built the
            # global prior over "whatever was still untested when the gate opened"
            # rather than over a fixed universe. Adding the already-acquired genes
            # back reconstructs the full shared pool, identically at every round.
            seen = {str(g) for st in history for g in (st.acquired_batch or [])}
            seen |= {str(o.candidate) for st in history
                     for o in (st.new_observations or [])}
            self._pknn_universe = tuple(sorted({str(c) for c in candidates} | seen))
            self._pknn_at = {g: i for i, g in enumerate(self._pknn_universe)}
        H, M, prior, names = self._pknn_state(self._pknn_universe)

        rev_idx, rev_val = [], []
        for step in history:
            for obs in step.new_observations or []:
                i = self._pknn_at.get(str(obs.candidate))
                if i is None:
                    continue
                lbl = obs.label
                hit = bool(lbl.get("hit")) if isinstance(lbl, dict) else bool(lbl)
                rev_idx.append(i); rev_val.append(1.0 if hit else 0.0)

        nrev = len(rev_idx)
        # kappa = 0 makes lambda 1, i.e. no shrinkage toward the global hit prior --
        # the switch that separates "the prior is the problem" from "the similarity is"
        lam = nrev / (nrev + self.pknn_kappa) if nrev else 0.0
        flat = np.exp(-1.0 / self.pknn_tau)
        R = np.asarray(rev_idx); r = np.asarray(rev_val, np.float32)
        hit_cols = R[r > 0] if nrev else R
        if nrev == 0 or (self.pknn_sim != "agreement" and len(hit_cols) == 0):
            w = np.full(H.shape[0], flat, np.float32)
        elif self.pknn_sim == "agreement":
            MR = M[:, R]
            cnt = MR.sum(1)
            sq = ((r[None, :] - H[:, R]) ** 2 * MR).sum(1)
            with np.errstate(divide="ignore", invalid="ignore"):
                dist = np.where(cnt > 0, sq / np.maximum(cnt, 1), np.nan)
            w = np.where(cnt > 0, np.exp(-dist / self.pknn_tau),
                         flat).astype(np.float32)
        else:
            inter = H[:, hit_cols].sum(1)                     # |H n H_s|
            n_s = H.sum(1)                                    # |H_s|
            n_h = float(len(hit_cols))                        # |H|
            if self.pknn_sim == "jaccard":
                denom = n_h + n_s - inter
                sim = np.where(denom > 0, inter / np.maximum(denom, 1e-12), 0.0)
            else:                                             # cosine
                denom = np.sqrt(np.maximum(n_s, 1.0) * max(n_h, 1.0))
                sim = inter / denom
            w = np.exp(-(1.0 - sim) / self.pknn_tau).astype(np.float32)
        if self.pknn_k > 0 and self.pknn_k < len(w):
            # zero everything outside the top k rather than reweighting, so the
            # surviving weights keep exactly the values their estimator gives them
            cut = np.argpartition(-w, self.pknn_k - 1)[self.pknn_k:]
            w = w.copy()
            w[cut] = 0.0
        # the screen under test must not vote on itself
        cur = str((task_context or {}).get("dataset_name", ""))
        if cur:
            w = w * np.array([0.0 if n == cur else 1.0 for n in names], np.float32)

        wM, wH = w @ M, w @ H
        with np.errstate(divide="ignore", invalid="ignore"):
            if self.pknn_shrink_m > 0:
                # m is scaled by the MEAN surviving weight, so the pseudo-count means
                # the same thing whatever the weight function's scale -- exp(-dist/tau)
                # runs around 1e-3 here, so an unscaled count of 1 would swamp it
                pos = w > 0
                wbar = float(w[pos].mean()) if pos.any() else 0.0
                m = self.pknn_shrink_m * wbar
                evid = (wH + m * prior) / np.maximum(wM + m, 1e-12)
                # a gene no neighbour measured now falls to the prior on its own
                full = (1 - lam) * prior + lam * evid
            else:
                evid = np.where(wM > 0, wH / np.maximum(wM, 1e-12), 0.0)
                full = np.where(wM > 0, (1 - lam) * prior + lam * evid, prior)
        out = np.full(len(candidates), np.nan)
        for k, c in enumerate(candidates):
            i = self._pknn_at.get(str(c))
            if i is not None:
                out[k] = full[i]
        return out

    def _gpt_scores(self, history: list[StepRecord],
                    candidates: list) -> np.ndarray | None:
        """The paper's kNN baseline, scored as a ranker.

        Faithful to `KNNGeneEmbedding`: a distance-weighted vote over the k nearest
        ALREADY-TESTED genes in GenePT embedding space, weight ``1/(1 - cos + eps)``,
        labels in {0, 1}. Unlike every other scorer here it needs NO hits -- a round
        of pure misses is still information to it -- so its gate is one label, not
        ``min_seeds`` confirmed hits.
        """
        if self._G is None:
            return None
        rows, y = [], []
        for step in history:
            for obs in step.new_observations or []:
                i = self._G_row.get(str(obs.candidate))
                if i is None:
                    continue
                lbl = obs.label
                rows.append(i)
                y.append(1.0 if (bool(lbl.get("hit")) if isinstance(lbl, dict)
                                 else bool(lbl)) else 0.0)
        if not rows:
            return None
        # NO max_labels truncation here. That cap exists for the DepMap scorer, where
        # it bounds the cost of averaging profiles; for a kNN it would silently drop
        # the earliest rounds from the neighbour pool and stop this being THEIR
        # estimator, which is the only reason to run it.
        rows = np.asarray(rows, dtype=np.int64)
        y = np.asarray(y, dtype=np.float64)
        k = min(self.gpt_k, len(rows))

        have = np.array([self._G_row.get(str(c), -1) for c in candidates])
        ok = have >= 0
        out = np.full(len(candidates), np.nan)
        if not ok.any():
            return None
        # rows are unit-norm, so the product IS the cosine matrix (n_cand x n_tested)
        S = self._G[have[ok]] @ self._G[rows].T
        top = np.argpartition(-S, k - 1, axis=1)[:, :k]
        sim = np.take_along_axis(S, top, axis=1).astype(np.float64)
        w = 1.0 / (1.0 - sim + 1e-6)
        ws = w.sum(axis=1)
        out[ok] = (w * y[top]).sum(axis=1) / np.where(ws > 0, ws, 1.0)
        return out

    def _knn_scores(self, sim: np.ndarray | None, candidates: list[Any]) -> np.ndarray | None:
        """Hit frequency among the training screens most similar to this one."""
        if sim is None:
            return None
        top = np.argsort(-sim)[:self.knn_k]
        w = sim[top]
        if w.sum() <= 0:
            return None
        prof = np.asarray(self._K_csr[top].T.dot(w)).ravel() / w.sum()
        out = np.full(len(candidates), np.nan)
        for i, c in enumerate(candidates):
            j = self._K_col.get(str(c))
            if j is not None:
                out[i] = prof[j]
        return out

    def _hsa_scores(self, history: list[StepRecord], candidates: list[Any],
                    task_context: dict[str, Any] | None) -> np.ndarray | None:
        """Per-candidate score from the caller-supplied 42-D historical scorer.

        The flat observation list is what the reference scorer expects -- it selects
        the measured ones itself via ``metadata['in_library']``, so filtering here
        would change its fit. Returns ``None`` when no scorer is configured, which
        leaves the ensemble exactly as it was.
        """
        if self._hsa_scorer is None:
            return None
        obs = [o for st in history for o in (st.new_observations or [])]
        v = self._hsa_scorer(obs, list(candidates), task_context)
        if isinstance(v, dict):
            out = np.array([v[c] for c in candidates], dtype=np.float64)
        else:
            out = np.asarray(v, dtype=np.float64)
            if out.shape != (len(candidates),):
                raise ValueError(
                    f"hsa_scorer returned {out.shape}, expected ({len(candidates)},)")
        if not np.isfinite(out).all():
            raise ValueError("hsa_scorer returned a non-finite score")
        return out

    @staticmethod
    def _ranks(order: list[int], n: int) -> np.ndarray:
        r = np.full(n, np.inf)
        for pos, i in enumerate(order):
            r[i] = pos
        return r

    def suggest(
        self,
        history: list[StepRecord],
        candidates: list[Any],
        batch_size: int,
        model_prediction: ModelPrediction | None = None,
        task_context: dict[str, Any] | None = None,
    ) -> list[Any]:
        k = min(batch_size, len(candidates))
        if not candidates:
            return []
        scores = (model_prediction.scores if model_prediction else None) or {}
        if not scores:
            self._last_trace = {"reason": "no_prediction"}
            rng = np.random.default_rng(self.seed)
            return [candidates[i] for i in rng.permutation(len(candidates))[:k]]

        if self.oracle_seeds:
            # THE ESSENTIAL HALF OF THE CONTROL: a gene handed over as a known hit must
            # not then be pickable, or the arm banks it as a hit and the "bound" is just
            # the free hits. Choose oracle seeds that the comparison arms never acquired
            # (see the §4.49 driver) and this costs the oracle arm nothing relative to
            # them.
            blocked = set(self.oracle_seeds)
            candidates = [c for c in candidates if str(c) not in blocked]
        af = np.array([scores.get(c, -np.inf) for c in candidates], dtype=np.float64)
        af_order = sorted(range(len(candidates)), key=lambda i: (-af[i], str(candidates[i])))

        q, n_seeds = self._seed_query(history)
        sim = (self._neighbor_sim(self._hit_symbols(history))
               if (self.knn_weight != 0.0 or self.shrink_alpha > 0.0) else None)
        knn = self._knn_scores(sim, candidates) if self.knn_weight != 0.0 else None
        ppi = (self._ppi_scores(history, candidates)
               if (self.ppi_weight != 0.0 and self._A is not None) else None)
        gpt = (self._gpt_scores(history, candidates)
               if (self.gpt_weight != 0.0 and self._G is not None) else None)
        # `self._hit_symbols(history)` is already materialised above for our own kNN.
        pknn = (self._pknn_scores(history, candidates, task_context)
                if (self.pknn_weight != 0.0
                    and (not self.pknn_gate
                         or len(self._hit_symbols(history)) >= self.min_seeds))
                else None)
        # Exactly once per round, and only when weighted: the scorer is allowed to
        # be stateful about the round index.
        hsa = (self._hsa_scores(history, candidates, task_context)
               if self.hsa_weight != 0.0 else None)

        # shrink the live direction toward the prior borrowed from similar
        # training screens. alpha=0 leaves q untouched, bit for bit.
        n_prior, cos_qp = 0, float("nan")
        if q is not None and self.shrink_alpha > 0.0:
            p, n_prior = self._prior_dir(sim)
            if p is not None:
                cos_qp = float(q @ p)
                blend = (1.0 - self.shrink_alpha) * q + self.shrink_alpha * p
                nrm = np.linalg.norm(blend)
                if nrm > 0:
                    q = blend / nrm

        use_codep = q is not None and self.weight != 0.0
        q2, use_codep2 = None, False
        if self._P2 is not None:
            q2, _n2 = self._seed_query(history, self._P2, self._row2)
            use_codep2 = q2 is not None and self.weight2 != 0.0
        if (not use_codep and not use_codep2 and knn is None and ppi is None
                and gpt is None and pknn is None and hsa is None):
            # nested baseline: exactly greedy
            picked, n_inj = self._apply_injection(af_order, candidates, history, k)
            self._last_trace = {"n_seeds": n_seeds, "codep_active": False,
                                "codep2_active": False, "knn_active": False,
                                "ppi_active": False, "gpt_active": False,
                                "pknn_active": False,
                                "n_injected": n_inj,
                                "inject_per_round": self.inject_per_round,
                                "reason": "no_seeds" if q is None else "weight_zero"}
            return [candidates[i] for i in picked]

        cd = np.full(len(candidates), np.nan)
        covered = np.zeros(len(candidates), dtype=bool)
        prox_corrected = 0
        if use_codep:
            have = np.array([self._row.get(str(c), -1) for c in candidates])
            covered = have >= 0
            if self._locus is None:
                cd[covered] = self._P[have[covered]] @ q
            else:
                cdp, _nh, prox_corrected = self._codep_prox_scores(history, candidates)
                if cdp is None:
                    covered = np.zeros(len(candidates), dtype=bool)
                else:
                    cd, covered = cdp, np.isfinite(cdp)
            if self.amp_power:
                # reliability weight: a flat profile's cosine is noise whatever its
                # value, so scale it down rather than letting it compete at the
                # top-100 cut where picks actually change
                cd[covered] *= self._amp[have[covered]] ** self.amp_power
        cd2 = np.full(len(candidates), np.nan)
        if use_codep2:
            have2 = np.array([self._row2.get(str(c), -1) for c in candidates])
            ok2 = have2 >= 0
            cd2[ok2] = self._P2[have2[ok2]] @ q2

        def _z(v):
            """Standardise the defined entries; undefined ones take the mean (zero),
            so that merely HAVING a score is not itself worth anything."""
            out = np.zeros(len(candidates))
            ok = np.isfinite(v)
            if ok.any():
                out[ok] = (v[ok] - v[ok].mean()) / (v[ok].std() or 1.0)
            return out

        def _rrf_term(v):
            """Reciprocal-rank contribution, with undefined entries at the mean of
            the defined ones -- the same neutrality rule as _z."""
            ok = np.isfinite(v)
            if not ok.any():
                return np.zeros(len(candidates))
            order = [int(i) for i in np.argsort(-np.where(ok, v, -np.inf)) if ok[i]]
            r = self._ranks(order, len(candidates))
            t = np.where(np.isfinite(r), 1.0 / (RRF_K + r), np.nan)
            t[~np.isfinite(t)] = np.nanmean(t)
            return t

        if self.fuse == "zsum":
            a = (af - af[np.isfinite(af)].mean()) / (af[np.isfinite(af)].std() or 1.0)
            fused = self.af_weight * a + self.weight * _z(cd)
            if use_codep2:
                fused = fused + self.weight2 * _z(cd2)
            if knn is not None:
                fused = fused + self.knn_weight * _z(knn)
            if ppi is not None:
                fused = fused + self.ppi_weight * _z(ppi)
            if gpt is not None:
                fused = fused + self.gpt_weight * _z(gpt)
            if pknn is not None:
                fused = fused + self.pknn_weight * _z(pknn)
            if hsa is not None:
                fused = fused + self.hsa_weight * _z(hsa)
            order = sorted(range(len(candidates)),
                           key=lambda i: (-fused[i], str(candidates[i])))
        else:
            # RRF. Uncovered genes take the MEAN covered contribution, so that
            # being in the atlas is not itself worth up to weight/RRF_K.
            r_af = self._ranks(af_order, len(candidates))
            fused = self.af_weight * (1.0 / (RRF_K + r_af))
            if use_codep:
                fused = fused + self.weight * _rrf_term(cd)
            if use_codep2:
                fused = fused + self.weight2 * _rrf_term(cd2)
            if knn is not None:
                fused = fused + self.knn_weight * _rrf_term(knn)
            if ppi is not None:
                fused = fused + self.ppi_weight * _rrf_term(ppi)
            if gpt is not None:
                fused = fused + self.gpt_weight * _rrf_term(gpt)
            if pknn is not None:
                fused = fused + self.pknn_weight * _rrf_term(pknn)
            if hsa is not None:
                fused = fused + self.hsa_weight * _rrf_term(hsa)
            order = sorted(range(len(candidates)),
                           key=lambda i: (-fused[i], str(candidates[i])))

        picked, n_inj = self._apply_injection(order, candidates, history, k)
        self._last_trace = {
            "n_seeds": n_seeds, "codep_active": bool(use_codep), "fuse": self.fuse,
            "weight": self.weight, "query": self.query,
            "af_weight": self.af_weight,
            "hsa_active": hsa is not None, "hsa_weight": self.hsa_weight,
            "knn_active": knn is not None, "knn_weight": self.knn_weight,
            "n_covered": int(covered.sum()), "n_candidates": len(candidates),
            "prox_level": self.prox_level if self._locus is not None else None,
            "prox_corrected": int(prox_corrected),
            "overlap_with_greedy": len(set(picked) & set(af_order[:k])),
            "codep2_active": bool(use_codep2), "weight2": self.weight2,
            "shrink_alpha": self.shrink_alpha, "n_prior_screens": n_prior,
            "n_injected": n_inj, "inject_per_round": self.inject_per_round,
            "drop_pcs": self.drop_pcs, "amp_power": self.amp_power,
            "ppi_active": ppi is not None, "ppi_weight": self.ppi_weight,
            "gpt_active": gpt is not None, "gpt_weight": self.gpt_weight,
            "pknn_active": pknn is not None, "pknn_weight": self.pknn_weight,
            "n_ppi_covered": int(np.isfinite(ppi).sum()) if ppi is not None else 0,
            # cosine between the live estimate and the borrowed prior: the
            # diagnostic that says whether shrinkage is nudging or overruling
            "cos_q_prior": cos_qp,
        }
        return [candidates[i] for i in picked]


__all__ = ["CodepSeedAcquisition"]
