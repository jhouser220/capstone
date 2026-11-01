"""Inference-time Model wrapper for the amortized gene ranker.

Loads a trained :class:`~assayloop.amortized.model.RankerNet` (from a
checkpoint directory) and exposes it through the assayloop
:class:`~assaybench.core.Model` interface so it can be driven by any
acquisition (greedy/UCB) inside the AL loop and evaluated like every other
run type.

It can also be built from in-memory components (used by the trainer to
evaluate the live model on the validation screens without round-tripping to
disk).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import torch

from assaybench.core.model import Model
from assaybench.core.types import ModelPrediction, Observation

log = logging.getLogger("assayloop.models.amortized_ranker")


def _pick_device(device: str | None) -> torch.device:
    if device and device != "auto":
        return torch.device(device)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _resolve_ckpt(ckpt: Path, ckpt_file: str) -> Path:
    """The requested weights file, or an error that says what is there instead.

    Deliberately not a fallback: a directory holding only ``model_last.pt`` is
    not a directory holding ``model.pt``, and loading the final epoch when the
    caller asked for the best-on-validation one would change a reported number
    without saying so. But the two names differ by a suffix and mean different
    things, and the released checkpoint on HuggingFace ships only the second, so
    the bare ``FileNotFoundError`` torch raises is a bad first experience. Name
    the alternatives and let the caller choose.
    """
    path = ckpt / ckpt_file
    if path.exists():
        return path
    others = sorted(p.name for p in ckpt.glob("*.pt"))
    raise FileNotFoundError(
        f"{path} does not exist. " + (
            f"This checkpoint directory holds {', '.join(others)} -- pass "
            f"ckpt_file='{others[0]}' if that is the one you want. model.pt is "
            "the best-on-validation epoch and model_last.pt the final one; the "
            "paper's AssayLoop rows, and the released Genentech/assayformer "
            "checkpoint, use model_last.pt."
            if others else
            f"There are no .pt files in {ckpt} at all."
        )
    )


class AmortizedRankerModel(Model):
    """Score candidates with a trained context-conditioned ranker.

    Construct either from a checkpoint dir (``checkpoint=...``) or from
    in-memory components (``net``, ``vocab``, ``embedder``).

    Args:
        checkpoint: path to a directory containing ``model.pt``, ``vocab.json``,
            ``config.json`` (as written by ``amortized.train``).
        ckpt_file: which weights file in that directory to load. ``model.pt`` is
            the best-on-validation epoch, ``model_last.pt`` the final one; the
            paper's AssayLoop rows and the released ``Genentech/assayformer``
            checkpoint use ``model_last.pt``, which has to be asked for.
        net / vocab / embedder: in-memory alternative to ``checkpoint``.
        desc_emb_by_name: optional ``{dataset_name: vector}`` to avoid
            recomputing description embeddings during evaluation.
        device: "auto" | "cpu" | "cuda".

        ctx_ensemble: number of context variants to encode per prediction (K).
            The screen latent û is a deterministic function of the observed-gene
            *set*, so bagging that set gives K different û for the same history;
            their rankings are then fused. ``1`` (the default) is the released
            model's exact behaviour -- the baseline is nested, so the worst case
            is parity. Costs K forward passes and no API calls. The checkpoint is
            ``variational: false`` with ``dropout: 0.0``, so perturbing the
            context is the only available source of û spread.
        ctx_subsample: fraction of observations kept in each variant when
            ``ctx_ensemble > 1``. The first variant always keeps the full
            context, so the ensemble contains the unmodified model.
        ctx_fuse: "mean" (average the scores) or "rrf" (reciprocal rank fusion
            over the K rankings, k=60).
        ctx_weak_hits: what to do with revealed hits whose ``relevance_score``
            is in the weakest ``ctx_weak_q`` quantile of this screen's revealed
            hits -- "keep" (default, unmodified), "drop" (remove the observation
            from the context) or "negative" (keep the gene but flip its label).
            ``relevance_score`` is revealed per acquisition by the task
            (gene_batch.py), so this uses no privileged information. Motivated by
            the measurement that the weak tail of the hit labels is close to
            noise on annotation coverage.
        ctx_weak_q: the quantile boundary for the above (0.25 = weakest quarter).
        ctx_seed: seed for the subsampling RNG. Variants are drawn from a
            generator seeded by (ctx_seed, number of observations), so a given
            round of a given screen is reproducible.
    """

    def __init__(
        self,
        *,
        checkpoint: str | Path | None = None,
        ckpt_file: str = "model.pt",
        net: Any = None,
        vocab: Any = None,
        embedder: Any = None,
        desc_emb_by_name: dict[str, np.ndarray] | None = None,
        desc_text_by_name: dict[str, str] | None = None,
        tokenizer: Any = None,
        ignore_context: bool = False,
        device: str | None = "auto",
        ctx_ensemble: int = 1,
        ctx_subsample: float = 0.8,
        ctx_fuse: str = "mean",
        ctx_weak_hits: str = "keep",
        ctx_weak_q: float = 0.25,
        ctx_seed: int = 0,
        score_mode: str = "full",
    ):
        from ..amortized.data import GeneVocab
        from ..amortized.model import RankerConfig, RankerNet

        self.device = _pick_device(device)
        self._desc_emb_by_name = desc_emb_by_name or {}
        self._desc_text_by_name = desc_text_by_name or {}
        self._emb_cache = None
        self._meta: dict[str, Any] = {}
        self.tokenizer = tokenizer
        # Diagnostic ablation: when True, the model ignores the accumulated AL
        # observations and always scores from the description alone (a static
        # one-shot ranking). Lets eval-ranker measure the marginal value of
        # the context.
        self.ignore_context = bool(ignore_context)

        # -- test-time context interventions (no retraining; K=1 + "keep" is
        # -- exactly the released model)
        self.ctx_ensemble = max(1, int(ctx_ensemble))
        self.ctx_subsample = float(ctx_subsample)
        self.ctx_fuse = str(ctx_fuse)
        self.ctx_weak_hits = str(ctx_weak_hits)
        self.ctx_weak_q = float(ctx_weak_q)
        self.ctx_seed = int(ctx_seed)
        # Ablation ladder over the two additive terms of a(g) = u.v_g + b_g:
        #   full       both terms (the released model)
        #   no_bias    drop b_g -- the purely screen-conditioned part
        #   bias_only  drop u.v_g -- a fixed global gene ranking that ignores
        #              the screen entirely, the floor any policy should beat
        self.score_mode = str(score_mode)
        if self.score_mode not in ("full", "no_bias", "bias_only"):
            raise ValueError(f"score_mode must be full|no_bias|bias_only, got {score_mode!r}")
        if self.ctx_fuse not in ("mean", "rrf"):
            raise ValueError(f"ctx_fuse must be 'mean' or 'rrf', got {ctx_fuse!r}")
        if self.ctx_weak_hits not in ("keep", "drop", "negative"):
            raise ValueError(
                f"ctx_weak_hits must be 'keep', 'drop' or 'negative', got {ctx_weak_hits!r}")
        if not 0.0 < self.ctx_subsample <= 1.0:
            raise ValueError(f"ctx_subsample must be in (0, 1], got {ctx_subsample}")

        if checkpoint is not None:
            ckpt = Path(checkpoint)
            cfg_d = json.loads((ckpt / "config.json").read_text())
            arch = cfg_d["arch"]
            self.vocab = GeneVocab.load(ckpt / "vocab.json")
            self.net = RankerNet(RankerConfig(**arch))
            state = torch.load(_resolve_ckpt(ckpt, ckpt_file),
                               map_location=self.device)
            self.net.load_state_dict(state)
            self._text_backend = cfg_d.get("text_backend", "auto")
            self._text_model = cfg_d.get("text_model")
            self.embedder = embedder
            self._ckpt = str(ckpt)
        else:
            if net is None or vocab is None:
                raise ValueError("Provide either checkpoint= or (net=, vocab=).")
            self.net = net
            self.vocab = vocab
            self.embedder = embedder
            self._text_backend = "auto"
            self._text_model = None
            self._ckpt = ""

        self.encoder_type = getattr(self.net.cfg, "encoder_type", "transformer")
        # The text encoder needs a tokenizer; load it lazily from the saved
        # backbone name if one was not supplied by the trainer.
        if self.encoder_type == "modernbert_text" and self.tokenizer is None:
            from transformers import AutoTokenizer

            self.tokenizer = AutoTokenizer.from_pretrained(self.net.cfg.bert_model)

        self.net.to(self.device)
        self.net.eval()

    # -- lazy text embedder (only needed when desc emb not precomputed) -----

    def _ensure_embedder(self):
        if self.embedder is None:
            from ..amortized import text_embed as te

            kwargs = {}
            if self._text_model:
                kwargs["model"] = self._text_model
            self.embedder = te.get_text_embedder(self._text_backend, **kwargs)
        return self.embedder

    def _desc_emb(self, ctx: dict[str, Any]) -> np.ndarray:
        name = ctx.get("dataset_name")
        if name and name in self._desc_emb_by_name:
            return np.asarray(self._desc_emb_by_name[name], dtype=np.float32)
        from ..amortized import text_embed as te

        if self._emb_cache is None:
            self._emb_cache = te.EmbeddingCache()
        text = te.screen_description_text(ctx)
        emb = te.embed_texts([text], self._ensure_embedder(), cache=self._emb_cache)
        return emb[0]

    def _desc_text(self, ctx: dict[str, Any]) -> str:
        name = ctx.get("dataset_name")
        if name and name in self._desc_text_by_name:
            return self._desc_text_by_name[name]
        from ..amortized import text_embed as te

        return te.screen_description_text(ctx)

    def _encode_text_context(
        self, ctx: dict[str, Any], observations: list[Observation]
    ) -> torch.Tensor:
        """Render description + observations to text, tokenize, encode -> (1, d_gene)."""
        from ..amortized.text_embed import render_context_text

        syms, hits = [], []
        for o in observations:
            lbl = o.label
            hit = bool(lbl.get("hit")) if isinstance(lbl, dict) else bool(lbl)
            syms.append(str(o.candidate))
            hits.append(1 if hit else 0)
        desc_text = self._desc_text(ctx) if self.net.cfg.use_description else ""
        text = render_context_text(desc_text, syms, hits)
        enc = self.tokenizer(
            [text], padding=True, truncation=True,
            max_length=self.net.cfg.text_max_tokens, return_tensors="pt",
        )
        return self.net.encode_text(
            enc["input_ids"].to(self.device), enc["attention_mask"].to(self.device)
        )

    # -- test-time context interventions -----------------------------------

    def _denoise(self, ids: list[int], hits: list[int],
                 rels: list[float]) -> tuple[list[int], list[int], int]:
        """Act on the weakest revealed hits. Returns (ids, hits, n_touched).

        ``rels`` holds the ``relevance_score`` the task revealed with each
        observation. Only hits are considered: a non-hit carries no strength to
        be weak about, and dropping one would also discard the "already tested"
        signal the context exists to convey.
        """
        if self.ctx_weak_hits == "keep":
            return ids, hits, 0
        hit_rels = [r for h, r in zip(hits, rels) if h and r == r]   # r == r drops NaN
        if len(hit_rels) < 4:          # too few to define a quantile meaningfully
            return ids, hits, 0
        thr = float(np.quantile(hit_rels, self.ctx_weak_q))
        out_ids, out_hits, touched = [], [], 0
        for i, h, r in zip(ids, hits, rels):
            weak = bool(h) and r == r and r <= thr
            if not weak:
                out_ids.append(i); out_hits.append(h); continue
            touched += 1
            if self.ctx_weak_hits == "drop":
                continue
            out_ids.append(i); out_hits.append(0)      # "negative"
        return out_ids, out_hits, touched

    def _context_variants(self, ids: list[int], hits: list[int]
                          ) -> list[tuple[list[int], list[int]]]:
        """K (ids, hits) context sets. Variant 0 is always the full context, so
        the ensemble strictly contains the unmodified model."""
        variants = [(ids, hits)]
        if self.ctx_ensemble == 1 or not ids:
            return variants
        keep = max(1, int(round(self.ctx_subsample * len(ids))))
        if keep >= len(ids):
            return variants * self.ctx_ensemble
        rng = np.random.default_rng((self.ctx_seed, len(ids)))
        for _ in range(self.ctx_ensemble - 1):
            sel = np.sort(rng.choice(len(ids), size=keep, replace=False))
            variants.append(([ids[i] for i in sel], [hits[i] for i in sel]))
        return variants

    def _fuse(self, scores: np.ndarray) -> np.ndarray:
        """(K, L) per-variant scores -> (L,) fused score, higher is better."""
        if scores.shape[0] == 1:
            return scores[0]
        if self.ctx_fuse == "mean":
            return scores.mean(axis=0)
        # RRF: rank 0 = best within each variant, then sum 1/(60 + rank).
        order = np.argsort(-scores, axis=1, kind="stable")
        ranks = np.empty_like(order)
        np.put_along_axis(ranks, order, np.arange(scores.shape[1])[None, :], axis=1)
        return (1.0 / (60.0 + ranks)).sum(axis=0)

    def name(self) -> str:
        return "amortized_ranker"

    def reset(self) -> None:
        self._meta = {}

    @torch.no_grad()
    def predict(
        self,
        observations: list[Observation],
        candidates: list[Any],
        task_context: dict[str, Any] | None = None,
    ) -> ModelPrediction:
        ctx = task_context or {}
        if not candidates:
            return ModelPrediction(scores={}, metadata={"name": self.name()})

        # Diagnostic ablation: drop the AL context entirely (static ranking).
        obs = [] if self.ignore_context else observations
        n_context = len(obs)
        n_weak = 0
        if self.encoder_type == "modernbert_text":
            desc_repr = self._encode_text_context(ctx, obs)  # (1, d_gene)
        else:
            if not self.net.cfg.use_description:
                # No-description ablation: the model ignores the desc embedding,
                # so skip the (potentially remote) embedding call and feed a zero.
                text_dim = self.net.desc_proj.in_features
                desc = torch.zeros(1, text_dim, device=self.device)
            else:
                desc = torch.from_numpy(self._desc_emb(ctx)).float().unsqueeze(0).to(self.device)

            # Build the observed-gene context from cumulative observations.
            ctx_ids: list[int] = []
            ctx_hits: list[int] = []
            ctx_rels: list[float] = []
            for o in obs:
                lbl = o.label
                is_d = isinstance(lbl, dict)
                hit = bool(lbl.get("hit")) if is_d else bool(lbl)
                ctx_ids.append(self.vocab.to_idx(str(o.candidate)))
                ctx_hits.append(1 if hit else 0)
                ctx_rels.append(float(lbl.get("relevance_score", float("nan")))
                                if is_d else float("nan"))

            ctx_ids, ctx_hits, n_weak = self._denoise(ctx_ids, ctx_hits, ctx_rels)
            variants = self._context_variants(ctx_ids, ctx_hits)
            K = len(variants)
            Lc = max(len(v[0]) for v in variants)

            if Lc:
                ci = torch.zeros(K, Lc, dtype=torch.int64, device=self.device)
                ch = torch.zeros(K, Lc, dtype=torch.int64, device=self.device)
                cpad = torch.ones(K, Lc, dtype=torch.bool, device=self.device)
                for b, (vi, vh) in enumerate(variants):
                    n = len(vi)
                    if not n:
                        continue
                    ci[b, :n] = torch.tensor(vi, dtype=torch.int64, device=self.device)
                    ch[b, :n] = torch.tensor(vh, dtype=torch.int64, device=self.device)
                    cpad[b, :n] = False
            else:
                ci = torch.zeros(K, 0, dtype=torch.int64, device=self.device)
                ch = torch.zeros(K, 0, dtype=torch.int64, device=self.device)
                cpad = torch.zeros(K, 0, dtype=torch.bool, device=self.device)

            desc_repr = self.net.encode(desc.expand(K, -1), ci, ch, cpad)  # (K, d_gene)

        K = desc_repr.shape[0]
        cand_ids = torch.tensor(
            [[self.vocab.to_idx(str(c)) for c in candidates]],
            dtype=torch.int64, device=self.device,
        ).expand(K, -1)
        if self.score_mode == "bias_only":
            per_variant = self.net.gene_bias[cand_ids].float().cpu().numpy()
        elif self.score_mode == "no_bias":
            per_variant = torch.einsum(
                "bd,bld->bl", desc_repr, self.net.gene_emb(cand_ids)
            ).float().cpu().numpy()
        else:
            per_variant = self.net.score_ids(desc_repr, cand_ids).float().cpu().numpy()
        scores = self._fuse(per_variant)

        score_map = {c: float(scores[i]) for i, c in enumerate(candidates)}
        self._meta = {"n_context": n_context, "n_candidates": len(candidates),
                      "ignore_context": self.ignore_context,
                      "ctx_ensemble": K, "ctx_fuse": self.ctx_fuse,
                      "ctx_weak_hits": self.ctx_weak_hits,
                      "ctx_weak_dropped": int(n_weak)}
        return ModelPrediction(
            scores=score_map,
            metadata={"name": self.name(), "checkpoint": self._ckpt, **self._meta},
        )


__all__ = ["AmortizedRankerModel"]
