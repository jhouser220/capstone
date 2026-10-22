
import numpy as np
from scipy.special import expit

def map_state_precision(z0, Sig_inv, W_o, b_o, y_o):
    """Damped-Newton MAP with general prior precision (mirrors core.map_state)."""
    z0 = np.asarray(z0, np.float64); x = z0.copy()
    W_o = np.asarray(W_o, np.float64); b_o = np.asarray(b_o, np.float64); y_o = np.asarray(y_o, np.float64)
    if not len(y_o):
        return x, np.linalg.cholesky(Sig_inv)
    def terms(q):
        a = W_o @ q + b_o; p = expit(a); d = q - z0
        return (float(np.logaddexp(0, a).sum() - y_o @ a + 0.5 * (d @ (Sig_inv @ d))),
                W_o.T @ (p - y_o) + Sig_inv @ d,
                (W_o.T * (p * (1.0 - p))) @ W_o + Sig_inv)
    for _ in range(30):
        v, g, H = terms(x)
        if np.max(np.abs(g)) < 1e-7: break
        step = np.linalg.solve(H, g); rate = 1.0
        while rate > 2**-25 and terms(x - rate * step)[0] > v - 1e-4 * rate * (g @ step): rate *= 0.5
        x -= rate * step
    _, _, H = terms(x)
    return x, np.linalg.cholesky(H)


def rank_fusion(source, corrected):
    from scipy.stats import rankdata
    return -.5 * (rankdata(-np.asarray(source), method="average") +
                  rankdata(-np.asarray(corrected), method="average"))


def posterior_scores(prior, precision, observed_design, observed_bias, labels,
                     candidate_design, candidate_bias, q, *, screen_id, draw, step):
    """E159 fit and 128 antithetic draws; q is applied exactly once."""
    import hashlib
    seeds = (2026101701, 2026101702)
    if draw not in (0, 1) or not 1 <= step <= 10:
        raise ValueError("The frozen comparison uses two repetitions and ten decisions")
    mu, chol = map_state_precision(prior, precision, observed_design, observed_bias, labels)
    gradient = observed_design.T @ (expit(observed_design @ mu + observed_bias) - labels) + precision @ (mu-prior)
    if np.max(np.abs(gradient)) >= 1e-5:
        raise ValueError("Posterior fit failed stationarity check")
    seed = int(hashlib.sha256(f"{seeds[draw]}|{screen_id}|fill|{step}|particles".encode()).hexdigest()[:16],16)
    eps = np.random.default_rng(seed).normal(size=(64,len(mu)))
    eps = np.concatenate((eps,-eps))
    particles = mu[None] + np.linalg.solve(chol.T,eps.T).T
    prob = np.clip(expit(particles @ candidate_design.T + candidate_bias),1e-9,1-1e-9).mean(0)
    q = np.asarray(q, dtype=float)
    if q.shape != prob.shape or not np.isfinite(q).all() or np.any((q<0)|(q>1)):
        raise ValueError("Invalid candidate availability probabilities")
    return q*prob, {"mode":mu.tolist(),"gradient_inf":float(np.max(np.abs(gradient))),
                   "measured_labels":len(labels),"particles":128,"step":step}
