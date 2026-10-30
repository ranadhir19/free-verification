"""Shared metric definitions for the rebuild (PREREGISTRATION_REBUILD.md s.5).

Used by every measurement script AND by validate.py, so the printed numbers
and their recomputation share one definition source. Inputs are always the
deposited artifacts (edge lists, descriptors, deposited model scores).
"""
import math

import numpy as np
from scipy.stats import rankdata


def filter_metrics(labels, keep):
    """labels: 0/1 array (1 = true/potent). keep: boolean array.
    rho = 1 - retained_wrong/total_wrong; sigma = retained_true/total_true;
    bits = log2(sigma/(1-rho)) from unrounded values."""
    labels = np.asarray(labels)
    keep = np.asarray(keep, dtype=bool)
    tot_t = int((labels == 1).sum())
    tot_w = int((labels == 0).sum())
    ret_t = int((labels[keep] == 1).sum())
    ret_w = int((labels[keep] == 0).sum())
    rho = 1.0 - ret_w / tot_w if tot_w else float("nan")
    sigma = ret_t / tot_t if tot_t else float("nan")
    bits = (math.log2(sigma / (1.0 - rho))
            if 0 < rho < 1 and sigma > 0 else float("nan"))
    return {"rho": rho, "sigma": sigma,
            "p_prime": ret_w / (ret_t + ret_w) if (ret_t + ret_w) else float("nan"),
            "bits": bits, "ret_true": ret_t, "ret_wrong": ret_w,
            "ret_total": ret_t + ret_w, "tot_true": tot_t, "tot_wrong": tot_w}


def auc_rank(scores, labels):
    """Mann-Whitney AUC with average ranks (ties handled)."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels)
    n1 = int((labels == 1).sum())
    n0 = int((labels == 0).sum())
    r = rankdata(scores)
    return float((r[labels == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def sweep(values, labels, min_sigma=0.50):
    """Frozen threshold sweep (pre-reg s.5.3): candidate cuts = unique
    observed values; directions '>=' scanned first, then '<=', cuts
    ascending; feasible iff sigma >= min_sigma and retained wrong >= 1;
    maximize bits, FIRST maximum in scan order wins.
    Returns dict with cut, direction and filter_metrics."""
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    cuts = np.unique(values)
    best = None
    for direction in (">=", "<="):
        for c in cuts:
            keep = values >= c if direction == ">=" else values <= c
            m = filter_metrics(labels, keep)
            if m["sigma"] >= min_sigma and m["ret_wrong"] >= 1 \
                    and math.isfinite(m["bits"]):
                if best is None or m["bits"] > best["bits"]:
                    best = dict(m, cut=float(c), direction=direction)
    return best


def wilson(x, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = x / n
    den = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (c - h, c + h)


def hanley_mcneil_se(auc, n_pos, n_neg):
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc * auc)
           + (n_neg - 1) * (q2 - auc * auc)) / (n_pos * n_neg)
    return math.sqrt(var)


def cluster_bootstrap_bits(labels, keep, clusters, B=2000, seed=20261030):
    """Percentile cluster bootstrap for bits of a fixed filter (pre-reg
    s.5.7). Resamples clusters with replacement; rho/sigma recomputed on
    the resampled multiset; degenerate resamples skipped and counted."""
    labels = np.asarray(labels)
    keep = np.asarray(keep, dtype=bool)
    clusters = np.asarray(clusters)
    uniq, inv = np.unique(clusters, return_inverse=True)
    nc = len(uniq)
    ct = np.bincount(inv, weights=(labels == 1), minlength=nc)
    cw = np.bincount(inv, weights=(labels == 0), minlength=nc)
    kt = np.bincount(inv, weights=((labels == 1) & keep), minlength=nc)
    kw = np.bincount(inv, weights=((labels == 0) & keep), minlength=nc)
    rng = np.random.default_rng(seed)
    vals = []
    skipped = 0
    for _ in range(B):
        draw = rng.integers(0, nc, size=nc)
        Ct = ct[draw].sum(); Cw = cw[draw].sum()
        Kt = kt[draw].sum(); Kw = kw[draw].sum()
        if Ct == 0 or Cw == 0 or Kt == 0 or Kw == 0:
            skipped += 1
            continue
        rho = 1.0 - Kw / Cw
        sigma = Kt / Ct
        if not (0 < rho < 1):
            skipped += 1
            continue
        vals.append(math.log2(sigma / (1.0 - rho)))
    vals = np.array(vals)
    if len(vals) == 0:
        # Every resample degenerate (e.g. the filter retains zero true
        # edges in the population, so sigma = 0 identically): the interval
        # is undefined; report nulls with valid = 0 (pre-reg s.5.7's
        # skip-and-count rule taken to its limit).
        return {"lo": None, "hi": None, "valid": 0,
                "skipped": int(skipped), "n_clusters": int(nc)}
    return {"lo": float(np.percentile(vals, 2.5)),
            "hi": float(np.percentile(vals, 97.5)),
            "valid": int(len(vals)), "skipped": int(skipped),
            "n_clusters": int(nc)}
