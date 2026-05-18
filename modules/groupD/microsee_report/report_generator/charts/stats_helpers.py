"""charts/stats_helpers.py — pure statistical test functions (no scipy dependency).

All implementations use only Python stdlib + numpy so the package runs on
minimal HPC environments where scipy is unavailable.

Previously these functions were scattered across stats.py, comparative.py, and
clinical.py — centralising them eliminates duplication and makes each test
independently testable.
"""
from __future__ import annotations

import math
import numpy as np


# ── Shared normal-distribution helper ────────────────────────────────────────

def _phi_complement(z: float) -> float:
    """Upper-tail probability of the standard normal: P(Z > |z|).

    Uses the Abramowitz & Stegun rational approximation 26.2.17 (max err 7.5e-8).
    """
    az = abs(z)
    t  = 1.0 / (1.0 + 0.2316419 * az)
    pd = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    return pd * math.exp(-0.5 * az * az) / math.sqrt(2.0 * math.pi)


# ── Wilcoxon signed-rank test ─────────────────────────────────────────────────

def wilcoxon_p(a: list[float], b: list[float]) -> float:
    """Two-tailed Wilcoxon signed-rank p-value.

    Exact enumeration for n ≤ 20; normal approximation with continuity
    correction for larger samples.
    """
    diffs = [float(x) - float(y) for x, y in zip(a, b) if float(x) != float(y)]
    n = len(diffs)
    if n < 2:
        return 1.0

    abs_d      = [abs(d) for d in diffs]
    sorted_idx = sorted(range(n), key=lambda i: abs_d[i])
    ranks      = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n - 1 and abs_d[sorted_idx[j]] == abs_d[sorted_idx[j + 1]]:
            j += 1
        avg = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[sorted_idx[k]] = avg
        i = j + 1

    w_plus  = sum(r for r, d in zip(ranks, diffs) if d > 0)
    total_w = n * (n + 1) / 2.0
    w_obs   = min(w_plus, total_w - w_plus)

    if n <= 20:
        cnt = sum(
            1 for mask in range(1 << n)
            if (s := sum(ranks[i] for i in range(n) if (mask >> i) & 1)) <= w_obs
            or s >= total_w - w_obs
        )
        return round(cnt / (1 << n), 4)

    mu    = total_w / 2.0
    sigma = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    z     = (w_obs - mu + 0.5) / sigma
    return round(min(1.0, 2.0 * _phi_complement(z)), 4)


# ── Mann-Whitney U test ───────────────────────────────────────────────────────

def mannwhitney_p(a: list[float], b: list[float]) -> float:
    """Two-tailed Mann-Whitney U p-value (normal approximation)."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 1.0

    combined = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    n        = len(combined)
    ranks    = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n - 1 and combined[j][0] == combined[j + 1][0]:
            j += 1
        avg = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1

    r_a   = sum(r for r, (_, g) in zip(ranks, combined) if g == 0)
    u_a   = r_a - na * (na + 1) / 2.0
    u     = min(u_a, na * nb - u_a)
    mu    = na * nb / 2.0
    sigma = math.sqrt(na * nb * (na + nb + 1) / 12.0)
    if sigma == 0:
        return 1.0

    z = (u - mu) / sigma
    return round(min(1.0, 2.0 * _phi_complement(z)), 4)


# ── Significance label ────────────────────────────────────────────────────────

def sig_label(p: float) -> str:
    """Human-readable significance annotation for a p-value."""
    if p < 0.001: return "p<0.001 ***"
    if p < 0.01:  return f"p={p:.3f} **"
    if p < 0.05:  return f"p={p:.3f} *"
    return f"p={p:.3f} ns"


# ── Welch's t-test (incomplete beta via continued fraction) ───────────────────

def welch_ttest_p(a: list[float], b: list[float]) -> float:
    """Two-tailed Welch's t-test p-value (unequal variances, numpy only)."""
    ar, br = np.array(a, dtype=float), np.array(b, dtype=float)
    na, nb = len(ar), len(br)
    if na < 2 or nb < 2:
        return 1.0

    va, vb = float(np.var(ar, ddof=1)), float(np.var(br, ddof=1))
    se = math.sqrt(va / na + vb / nb)
    if se == 0:
        return 1.0

    t  = abs(float(np.mean(ar)) - float(np.mean(br))) / se
    df = (va / na + vb / nb) ** 2 / (
        (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    )
    x = df / (df + t * t)
    if x <= 0: return 0.0
    if x >= 1: return 1.0

    def _ibeta(xv: float, av: float, bv: float) -> float:
        if xv <= 0: return 0.0
        if xv >= 1: return 1.0
        lbt = (av * math.log(xv) + bv * math.log(1 - xv)
               - math.lgamma(av) - math.lgamma(bv) + math.lgamma(av + bv))

        def _cf(xx: float, aa: float, bb: float) -> float:
            fpmin, eps = 1e-30, 3e-7
            c2, d2 = 1.0, 1.0 - (aa + bb) * xx / (aa + 1.0)
            d2 = fpmin if abs(d2) < fpmin else d2
            d2 = 1.0 / d2
            h  = d2
            for m in range(1, 201):
                m2  = 2 * m
                aa2 = m * (bb - m) * xx / ((aa + m2 - 1) * (aa + m2))
                d2  = 1.0 + aa2 * d2; d2 = fpmin if abs(d2) < fpmin else d2
                c2  = 1.0 + aa2 / c2; c2 = fpmin if abs(c2) < fpmin else c2
                d2  = 1.0 / d2; h *= d2 * c2
                aa2 = -(aa + m) * (aa + bb + m) * xx / ((aa + m2) * (aa + m2 + 1))
                d2  = 1.0 + aa2 * d2; d2 = fpmin if abs(d2) < fpmin else d2
                c2  = 1.0 + aa2 / c2; c2 = fpmin if abs(c2) < fpmin else c2
                d2  = 1.0 / d2; delta = d2 * c2; h *= delta
                if abs(delta - 1.0) < eps:
                    break
            return h

        if xv < (av + 1) / (av + bv + 2):
            return math.exp(lbt) * _cf(xv, av, bv) / av
        return 1.0 - math.exp(lbt) * _cf(1.0 - xv, bv, av) / bv

    return max(0.0, min(1.0, _ibeta(x, df / 2.0, 0.5)))


# ── Pearson correlation p-value ───────────────────────────────────────────────

def pearson_p(r_val: float, n: int) -> float:
    """Two-tailed p-value for Pearson r (Fisher z-transformation)."""
    if n < 4:
        return 1.0
    z = math.atanh(min(max(r_val, -0.9999), 0.9999)) * math.sqrt(n - 3)
    return min(1.0, 2.0 * _phi_complement(z))


# ── Spearman rank correlation ─────────────────────────────────────────────────

def spearman_r(a: list[float], b: list[float]) -> tuple[float, float]:
    """Spearman ρ and two-tailed p-value."""
    n = len(a)
    if n < 3:
        return 0.0, 1.0

    def _ranks(v: list[float]) -> list[float]:
        sv = sorted(enumerate(v), key=lambda x: x[1])
        r  = [0.0] * n
        i  = 0
        while i < n:
            j = i
            while j < n - 1 and sv[j][1] == sv[j + 1][1]:
                j += 1
            avg = (i + j + 2) / 2.0
            for k in range(i, j + 1):
                r[sv[k][0]] = avg
            i = j + 1
        return r

    ra, rb = _ranks(a), _ranks(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    num    = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da     = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db     = math.sqrt(sum((y - mb) ** 2 for y in rb))
    rho    = round(num / (da * db), 3) if da and db else 0.0
    return rho, pearson_p(rho, n)


# ── Benjamini-Hochberg FDR ────────────────────────────────────────────────────

def bh_fdr(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction. Returns q-values same length as input."""
    n = len(p_values)
    if n == 0:
        return []

    order = sorted(range(n), key=lambda i: p_values[i])
    q = [0.0] * n
    for rank, i in enumerate(order, 1):
        q[i] = min(1.0, p_values[i] * n / rank)

    # Enforce monotonicity: q non-decreasing in p-value order
    min_q = 1.0
    for i in sorted(range(n), key=lambda i: p_values[i], reverse=True):
        min_q = min(min_q, q[i])
        q[i]  = min_q

    return [round(v, 3) for v in q]
