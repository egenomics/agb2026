"""Unit tests for charts/stats_helpers.py — pure statistical functions.

All functions here use only stdlib + numpy, so these tests run with no
additional dependencies beyond the base package install.
"""

import numpy as np
import pytest

from report_generator.charts.stats_helpers import (
    bh_fdr,
    mannwhitney_p,
    pearson_p,
    sig_label,
    spearman_r,
    welch_ttest_p,
    wilcoxon_p,
)

# ── wilcoxon_p ────────────────────────────────────────────────────────────────

class TestWilcoxonP:
    def test_identical_returns_one(self):
        a = [1.0, 2.0, 3.0]
        assert wilcoxon_p(a, a) == pytest.approx(1.0)

    def test_clearly_different_is_small(self):
        # Need n≥6 to reach p<0.05 with exact Wilcoxon enumeration
        a = [1.0, 1.1, 0.9, 1.05, 0.95, 1.0]
        b = [10.0, 10.1, 9.9, 10.05, 9.95, 10.0]
        p = wilcoxon_p(a, b)
        assert p < 0.05, f"Expected p<0.05 for clearly different groups, got {p}"

    def test_symmetric(self):
        a = [1.0, 2.0, 3.0, 4.0]
        b = [5.0, 6.0, 7.0, 8.0]
        assert wilcoxon_p(a, b) == pytest.approx(wilcoxon_p(b, a), abs=0.01)

    def test_too_few_pairs_returns_one(self):
        assert wilcoxon_p([1.0], [2.0]) == 1.0

    def test_bounded_zero_to_one(self):
        rng = np.random.default_rng(42)
        a = rng.normal(0, 1, 15).tolist()
        b = rng.normal(1, 1, 15).tolist()
        p = wilcoxon_p(a, b)
        assert 0.0 <= p <= 1.0

    def test_large_n_uses_normal_approximation(self):
        a = list(range(21))
        b = [x + 5 for x in a]
        p = wilcoxon_p(a, b)
        assert 0.0 <= p <= 1.0


# ── mannwhitney_p ─────────────────────────────────────────────────────────────

class TestMannWhitneyP:
    def test_identical_distributions_near_one(self):
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        b = [1.0, 2.0, 3.0, 4.0, 5.0]
        p = mannwhitney_p(a, b)
        assert p > 0.5, f"Identical distributions should give high p, got {p}"

    def test_clearly_separated_is_small(self):
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        b = [100.0, 200.0, 300.0, 400.0, 500.0]
        p = mannwhitney_p(a, b)
        assert p < 0.05, f"Clearly separated groups should give p<0.05, got {p}"

    def test_too_small_returns_one(self):
        assert mannwhitney_p([1.0], [2.0, 3.0]) == 1.0
        assert mannwhitney_p([1.0, 2.0], [3.0]) == 1.0

    def test_bounded(self):
        rng = np.random.default_rng(7)
        a = rng.normal(0, 1, 10).tolist()
        b = rng.normal(2, 1, 10).tolist()
        p = mannwhitney_p(a, b)
        assert 0.0 <= p <= 1.0

    def test_symmetric(self):
        a = [1.0, 3.0, 5.0]
        b = [2.0, 4.0, 6.0]
        assert mannwhitney_p(a, b) == pytest.approx(mannwhitney_p(b, a), abs=1e-6)


# ── sig_label ─────────────────────────────────────────────────────────────────

class TestSigLabel:
    def test_very_significant(self):
        assert "***" in sig_label(0.0001)

    def test_significant(self):
        assert "**" in sig_label(0.005)

    def test_marginal(self):
        label = sig_label(0.03)
        assert "*" in label and "***" not in label and "**" not in label

    def test_not_significant(self):
        assert "ns" in sig_label(0.5)

    def test_boundary_0_05(self):
        # p=0.05 is NOT significant (strict <0.05)
        assert "ns" in sig_label(0.05)

    def test_boundary_0_001(self):
        # p=0.001 is NOT in the *** tier (strict <0.001)
        label = sig_label(0.001)
        assert "**" in label and "***" not in label


# ── welch_ttest_p ─────────────────────────────────────────────────────────────

class TestWelchTtestP:
    def test_same_mean_high_p(self):
        a = [5.0, 5.0, 5.0, 5.0, 5.0]
        b = [5.0, 5.0, 5.0, 5.0, 5.0]
        assert welch_ttest_p(a, b) == pytest.approx(1.0)

    def test_clearly_different_low_p(self):
        # Groups must have non-zero variance (se=0 when std=0 -> function returns 1.0 by design)
        rng = np.random.default_rng(1)
        a = (rng.normal(1.0,   0.1, 20)).tolist()
        b = (rng.normal(100.0, 0.1, 20)).tolist()
        p = welch_ttest_p(a, b)
        assert p < 0.001, f"Expected p<0.001, got {p}"

    def test_too_small_returns_one(self):
        assert welch_ttest_p([1.0], [2.0, 3.0]) == 1.0

    def test_bounded(self):
        rng = np.random.default_rng(99)
        a = rng.normal(0, 1, 15).tolist()
        b = rng.normal(0.5, 1, 15).tolist()
        p = welch_ttest_p(a, b)
        assert 0.0 <= p <= 1.0


# ── pearson_p ─────────────────────────────────────────────────────────────────

class TestPearsonP:
    def test_perfect_correlation_low_p(self):
        p = pearson_p(0.999, n=30)
        assert p < 0.001

    def test_zero_correlation_high_p(self):
        p = pearson_p(0.0, n=30)
        assert p > 0.9

    def test_too_few_samples_returns_one(self):
        assert pearson_p(0.9, n=3) == 1.0

    def test_bounded(self):
        for r in [-0.9, -0.5, 0.0, 0.5, 0.9]:
            p = pearson_p(r, n=20)
            assert 0.0 <= p <= 1.0, f"r={r} gave p={p} outside [0,1]"


# ── spearman_r ────────────────────────────────────────────────────────────────

class TestSpearmanR:
    def test_perfectly_correlated(self):
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        rho, p = spearman_r(a, a)
        assert rho == pytest.approx(1.0)
        assert p < 0.05

    def test_perfectly_anticorrelated(self):
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        b = [5.0, 4.0, 3.0, 2.0, 1.0]
        rho, p = spearman_r(a, b)
        assert rho == pytest.approx(-1.0)
        assert p < 0.05

    def test_too_few_returns_zero_one(self):
        rho, p = spearman_r([1.0, 2.0], [3.0, 4.0])
        assert rho == 0.0
        assert p == 1.0

    def test_rho_bounded(self):
        rng = np.random.default_rng(17)
        a = rng.normal(0, 1, 20).tolist()
        b = rng.normal(0, 1, 20).tolist()
        rho, p = spearman_r(a, b)
        assert -1.0 <= rho <= 1.0
        assert 0.0 <= p <= 1.0

    def test_known_value(self):
        # Monotone ranking: rho should be close to 1
        a = [1.0, 3.0, 5.0, 7.0, 9.0]
        b = [2.0, 4.0, 6.0, 8.0, 10.0]
        rho, _ = spearman_r(a, b)
        assert rho == pytest.approx(1.0)


# ── bh_fdr ────────────────────────────────────────────────────────────────────

class TestBhFdr:
    def test_empty_input(self):
        assert bh_fdr([]) == []

    def test_single_value_unchanged(self):
        assert bh_fdr([0.03]) == [0.03]

    def test_all_significant_stay_bounded(self):
        qs = bh_fdr([0.001, 0.002, 0.003])
        assert all(0.0 <= q <= 1.0 for q in qs)

    def test_monotone_non_decreasing(self):
        ps = [0.01, 0.04, 0.03, 0.20, 0.06]
        qs = bh_fdr(ps)
        # q-values in the original order — the smallest p should have the smallest q
        order = sorted(range(len(ps)), key=lambda i: ps[i])
        q_sorted = [qs[i] for i in order]
        assert q_sorted == sorted(q_sorted), "BH q-values should be non-decreasing in p-order"

    def test_length_preserved(self):
        ps = [0.05, 0.01, 0.20, 0.001, 0.10]
        assert len(bh_fdr(ps)) == len(ps)

    def test_all_one_stays_one(self):
        qs = bh_fdr([1.0, 1.0, 1.0])
        assert all(q == 1.0 for q in qs)

    def test_known_example(self):
        # 5 p-values — BH at alpha=0.05 rejects the two smallest
        # p = [0.001, 0.01, 0.05, 0.10, 0.20]  sorted
        # q[0] = 0.001*5/1=0.005, q[1]=0.01*5/2=0.025, q[2]=0.05*5/3≈0.083, ...
        ps = [0.10, 0.001, 0.05, 0.20, 0.01]
        qs = bh_fdr(ps)
        # Index 1 (p=0.001) should have the smallest q
        assert qs[1] == min(qs)
        assert all(0.0 <= q <= 1.0 for q in qs)

    def test_large_n_still_bounded(self):
        rng = np.random.default_rng(0)
        ps = rng.uniform(0, 1, 100).tolist()
        qs = bh_fdr(ps)
        assert len(qs) == 100
        assert all(0.0 <= q <= 1.0 for q in qs)


# ── cross-function consistency ────────────────────────────────────────────────

class TestCrossFunction:
    """Sanity-check that stat functions agree directionally on the same data."""

    def test_mw_and_welch_agree_direction(self):
        a = [1.0, 1.5, 2.0, 1.2, 1.8]
        b = [5.0, 6.0, 5.5, 6.5, 5.2]
        p_mw     = mannwhitney_p(a, b)
        p_welch  = welch_ttest_p(a, b)
        # Both should detect a large difference (p < 0.05)
        assert p_mw < 0.05,    f"Mann-Whitney p={p_mw} should be <0.05"
        assert p_welch < 0.05, f"Welch p={p_welch} should be <0.05"

    def test_spearman_and_pearson_agree_on_linear(self):
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        y = [x_i * 2 + 1 for x_i in x]
        rho, p_sp = spearman_r(x, y)
        p_pe      = pearson_p(rho, n=len(x))
        assert rho > 0.95
        assert p_sp < 0.01
        assert p_pe < 0.01
