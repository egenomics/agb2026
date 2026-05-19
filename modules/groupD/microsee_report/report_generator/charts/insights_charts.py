"""charts/insights_charts.py — per-chart explanation blocks (what / finding / pills).

generate_chart_explanations() returns a dict keyed by chart ID.  The template JS
auto-injects an ℹ button + collapsible panel for every chart that has an entry.

Kept separate from insights.py (section-level one-liners) because this module is
~700 lines and the two concerns are independent.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .preprocessing import (
    get_base_groups,
    get_patient_timepoints,
    get_unique_patients,
    sorted_timepoints,
)


def generate_chart_explanations(
    result: Any,
    chart_data: dict[str, Any],
    rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Return a what / finding / pills block for every chart in the report.

    Keys correspond to chart element IDs in the template (minus the "chart-" prefix
    where applicable).  The template JS auto-injects info buttons and panels.
    """
    taxa = result.taxa
    base_groups = get_base_groups(rows)
    patients = get_unique_patients(rows)
    timepoints = sorted_timepoints(rows)
    ex: dict[str, dict[str, Any]] = {}

    # ── Taxonomy ──────────────────────────────────────────────────────────────

    if taxa:
        means = {t: float(np.mean([float(r.get(t) or 0) for r in rows])) for t in taxa}
        sorted_t = sorted(means, key=lambda t: means[t], reverse=True)
        top3 = sorted_t[:3]
        top3_str = ", ".join(f"{t} ({means[t]:.1f}%)" for t in top3)
        top3_sum = sum(means[t] for t in top3)

        ex["composition"] = {
            "what": (
                "Stacked bar chart — each bar is one sample, each coloured segment is one bacterial family. "
                "Height shows relative abundance (%), so all bars reach 100%. "
                "Use the Timepoint buttons to isolate T0 (baseline) or T84 (post-intervention), "
                "the Group buttons to compare EAA vs Whey, and Top taxa to focus on the most abundant families. "
                "Families that shift between T0 and T84 may be responding to supplementation."
            ),
            "finding": (
                f"The three most abundant families are {top3_str}, accounting for "
                f"{top3_sum:.1f}% of the total microbiome on average across all {result.n_samples} samples. "
                "Use the T0/T84 filter to check whether their proportions shift after the intervention."
            ),
            "pills": [
                f"Top: {top3[0]} {means[top3[0]]:.1f}%",
                f"{result.n_taxa} families",
                f"{result.n_samples} samples",
            ],
        }

        ex["top_taxa"] = {
            "what": (
                "Horizontal bar chart showing mean relative abundance of each family across all samples, "
                "ranked from most to least dominant. This gives a quick snapshot of which microbes define "
                "this cohort's microbiome. Longer bars = more abundant on average."
            ),
            "finding": (
                f"{top3[0]} dominates at {means[top3[0]]:.1f}%, nearly "
                f"{'double' if means[top3[0]] > 2 * means[top3[1]] else 'more than'} "
                f"{top3[1]} ({means[top3[1]]:.1f}%). "
                f"The least abundant family ({sorted_t[-1]}) averages only {means[sorted_t[-1]]:.1f}%."
            ),
            "pills": [f"#{i + 1}: {t} {means[t]:.1f}%" for i, t in enumerate(top3)],
        }

        ex["donut"] = {
            "what": (
                "Donut chart showing the average relative abundance per group. "
                "Each arc segment represents one bacterial family. "
                "Hover over a segment to see the exact percentage. "
                "Compare the EAA and Whey donuts side-by-side to identify group-level differences."
            ),
            "finding": (
                f"{top3[0]} is the largest segment in all groups. "
                "Subtle arc-size differences between groups point to families worth investigating "
                "further in the Comparative section."
            ),
            "pills": ["Per-group averages", "Relative abundance (%)", "Hover for values"],
        }

        ex["sunburst"] = {
            "what": (
                "Sunburst chart — inner ring = treatment group, outer ring = bacterial family. "
                "Arc width is proportional to mean relative abundance. "
                "Click an inner segment to zoom into that group and see its family breakdown. "
                "Click the centre to zoom out."
            ),
            "finding": (
                f"The outer ring confirms {top3[0]} as the dominant family across groups. "
                "Look for arc-width asymmetry between groups in the outer ring — "
                "these indicate families that differ between EAA and Whey participants."
            ),
            "pills": ["Interactive — click to zoom", "Group → Family hierarchy", "Mean abundance"],
        }

    # ── Alpha diversity ───────────────────────────────────────────────────────

    bg_shan: dict[str, float] = {}
    for bg in base_groups:
        vals = [float(r.get("shannon") or 0) for r in rows if r.get("base_group", r["group"]) == bg]
        if vals:
            bg_shan[bg] = round(float(np.mean(vals)), 3)

    shan_desc = "; ".join(f"{bg}={v:.2f}" for bg, v in bg_shan.items())

    ex["alpha_strip"] = {
        "what": (
            "Strip chart (dot plot) showing the distribution of alpha diversity for each group. "
            "Every dot is one sample. The horizontal bar is the group mean. "
            "Alpha diversity measures how many different species live within a single person's sample "
            "and how evenly distributed they are. "
            "Shannon H′ is the most common metric — higher = more diverse and even community."
        ),
        "finding": (
            f"Group means: {shan_desc}. "
            "Switch the Metric button to compare Shannon, Simpson, Pielou evenness, Observed taxa, "
            "and Faith PD — consistent patterns across metrics are more reliable than any single measure."
        ),
        "pills": [f"{bg} mean {v:.2f}" for bg, v in list(bg_shan.items())[:3]]
        + ["Shannon H′ shown"],
    }

    ex["alpha_box"] = {
        "what": (
            "Box plot showing the distribution of alpha diversity per group. "
            "The box spans the interquartile range (25th–75th percentile); the line inside is the median; "
            "whiskers extend to 1.5× IQR; dots beyond whiskers are outliers. "
            "Significance brackets (if shown) give Wilcoxon paired p-values (T0 vs T84 within group) "
            "and Mann-Whitney p-values (between groups at T84)."
        ),
        "finding": (
            f"Group means: {shan_desc}. "
            "Check whether the confidence intervals overlap — non-overlapping boxes suggest a real difference. "
            "Brackets marked * (p < 0.05) or ** (p < 0.01) indicate statistically supported differences."
        ),
        "pills": ["IQR box", "Wilcoxon brackets", "Mann-Whitney between groups"],
    }

    ex["alpha_violin"] = {
        "what": (
            "Violin plot combining a box plot with kernel density estimation (the smooth shape). "
            "Wider sections of the violin = more samples at that diversity value. "
            "A narrow waist with a wide top indicates most samples have high diversity. "
            "Useful for detecting bimodal distributions (two peaks) that a box plot would miss."
        ),
        "finding": (
            f"Group means: {shan_desc}. "
            "Look for groups with wide, flat violins (high variability between patients) vs "
            "narrow, tall violins (patients with similar diversity levels)."
        ),
        "pills": ["Kernel density", "Median + IQR box", "Full distribution shape"],
    }

    ex["rarefaction"] = {
        "what": (
            "Rarefaction curve showing how many bacterial taxa are expected to be observed "
            "as sequencing depth (number of reads) increases. "
            "Curves that plateau indicate sufficient sequencing depth — all species have been captured. "
            "Curves that are still rising at the right edge suggest deeper sequencing would reveal more taxa. "
            "The shaded band shows ±1 standard deviation."
        ),
        "finding": (
            "If the curves plateau before the right edge, sequencing was sufficient for reliable comparisons. "
            "Groups with higher plateaus have more taxa on average — "
            "consistent with their alpha diversity scores."
        ),
        "pills": ["Log-scale x-axis", "Mean ± 1 SD", "Plateau = sufficient depth"],
    }

    ex["multimet"] = {
        "what": (
            "Combined chart showing two complementary alpha diversity metrics on the same plot. "
            "Bars (left axis) = Observed taxa count (raw richness). "
            "Diamonds (right axis) = Pielou J′ evenness (how equally distributed the reads are, 0–1). "
            "A sample can have many taxa (high richness) but still be dominated by one family (low evenness)."
        ),
        "finding": (
            "Samples with tall bars but low diamonds are dominated by a few taxa. "
            "Samples with short bars but high diamonds are species-poor but well-balanced. "
            "The ideal healthy microbiome tends to show both high richness and high evenness."
        ),
        "pills": ["Bars = Observed taxa", "Diamonds = Pielou J′", "Dual axis"],
    }

    # ── Beta diversity ────────────────────────────────────────────────────────

    bray: dict[str, Any] = chart_data.get("pcoa_bray") or {}
    jacc: dict[str, Any] = chart_data.get("pcoa_jaccard") or {}
    nmds: dict[str, Any] = chart_data.get("nmds") or {}
    pm: dict[str, Any] = chart_data.get("permanova") or {}

    if bray.get("pct1") is not None:
        top_driver: str = str(pm.get("top_name") or "individual identity")
        top_R2: float = float(pm.get("top_R2") or 0.0)

        ex["pcoa_bray"] = {
            "what": (
                "Principal Coordinates Analysis (PCoA) using Bray-Curtis dissimilarity. "
                "Each dot is one sample — dots closer together have more similar microbial communities. "
                "PC1 (horizontal) and PC2 (vertical) capture the most variation in the dataset. "
                "If samples cluster by colour (group), supplementation affected composition."
            ),
            "finding": (
                f"PC1 explains {bray['pct1']:.1f}% and PC2 {bray.get('pct2', 0):.1f}% of total variance. "
                f"The PERMANOVA confirms that {top_driver} is the primary driver of community structure "
                f"(R²={top_R2:.3f}). "
                "This is typical in microbiome studies — each person's gut is highly individual."
            ),
            "pills": [f"PC1 {bray['pct1']:.1f}%", f"PC2 {bray.get('pct2', 0):.1f}%", "Bray-Curtis"],
        }

    if jacc.get("pct1") is not None:
        ex["pcoa_jaccard"] = {
            "what": (
                "PCoA using Jaccard dissimilarity — presence/absence only, ignores abundance. "
                "Compare with Bray-Curtis: if patterns are similar, community structure is driven "
                "by which families are present; if different, abundance differences matter more."
            ),
            "finding": (
                f"PC1 explains {jacc['pct1']:.1f}% and PC2 {jacc.get('pct2', 0):.1f}% of variance. "
                "A 5% presence threshold was applied (families < 5% reads treated as absent). "
                "Compare clustering patterns with the Bray-Curtis plot above."
            ),
            "pills": [
                f"PC1 {jacc['pct1']:.1f}%",
                f"PC2 {jacc.get('pct2', 0):.1f}%",
                "Presence/absence",
                "5% threshold",
            ],
        }

    if nmds.get("pct1") is not None:
        ex["nmds"] = {
            "what": (
                "Non-metric Multidimensional Scaling (NMDS) on Bray-Curtis distances. "
                "Unlike PCoA, NMDS preserves rank-order relationships rather than exact distances — "
                "it is better at capturing non-linear community variation. "
                "Axes have no direct biological meaning; only relative positions matter."
            ),
            "finding": (
                f"NMDS1 explains {nmds['pct1']:.1f}% and NMDS2 {nmds.get('pct2', 0):.1f}% of variance. "
                "If NMDS and PCoA show similar groupings, the community structure is robust. "
                "Discrepancies suggest non-linear variation that PCoA misses."
            ),
            "pills": [
                f"NMDS1 {nmds['pct1']:.1f}%",
                f"NMDS2 {nmds.get('pct2', 0):.1f}%",
                "Rank-order",
            ],
        }

    ex["dendrogram"] = {
        "what": (
            "Hierarchical clustering dendrogram using average-linkage on Bray-Curtis distances. "
            "Each leaf (right end) is one sample. Samples that branch together early "
            "(branches close to the left = more similar) have more similar microbial communities. "
            "If T0 and T84 samples from the same patient cluster together, "
            "individual identity dominates over time."
        ),
        "finding": (
            "Look for same-patient pairs (e.g., EAA01_T0 and EAA01_T84) branching near each other — "
            "this confirms the personal microbiome fingerprint is stable over the study period. "
            "Cross-group clustering at T84 would suggest supplementation converged community composition."
        ),
        "pills": ["Average linkage", "Bray-Curtis", "Leaf = 1 sample"],
    }

    if taxa and patients:
        max_delta, max_taxon = 0.0, ""
        for p in patients:
            r0, r84 = get_patient_timepoints(rows, p)
            if r0 is None or r84 is None:
                continue
            tot0 = sum(float(r0.get(t) or 0) for t in taxa) or 1.0
            tot84 = sum(float(r84.get(t) or 0) for t in taxa) or 1.0
            for t in taxa:
                d = float(r84.get(t) or 0) / tot84 * 100 - float(r0.get(t) or 0) / tot0 * 100
                if abs(d) > abs(max_delta):
                    max_delta, max_taxon = d, t

        ex["delta_heatmap"] = {
            "what": (
                "Heatmap of Δ relative abundance (T84 − T0) for each patient × family cell. "
                "Red = family increased after supplementation; blue = decreased; white = no change. "
                "Rows (families) are sorted by largest absolute change across all patients. "
                "This reveals whether shifts are consistent (whole row red/blue) or patient-specific."
            ),
            "finding": (
                f"The family with the largest overall shift is {max_taxon} (Δ={max_delta:+.1f}%). "
                "Rows that are consistently red or blue across all patients suggest a treatment-wide effect. "
                "Patchy rows (mixed red/blue) indicate high inter-individual variability in response."
            ),
            "pills": [f"Max Δ: {max_taxon}", f"{max_delta:+.1f}%", "Red=increase · Blue=decrease"],
        }

    # ── Individual analysis ───────────────────────────────────────────────────

    if taxa and patients and len(timepoints) >= 2:
        increases, decreases = 0, 0
        for p in patients:
            r0, r84 = get_patient_timepoints(rows, p)
            if r0 and r84:
                if float(r84.get("shannon") or 0) > float(r0.get("shannon") or 0):
                    increases += 1
                else:
                    decreases += 1
        total_p = increases + decreases

        ex["paired_slope"] = {
            "what": (
                "Each line connects one patient's Shannon diversity at T0 (baseline) and T84 (end of study). "
                "Lines sloping upward = diversity increased; downward = decreased. "
                "The dashed line is the group mean trajectory. "
                "Shannon H′ ranges from 0 (single species) up — higher values = richer, more balanced community. "
                "Use the metric toggle to check whether Simpson diversity follows the same pattern."
            ),
            "finding": (
                f"{increases} of {total_p} patients showed increased Shannon diversity by T84; "
                f"{decreases} showed a decrease. "
                f"{'The majority gained diversity' if increases > decreases else 'The majority lost diversity'}, "
                "but check whether the group mean line (dashed) differs between EAA and Whey."
            ),
            "pills": [
                f"{increases} increased",
                f"{decreases} decreased",
                f"{total_p} patients",
                "Toggle Shannon/Simpson",
            ],
        }

    stab = chart_data.get("stability_bar", [])
    if stab and stab[0].get("x") and stab[0].get("y"):
        bc_vals = [float(v) for v in stab[0]["x"]]
        pts = stab[0]["y"]
        median_bc = float(np.median(bc_vals))
        n_stable = sum(1 for v in bc_vals if v < 0.2)

        ex["stability"] = {
            "what": (
                "Horizontal bar chart showing Bray-Curtis dissimilarity between each patient's "
                "T0 (baseline) and T84 (post-intervention) sample. "
                "BC dissimilarity ranges from 0 (identical communities) to 1 (completely different). "
                "Shorter bars = more stable microbiome over the study period."
            ),
            "finding": (
                f"Median stability score: {median_bc:.3f}. "
                f"{pts[0]} was most stable (BC={bc_vals[0]:.3f}); "
                f"{pts[-1]} showed the most change (BC={bc_vals[-1]:.3f}). "
                f"{n_stable}/{len(bc_vals)} patients maintained stable composition (BC < 0.2)."
            ),
            "pills": [
                f"Median BC {median_bc:.3f}",
                f"{n_stable}/{len(bc_vals)} stable",
                "0 = identical · 1 = different",
            ],
        }

    ex["diversity_rank"] = {
        "what": (
            "Samples are ranked left to right from lowest to highest alpha diversity. "
            "Circles (○) are T0 (baseline) samples; diamonds (◆) are T84 (post-intervention) samples. "
            "Colour indicates treatment group. "
            "If T84 diamonds consistently sit to the right of T0 circles for the same group, "
            "diversity increased after supplementation. "
            "Use the metric toggle to compare Shannon vs Simpson rankings."
        ),
        "finding": (
            "Look for whether T84 samples (diamonds) shift right relative to T0 (circles) within each group. "
            "A systematic rightward shift indicates supplementation increased diversity. "
            "Mixing of group colours along the rank axis suggests groups had similar baseline diversity."
        ),
        "pills": ["○ = T0 · ◆ = T84", "Ranked low→high", "Toggle metric"],
    }

    ex["radar"] = {
        "what": (
            "Radar (spider) chart showing the compositional profile of the selected patient. "
            "Each axis is one bacterial family. The filled shape = T0 (baseline); "
            "the dashed line = T84 (post-intervention); the dotted line = group mean T0 (reference). "
            "A larger filled area = more even composition across families. "
            "Select a different patient from the dropdown to compare individuals."
        ),
        "finding": (
            "The difference between the filled (T0) and dashed (T84) shapes shows how each patient's "
            "microbiome shifted during the study. "
            "The table on the right gives exact percentages and Δ values for each family — "
            "orange rows increased >2%, blue rows decreased >2%."
        ),
        "pills": ["Filled = T0", "Dashed = T84", "Dotted = group mean", "Per-patient selector"],
    }

    ex["nmds_traj"] = {
        "what": (
            "Ordination plot (PCoA on Bray-Curtis) showing each patient's trajectory from T0 to T84. "
            "Open circle (○) = baseline sample (T0); filled circle (●) = post-intervention (T84). "
            "The line connects the two timepoints for the same patient. "
            "Short lines = stable microbiome; long lines = large community shift."
        ),
        "finding": (
            "Short arrows indicate the microbiome was resistant to change over the study period. "
            "Arrows pointing toward a group cluster at T84 would suggest supplementation converged "
            "community composition. Look for whether EAA and Whey arrows point in similar directions."
        ),
        "pills": ["○ = T0 · ● = T84", "Line length = amount of change", "Per-patient trajectories"],
    }

    ex["faceted"] = {
        "what": (
            "Small multiples — one stacked bar per patient showing T0 and T84 side by side. "
            "Each coloured segment is one bacterial family; bar height = relative abundance (%). "
            "This allows direct visual comparison of each individual's microbiome before and after "
            "supplementation, and identifies outlier patients with unusual compositions."
        ),
        "finding": (
            "Compare T0 and T84 bars within each patient panel. "
            "Patients whose bar colours shift noticeably are the responders. "
            "Patients whose bars look nearly identical were microbiome-stable "
            "(consistent with low BC scores)."
        ),
        "pills": ["T0 vs T84 per patient", "Colour = family", "Click legend to isolate family"],
    }

    # ── Comparative ───────────────────────────────────────────────────────────

    # Extract top LFC hits from diff_abundance traces
    _da_traces = chart_data.get("diff_abundance") or []
    _lfc_up: tuple[str, float, str] = ("", 0.0, "")
    _lfc_dn: tuple[str, float, str] = ("", 0.0, "")
    for _tr in _da_traces:
        _bg = _tr.get("name", "")
        for _taxon, _v in zip(_tr.get("y", []), _tr.get("x", []), strict=False):
            _lfc_v = float(_v)
            if _lfc_v > _lfc_up[1]:
                _lfc_up = (str(_taxon), _lfc_v, _bg)
            if _lfc_v < _lfc_dn[1]:
                _lfc_dn = (str(_taxon), _lfc_v, _bg)

    ex["diff_abundance"] = {
        "what": (
            "Log2 Fold Change (LFC) bar chart comparing T84 vs T0 abundance for each family. "
            "Positive LFC = family increased after supplementation; negative = decreased. "
            "Log2 scale: LFC=1 means doubled, LFC=−1 means halved. "
            "Bars are grouped by treatment group so EAA and Whey responses can be compared side by side."
        ),
        "finding": (
            f"Largest increase: {_lfc_up[0]} (Log2FC={_lfc_up[1]:+.2f} in {_lfc_up[2]}); "
            f"largest decrease: {_lfc_dn[0]} (Log2FC={_lfc_dn[1]:+.2f} in {_lfc_dn[2]}). "
            "A large difference between EAA and Whey bars for the same family indicates "
            "a group-specific response — cross-check in the ANCOM and Volcano charts."
        )
        if _lfc_up[0]
        else (
            "Families on the right of zero increased during the study; those on the left decreased. "
            "A large difference in bar height between EAA and Whey for the same family indicates "
            "a group-specific response worth further investigation."
        ),
        "pills": (
            [f"Top +: {_lfc_up[0]}", f"Log2FC={_lfc_up[1]:+.2f}", "EAA vs Whey grouped"]
            if _lfc_up[0]
            else ["Log2 scale", "T84 vs T0", "EAA vs Whey grouped"]
        ),
    }

    # Count significant volcano hits (marker color "#D84E6A" = passed both thresholds)
    _vol_traces = chart_data.get("volcano") or []
    _sig_names_vol: list[str] = []
    for _tr in _vol_traces:
        _colors = (_tr.get("marker") or {}).get("color", [])
        _texts = _tr.get("text", [])
        for _color, _text in zip(_colors, _texts, strict=False):
            if _color == "#D84E6A":
                _sig_names_vol.append(f"{_text} ({_tr.get('name', '')})")
    _n_sig_vol = len(_sig_names_vol)
    _vol_top = ", ".join(_sig_names_vol[:3]) + (" …" if _n_sig_vol > 3 else "")

    ex["volcano"] = {
        "what": (
            "Volcano plot combining effect size (x = Log2 Fold Change) with statistical significance "
            "(y = −log10 p-value, higher = more significant). "
            "Each point is one bacterial family. "
            "Points in the upper right: significantly increased. Upper left: significantly decreased. "
            "Red points pass both thresholds — they changed a lot AND the change is statistically reliable."
        ),
        "finding": (
            f"{_n_sig_vol} {'family' if _n_sig_vol == 1 else 'families'} passed both thresholds "
            f"(|LFC| > 0.5 and FDR q < 0.1): {_vol_top}. "
            "These are the most biologically meaningful hits — large fold-change AND statistically robust. "
            "Cross-check against the ANCOM chart for CLR-corrected confirmation."
        )
        if _n_sig_vol > 0
        else (
            "No family passed both thresholds — changes were either small or high-variance. "
            "Broaden the LFC or FDR cut-offs, or inspect the ANCOM chart for CLR-corrected results."
        ),
        "pills": [
            f"{_n_sig_vol} significant {'hit' if _n_sig_vol == 1 else 'hits'}",
            "FDR q < 0.1 · |LFC| > 0.5",
            "Red = both thresholds",
        ],
    }

    # Count significant ANCOM hits (red = "#D84E6A" up, blue = "#4A7ED4" down)
    _ancom_traces = chart_data.get("ancom_style") or []
    _n_sig_ancom = sum(
        1
        for _tr in _ancom_traces
        for _c in (_tr.get("marker") or {}).get("color", [])
        if _c in ("#D84E6A", "#4A7ED4")
    )

    ex["ancom"] = {
        "what": (
            "CLR-transformed paired Wilcoxon test (an ANCOM-style analysis). "
            "Centred Log-Ratio (CLR) transformation corrects for the compositional nature of microbiome data — "
            "raw percentages are misleading because increasing one family mathematically decreases all others. "
            "Red bars = significantly increased (FDR q < 0.1); blue = significantly decreased."
        ),
        "finding": (
            f"{_n_sig_ancom} {'family reached' if _n_sig_ancom == 1 else 'families reached'} FDR q < 0.1 "
            f"after CLR transformation. "
            f"{'Overlap with the Volcano hits is the strongest evidence of a real shift. ' if _n_sig_ancom > 0 else ''}"
            "Discrepancies between CLR and LFC charts reveal composition-driven artefacts in the raw fold-changes."
        ),
        "pills": [
            f"{_n_sig_ancom} significant {'family' if _n_sig_ancom == 1 else 'families'}",
            "CLR-transformed",
            "Paired Wilcoxon · BH-FDR",
            "More robust than LFC",
        ],
    }

    ex["heatmap"] = {
        "what": (
            "Sample × family abundance heatmap. Rows = samples, columns = families. "
            "Colour intensity = relative abundance (%) — darker = more abundant. "
            "Samples and families are ordered to reveal clustering patterns. "
            "Blocks of similar colour indicate samples with consistent composition."
        ),
        "finding": (
            "Look for column-wise patterns (families abundant across a subgroup of patients) and "
            "row-wise patterns (patients with similar overall composition). "
            "A column that is dark only in T0 or only in T84 rows suggests a temporal shift."
        ),
        "pills": [
            "Sample × taxon matrix",
            "Colour = relative abundance (%)",
            "Hover for exact value",
        ],
    }

    ex["corr_matrix"] = {
        "what": (
            "Pairwise Pearson correlation matrix between bacterial family abundances. "
            "Each cell shows how strongly two families co-vary across all samples. "
            "Dark red = strong positive correlation (co-occurring families). "
            "Dark blue = strong negative correlation (competing families). "
            "White = no correlation."
        ),
        "finding": (
            "Positively correlated families tend to grow together — they may share ecological niches or "
            "benefit from the same dietary substrates. "
            "Negatively correlated families compete. "
            "NOTE: correlations in compositional data can be spurious; interpret with caution."
        ),
        "pills": [
            "Pearson r",
            "Red = co-occurring",
            "Blue = competing",
            "⚠ Compositional artefacts possible",
        ],
    }

    # ── Clinical ──────────────────────────────────────────────────────────────

    if result.has_clinical:
        corr_mwt: dict[str, Any] = chart_data.get("corr_mwt") or {}
        corr_il18: dict[str, Any] = chart_data.get("corr_il18") or {}

        # Count patients who improved on each clinical measure (T84 vs T0)
        _mwt_per_p: dict[str, dict[str, float]] = {}
        _il18_per_p: dict[str, dict[str, float]] = {}
        for _r in rows:
            _p = _r["patient"]
            _tp = _r.get("timepoint", "")
            _mwt_v = float(_r.get("sixmwt") or 0)
            _il18_v = float(_r.get("il18") or 0)
            if _mwt_v > 0:
                _mwt_per_p.setdefault(_p, {})[_tp] = _mwt_v
            if _il18_v > 0:
                _il18_per_p.setdefault(_p, {})[_tp] = _il18_v

        _n_mwt_tot = sum(1 for _pd in _mwt_per_p.values() if "T0" in _pd and "T84" in _pd)
        _n_mwt_imp = sum(
            1
            for _pd in _mwt_per_p.values()
            if "T0" in _pd and "T84" in _pd and _pd["T84"] > _pd["T0"]
        )
        _n_il18_tot = sum(1 for _pd in _il18_per_p.values() if "T0" in _pd and "T84" in _pd)
        _n_il18_imp = sum(
            1
            for _pd in _il18_per_p.values()
            if "T0" in _pd and "T84" in _pd and _pd["T84"] < _pd["T0"]
        )

        ex["clinical_mwt"] = {
            "what": (
                "Slopegraph showing each patient's 6-Minute Walk Test (6MWT) distance at T0 and T84. "
                "6MWT measures physical function — how far a person can walk in 6 minutes. "
                "Lines sloping upward = physical function improved; downward = declined. "
                "The dashed line shows the group mean. Higher 6MWT = better physical function."
            ),
            "finding": (
                f"{_n_mwt_imp} of {_n_mwt_tot} patients improved their 6MWT distance by T84. "
                "Compare EAA vs Whey mean lines (dashed) — a larger upward shift in one group "
                "suggests that supplement may better support physical function."
            )
            if _n_mwt_tot > 0
            else (
                "Upward-sloping lines indicate patients who improved their walking capacity during the study. "
                "Compare EAA vs Whey mean lines — a larger upward shift in one group suggests "
                "that supplement may better support physical function."
            ),
            "pills": [
                f"{_n_mwt_imp}/{_n_mwt_tot} improved" if _n_mwt_tot > 0 else "6MWT",
                "6MWT = 6-minute walk distance",
                "Higher = better function",
            ],
        }

        ex["clinical_il18"] = {
            "what": (
                "Slopegraph showing each patient's IL-18 cytokine level (pg/mL) at T0 and T84. "
                "IL-18 is a pro-inflammatory cytokine — lower levels indicate less inflammation. "
                "Lines sloping downward = inflammation reduced (beneficial); upward = increased. "
                "The dashed line shows the group mean."
            ),
            "finding": (
                f"{_n_il18_imp} of {_n_il18_tot} patients reduced their IL-18 level by T84. "
                "Compare EAA vs Whey mean lines (dashed) — a greater downward shift in one group "
                "suggests that supplement had a larger anti-inflammatory effect."
            )
            if _n_il18_tot > 0
            else (
                "Downward-sloping lines indicate patients who reduced systemic inflammation during the study. "
                "Compare EAA vs Whey to see which supplement had a larger anti-inflammatory effect. "
                "IL-18 reductions may reflect changes in gut permeability or microbiome composition."
            ),
            "pills": [
                f"{_n_il18_imp}/{_n_il18_tot} reduced" if _n_il18_tot > 0 else "IL-18",
                "IL-18 pro-inflammatory marker",
                "Lower = less inflammation",
            ],
        }

        if corr_mwt.get("r") is not None:
            r_mwt: float = float(corr_mwt["r"])
            p_mwt: float = float(corr_mwt.get("p") or 1.0)
            ex["corr_mwt"] = {
                "what": (
                    "Scatter plot of Shannon diversity (x) vs 6MWT distance (y) across all samples. "
                    "Each dot = one sample. The dashed line is the Pearson regression line. "
                    "A positive slope suggests patients with richer microbiomes tend to walk farther — "
                    "consistent with the gut-muscle axis hypothesis."
                ),
                "finding": (
                    f"Pearson r = {r_mwt:.2f}, p = {p_mwt:.3f}. "
                    f"{'Statistically significant (p < 0.05).' if p_mwt < 0.05 else 'Not statistically significant (p ≥ 0.05).'} "
                    f"{'A positive correlation suggests higher gut diversity associates with better physical function.' if r_mwt > 0 else 'A negative correlation was unexpected — inspect for outlier samples.'}"
                ),
                "pills": [
                    f"r = {r_mwt:.2f}",
                    f"p = {p_mwt:.3f}",
                    "✱ Significant" if p_mwt < 0.05 else "Not significant",
                    "Pearson",
                ],
            }

        if corr_il18.get("r") is not None:
            r_il18: float = float(corr_il18["r"])
            p_il18: float = float(corr_il18.get("p") or 1.0)
            ex["corr_il18"] = {
                "what": (
                    "Scatter plot of Shannon diversity (x) vs IL-18 cytokine level (y) across all samples. "
                    "A negative slope (r < 0) would mean patients with higher gut diversity have lower "
                    "inflammation — consistent with the microbiome's role in immune regulation."
                ),
                "finding": (
                    f"Pearson r = {r_il18:.2f}, p = {p_il18:.3f}. "
                    f"{'Statistically significant (p < 0.05).' if p_il18 < 0.05 else 'Not statistically significant (p ≥ 0.05).'} "
                    f"{'A negative r suggests higher diversity associates with lower inflammation.' if r_il18 < 0 else 'A positive r was unexpected — check for confounders or outliers.'}"
                ),
                "pills": [
                    f"r = {r_il18:.2f}",
                    f"p = {p_il18:.3f}",
                    "✱ Significant" if p_il18 < 0.05 else "Not significant",
                    "Pearson",
                ],
            }

        ex["taxa_clinical"] = {
            "what": (
                "Spearman correlation heatmap between the change in each bacterial family (Δ relative abundance, "
                "T84 − T0) and the change in each clinical outcome (Δ 6MWT and Δ IL-18). "
                "Red = family increase associates with clinical improvement; blue = opposite. "
                "Stars indicate statistical significance (* p < 0.05, ** p < 0.01) after BH-FDR correction."
            ),
            "finding": (
                "Families with starred red cells in the 6MWT column are the ones whose increase "
                "associates with improved physical function. "
                "Families with starred blue cells in the IL-18 column associate with reduced inflammation. "
                "These are the most actionable findings for future mechanistic studies."
            ),
            "pills": [
                "Spearman ρ",
                "BH-FDR corrected",
                "* p<0.05 · ** p<0.01",
                "Δ taxon × Δ clinical",
            ],
        }

    # ── Longitudinal ──────────────────────────────────────────────────────────

    # Extract per-group Shannon change (T84 - T0) from longitudinal traces
    _long_traces = chart_data.get("longitudinal") or []
    _long_changes: list[tuple[str, float]] = []
    for _tr in _long_traces:
        _xs = _tr.get("x", [])
        _ys = _tr.get("y", [])
        if "T0" in _xs and "T84" in _xs:
            _delta = float(_ys[_xs.index("T84")]) - float(_ys[_xs.index("T0")])
            _long_changes.append((str(_tr.get("name", "")), _delta))
    _long_str = "; ".join(f"{_g} {_d:+.3f}" for _g, _d in _long_changes)

    ex["longitudinal"] = {
        "what": (
            "Line chart showing mean Shannon diversity per group at each timepoint. "
            "Lines connect the group mean across time. "
            "An upward slope indicates the group's average diversity increased during the study; "
            "flat lines indicate stability. "
            "This is a summary view — the LME chart adds confidence intervals and individual patient lines."
        ),
        "finding": (
            f"Shannon H′ change T0→T84: {_long_str}. "
            "A positive value means average diversity increased; a diverging pattern between groups "
            "suggests a supplement-specific effect on the gut microbiome."
        )
        if _long_changes
        else (
            "Compare the slope and endpoint height of EAA vs Whey lines. "
            "A diverging pattern (one group increasing while the other stays flat) would suggest "
            "the supplements have different effects on gut microbiome diversity."
        ),
        "pills": ([f"{_g}: {_d:+.3f}" for _g, _d in _long_changes[:2]] + ["Shannon H′"])
        if _long_changes
        else ["Group mean per timepoint", "Shannon H′", "Summary view"],
    }

    ex["lme"] = {
        "what": (
            "LME-style trajectory plot: mean diversity ± 95% confidence interval per group (thick lines + band), "
            "with individual patient lines underneath (thin lines). "
            "If bands for two groups don't overlap, the difference is statistically robust. "
            "The Wilcoxon p-value tests whether T0 and T84 distributions differ within each group."
        ),
        "finding": (
            "Wide confidence bands indicate high inter-individual variability — patients respond differently. "
            "Overlapping bands mean the groups cannot be reliably distinguished at this sample size. "
            "Use the metric toggle to compare Shannon, Simpson, and other alpha metrics."
        ),
        "pills": ["Mean ± 95% CI", "Individual patient lines", "Wilcoxon p annotated"],
    }

    # ── Statistics ────────────────────────────────────────────────────────────

    if pm.get("rows"):
        top_name: str = str(pm.get("top_name") or "Unknown")
        _perm_r2: float = float(pm.get("top_R2") or 0.0)
        pm_rows: list[Any] = pm.get("rows") or []
        supp_row: list[Any] | None = next(
            (r for r in pm_rows if "Group" in str(r[0]) or "Suppl" in str(r[0])),
            None,
        )
        supp_p: float = float(supp_row[3]) if supp_row is not None else 1.0
        supp_R2: float = float(supp_row[1]) if supp_row is not None else 0.0

        ex["permanova"] = {
            "what": (
                "PERMANOVA (Permutational MANOVA) tests which experimental factor explains the most variation "
                "in overall microbial community composition. "
                "R² = fraction of total variance explained. "
                "p-value = probability the observed R² arose by chance (99 permutations). "
                "Three factors are tested: supplementation group, timepoint, and individual patient identity."
            ),
            "finding": (
                f"{top_name} explains the most variance (R²={_perm_r2:.3f}). "
                f"Supplementation group accounts for R²={supp_R2:.3f} "
                f"(p={supp_p:.3f} — {'significant ✓' if supp_p < 0.05 else 'not significant'}). "
                "Individual identity dominating over group is expected and does not mean the intervention failed."
            ),
            "pills": [
                f"{top_name} R²={_perm_r2:.3f}",
                f"Group R²={supp_R2:.3f}",
                f"Group p={supp_p:.3f}",
                "99 permutations",
            ],
        }

    ex["stats_summary"] = {
        "what": (
            "Summary table of mean ± standard deviation for each alpha diversity metric, "
            "broken down by group. "
            "This is the standard 'Table 1' style summary used in microbiome publications. "
            "Use it to compare group-level diversity at a glance and to check whether groups "
            "had similar baseline diversity (T0 rows should be comparable)."
        ),
        "finding": (
            "Groups with similar T0 means are well-matched at baseline. "
            "A larger difference between T0 and T84 rows within the same group suggests "
            "the intervention had a measurable effect on that diversity metric. "
            "All five metrics are shown; consistent patterns across metrics are more reliable."
        ),
        "pills": ["Mean ± SD per group", "All 5 alpha metrics", "Groups × timepoints"],
    }

    return ex
