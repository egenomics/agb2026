"""charts/insights_charts.py — per-chart explanation blocks (what / finding / pills).

generate_chart_explanations() returns a dict keyed by chart ID.  The template JS
auto-injects an ℹ button + collapsible panel for every chart that has an entry.

Each section is a private _explain_* function so this file stays navigable.
Static "what" descriptions live as module-level constants; functions only compute
the dynamic "finding" and "pills" from the data.
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

# ── Tiny helper ───────────────────────────────────────────────────────────────


def _entry(what: str, finding: str, pills: list[str]) -> dict[str, Any]:
    return {"what": what, "finding": finding, "pills": pills}


# ── Static "what" descriptions (one per chart) ────────────────────────────────

_WHAT_COMPOSITION = (
    "Stacked bar chart — each bar is one sample, each coloured segment is one bacterial family. "
    "Height shows relative abundance (%), so all bars reach 100%. "
    "Use the Timepoint buttons to isolate T0 (baseline) or T84 (post-intervention), "
    "the Group buttons to compare EAA vs Whey, and Top taxa to focus on the most abundant families. "
    "Families that shift between T0 and T84 may be responding to supplementation."
)
_WHAT_TOP_TAXA = (
    "Horizontal bar chart showing mean relative abundance of each family across all samples, "
    "ranked from most to least dominant. This gives a quick snapshot of which microbes define "
    "this cohort's microbiome. Longer bars = more abundant on average."
)
_WHAT_DONUT = (
    "Donut chart showing the average relative abundance per group. "
    "Each arc segment represents one bacterial family. "
    "Hover over a segment to see the exact percentage. "
    "Compare the EAA and Whey donuts side-by-side to identify group-level differences."
)
_WHAT_SUNBURST = (
    "Sunburst chart — inner ring = treatment group, outer ring = bacterial family. "
    "Arc width is proportional to mean relative abundance. "
    "Click an inner segment to zoom into that group and see its family breakdown. "
    "Click the centre to zoom out."
)
_WHAT_ALPHA_STRIP = (
    "Strip chart (dot plot) showing the distribution of alpha diversity for each group. "
    "Every dot is one sample. The horizontal bar is the group mean. "
    "Alpha diversity measures how many different species live within a single person's sample "
    "and how evenly distributed they are. "
    "Shannon H′ is the most common metric — higher = more diverse and even community."
)
_WHAT_ALPHA_BOX = (
    "Box plot showing the distribution of alpha diversity per group. "
    "The box spans the interquartile range (25th–75th percentile); the line inside is the median; "
    "whiskers extend to 1.5× IQR; dots beyond whiskers are outliers. "
    "Significance brackets (if shown) give Wilcoxon paired p-values (T0 vs T84 within group) "
    "and Mann-Whitney p-values (between groups at T84)."
)
_WHAT_ALPHA_VIOLIN = (
    "Violin plot combining a box plot with kernel density estimation (the smooth shape). "
    "Wider sections of the violin = more samples at that diversity value. "
    "A narrow waist with a wide top indicates most samples have high diversity. "
    "Useful for detecting bimodal distributions (two peaks) that a box plot would miss."
)
_WHAT_RAREFACTION = (
    "Rarefaction curve showing how many bacterial taxa are expected to be observed "
    "as sequencing depth (number of reads) increases. "
    "Curves that plateau indicate sufficient sequencing depth — all species have been captured. "
    "Curves that are still rising at the right edge suggest deeper sequencing would reveal more taxa. "
    "The shaded band shows ±1 standard deviation."
)
_WHAT_MULTIMET = (
    "Combined chart showing two complementary alpha diversity metrics on the same plot. "
    "Bars (left axis) = Observed taxa count (raw richness). "
    "Diamonds (right axis) = Pielou J′ evenness (how equally distributed the reads are, 0–1). "
    "A sample can have many taxa (high richness) but still be dominated by one family (low evenness)."
)
_WHAT_PCOA_BRAY = (
    "Principal Coordinates Analysis (PCoA) using Bray-Curtis dissimilarity. "
    "Each dot is one sample — dots closer together have more similar microbial communities. "
    "PC1 (horizontal) and PC2 (vertical) capture the most variation in the dataset. "
    "If samples cluster by colour (group), supplementation affected composition."
)
_WHAT_PCOA_JACCARD = (
    "PCoA using Jaccard dissimilarity — presence/absence only, ignores abundance. "
    "Compare with Bray-Curtis: if patterns are similar, community structure is driven "
    "by which families are present; if different, abundance differences matter more."
)
_WHAT_NMDS = (
    "Non-metric Multidimensional Scaling (NMDS) on Bray-Curtis distances. "
    "Unlike PCoA, NMDS preserves rank-order relationships rather than exact distances — "
    "it is better at capturing non-linear community variation. "
    "Axes have no direct biological meaning; only relative positions matter."
)
_WHAT_DENDROGRAM = (
    "Hierarchical clustering dendrogram using average-linkage on Bray-Curtis distances. "
    "Each leaf (right end) is one sample. Samples that branch together early "
    "(branches close to the left = more similar) have more similar microbial communities. "
    "If T0 and T84 samples from the same patient cluster together, "
    "individual identity dominates over time."
)
_WHAT_DELTA_HEATMAP = (
    "Heatmap of Δ relative abundance (T84 − T0) for each patient × family cell. "
    "Red = family increased after supplementation; blue = decreased; white = no change. "
    "Rows (families) are sorted by largest absolute change across all patients. "
    "This reveals whether shifts are consistent (whole row red/blue) or patient-specific."
)
_WHAT_PAIRED_SLOPE = (
    "Each line connects one patient's Shannon diversity at T0 (baseline) and T84 (end of study). "
    "Lines sloping upward = diversity increased; downward = decreased. "
    "The dashed line is the group mean trajectory. "
    "Shannon H′ ranges from 0 (single species) up — higher values = richer, more balanced community. "
    "Use the metric toggle to check whether Simpson diversity follows the same pattern."
)
_WHAT_STABILITY = (
    "Horizontal bar chart showing Bray-Curtis dissimilarity between each patient's "
    "T0 (baseline) and T84 (post-intervention) sample. "
    "BC dissimilarity ranges from 0 (identical communities) to 1 (completely different). "
    "Shorter bars = more stable microbiome over the study period."
)
_WHAT_DIVERSITY_RANK = (
    "Samples are ranked left to right from lowest to highest alpha diversity. "
    "Circles (○) are T0 (baseline) samples; diamonds (◆) are T84 (post-intervention) samples. "
    "Colour indicates treatment group. "
    "If T84 diamonds consistently sit to the right of T0 circles for the same group, "
    "diversity increased after supplementation. "
    "Use the metric toggle to compare Shannon vs Simpson rankings."
)
_WHAT_RADAR = (
    "Radar (spider) chart showing the compositional profile of the selected patient. "
    "Each axis is one bacterial family. The filled shape = T0 (baseline); "
    "the dashed line = T84 (post-intervention); the dotted line = group mean T0 (reference). "
    "A larger filled area = more even composition across families. "
    "Select a different patient from the dropdown to compare individuals."
)
_WHAT_NMDS_TRAJ = (
    "Ordination plot (PCoA on Bray-Curtis) showing each patient's trajectory from T0 to T84. "
    "Open circle (○) = baseline sample (T0); filled circle (●) = post-intervention (T84). "
    "The line connects the two timepoints for the same patient. "
    "Short lines = stable microbiome; long lines = large community shift."
)
_WHAT_FACETED = (
    "Small multiples — one stacked bar per patient showing T0 and T84 side by side. "
    "Each coloured segment is one bacterial family; bar height = relative abundance (%). "
    "This allows direct visual comparison of each individual's microbiome before and after "
    "supplementation, and identifies outlier patients with unusual compositions."
)
_WHAT_DIFF_ABUNDANCE = (
    "Log2 Fold Change (LFC) bar chart comparing T84 vs T0 abundance for each family. "
    "Positive LFC = family increased after supplementation; negative = decreased. "
    "Log2 scale: LFC=1 means doubled, LFC=−1 means halved. "
    "Bars are grouped by treatment group so EAA and Whey responses can be compared side by side."
)
_WHAT_VOLCANO = (
    "Volcano plot combining effect size (x = Log2 Fold Change) with statistical significance "
    "(y = −log10 p-value, higher = more significant). "
    "Each point is one bacterial family. "
    "Points in the upper right: significantly increased. Upper left: significantly decreased. "
    "Red points pass both thresholds — they changed a lot AND the change is statistically reliable."
)
_WHAT_ANCOM = (
    "CLR-transformed paired Wilcoxon test (an ANCOM-style analysis). "
    "Centred Log-Ratio (CLR) transformation corrects for the compositional nature of microbiome data — "
    "raw percentages are misleading because increasing one family mathematically decreases all others. "
    "Red bars = significantly increased (FDR q < 0.1); blue = significantly decreased."
)
_WHAT_HEATMAP = (
    "Sample × family abundance heatmap. Rows = samples, columns = families. "
    "Colour intensity = relative abundance (%) — darker = more abundant. "
    "Samples and families are ordered to reveal clustering patterns. "
    "Blocks of similar colour indicate samples with consistent composition."
)
_WHAT_CORR_MATRIX = (
    "Pairwise Pearson correlation matrix between bacterial family abundances. "
    "Each cell shows how strongly two families co-vary across all samples. "
    "Dark red = strong positive correlation (co-occurring families). "
    "Dark blue = strong negative correlation (competing families). "
    "White = no correlation."
)
_WHAT_CLINICAL_MWT = (
    "Slopegraph showing each patient's 6-Minute Walk Test (6MWT) distance at T0 and T84. "
    "6MWT measures physical function — how far a person can walk in 6 minutes. "
    "Lines sloping upward = physical function improved; downward = declined. "
    "The dashed line shows the group mean. Higher 6MWT = better physical function."
)
_WHAT_CLINICAL_IL18 = (
    "Slopegraph showing each patient's IL-18 cytokine level (pg/mL) at T0 and T84. "
    "IL-18 is a pro-inflammatory cytokine — lower levels indicate less inflammation. "
    "Lines sloping downward = inflammation reduced (beneficial); upward = increased. "
    "The dashed line shows the group mean."
)
_WHAT_CORR_MWT = (
    "Scatter plot of Shannon diversity (x) vs 6MWT distance (y) across all samples. "
    "Each dot = one sample. The dashed line is the Pearson regression line. "
    "A positive slope suggests patients with richer microbiomes tend to walk farther — "
    "consistent with the gut-muscle axis hypothesis."
)
_WHAT_CORR_IL18 = (
    "Scatter plot of Shannon diversity (x) vs IL-18 cytokine level (y) across all samples. "
    "A negative slope (r < 0) would mean patients with higher gut diversity have lower "
    "inflammation — consistent with the microbiome's role in immune regulation."
)
_WHAT_TAXA_CLINICAL = (
    "Spearman correlation heatmap between the change in each bacterial family (Δ relative abundance, "
    "T84 − T0) and the change in each clinical outcome (Δ 6MWT and Δ IL-18). "
    "Red = family increase associates with clinical improvement; blue = opposite. "
    "Stars indicate statistical significance (* p < 0.05, ** p < 0.01) after BH-FDR correction."
)
_WHAT_LONGITUDINAL = (
    "Line chart showing mean Shannon diversity per group at each timepoint. "
    "Lines connect the group mean across time. "
    "An upward slope indicates the group's average diversity increased during the study; "
    "flat lines indicate stability. "
    "This is a summary view — the LME chart adds confidence intervals and individual patient lines."
)
_WHAT_LME = (
    "LME-style trajectory plot: mean diversity ± 95% confidence interval per group (thick lines + band), "
    "with individual patient lines underneath (thin lines). "
    "If bands for two groups don't overlap, the difference is statistically robust. "
    "The Wilcoxon p-value tests whether T0 and T84 distributions differ within each group."
)
_WHAT_PERMANOVA = (
    "PERMANOVA (Permutational MANOVA) tests which experimental factor explains the most variation "
    "in overall microbial community composition. "
    "R² = fraction of total variance explained. "
    "p-value = probability the observed R² arose by chance (99 permutations). "
    "Three factors are tested: supplementation group, timepoint, and individual patient identity."
)
_WHAT_STATS_SUMMARY = (
    "Summary table of mean ± standard deviation for each alpha diversity metric, "
    "broken down by group. "
    "This is the standard 'Table 1' style summary used in microbiome publications. "
    "Use it to compare group-level diversity at a glance and to check whether groups "
    "had similar baseline diversity (T0 rows should be comparable)."
)

# ── Section helpers ───────────────────────────────────────────────────────────


def _explain_taxonomy(
    result: Any, taxa: list[str], rows: list[dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    if not taxa:
        return {}
    means = {t: float(np.mean([float(r.get(t) or 0) for r in rows])) for t in taxa}
    sorted_t = sorted(means, key=lambda t: means[t], reverse=True)
    top3 = sorted_t[:3]
    top3_str = ", ".join(f"{t} ({means[t]:.1f}%)" for t in top3)
    top3_sum = sum(means[t] for t in top3)

    return {
        "composition": _entry(
            _WHAT_COMPOSITION,
            f"The three most abundant families are {top3_str}, accounting for "
            f"{top3_sum:.1f}% of the total microbiome on average across all {result.n_samples} samples. "
            "Use the T0/T84 filter to check whether their proportions shift after the intervention.",
            [
                f"Top: {top3[0]} {means[top3[0]]:.1f}%",
                f"{result.n_taxa} families",
                f"{result.n_samples} samples",
            ],
        ),
        "top_taxa": _entry(
            _WHAT_TOP_TAXA,
            f"{top3[0]} dominates at {means[top3[0]]:.1f}%, nearly "
            f"{'double' if means[top3[0]] > 2 * means[top3[1]] else 'more than'} "
            f"{top3[1]} ({means[top3[1]]:.1f}%). "
            f"The least abundant family ({sorted_t[-1]}) averages only {means[sorted_t[-1]]:.1f}%.",
            [f"#{i + 1}: {t} {means[t]:.1f}%" for i, t in enumerate(top3)],
        ),
        "donut": _entry(
            _WHAT_DONUT,
            f"{top3[0]} is the largest segment in all groups. "
            "Subtle arc-size differences between groups point to families worth investigating "
            "further in the Comparative section.",
            ["Per-group averages", "Relative abundance (%)", "Hover for values"],
        ),
        "sunburst": _entry(
            _WHAT_SUNBURST,
            f"The outer ring confirms {top3[0]} as the dominant family across groups. "
            "Look for arc-width asymmetry between groups in the outer ring — "
            "these indicate families that differ between EAA and Whey participants.",
            ["Interactive — click to zoom", "Group → Family hierarchy", "Mean abundance"],
        ),
    }


def _explain_alpha(base_groups: list[str], rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    bg_shan: dict[str, float] = {}
    for bg in base_groups:
        vals = [float(r.get("shannon") or 0) for r in rows if r.get("base_group", r["group"]) == bg]
        if vals:
            bg_shan[bg] = round(float(np.mean(vals)), 3)
    shan_desc = "; ".join(f"{bg}={v:.2f}" for bg, v in bg_shan.items())
    shan_pills = [f"{bg} mean {v:.2f}" for bg, v in list(bg_shan.items())[:3]] + [
        "Shannon H′ shown"
    ]

    return {
        "alpha_strip": _entry(
            _WHAT_ALPHA_STRIP,
            f"Group means: {shan_desc}. "
            "Switch the Metric button to compare Shannon, Simpson, Pielou evenness, Observed taxa, "
            "and Faith PD — consistent patterns across metrics are more reliable than any single measure.",
            shan_pills,
        ),
        "alpha_box": _entry(
            _WHAT_ALPHA_BOX,
            f"Group means: {shan_desc}. "
            "Check whether the confidence intervals overlap — non-overlapping boxes suggest a real difference. "
            "Brackets marked * (p < 0.05) or ** (p < 0.01) indicate statistically supported differences.",
            ["IQR box", "Wilcoxon brackets", "Mann-Whitney between groups"],
        ),
        "alpha_violin": _entry(
            _WHAT_ALPHA_VIOLIN,
            f"Group means: {shan_desc}. "
            "Look for groups with wide, flat violins (high variability between patients) vs "
            "narrow, tall violins (patients with similar diversity levels).",
            ["Kernel density", "Median + IQR box", "Full distribution shape"],
        ),
        "rarefaction": _entry(
            _WHAT_RAREFACTION,
            "If the curves plateau before the right edge, sequencing was sufficient for reliable comparisons. "
            "Groups with higher plateaus have more taxa on average — "
            "consistent with their alpha diversity scores.",
            ["Log-scale x-axis", "Mean ± 1 SD", "Plateau = sufficient depth"],
        ),
        "multimet": _entry(
            _WHAT_MULTIMET,
            "Samples with tall bars but low diamonds are dominated by a few taxa. "
            "Samples with short bars but high diamonds are species-poor but well-balanced. "
            "The ideal healthy microbiome tends to show both high richness and high evenness.",
            ["Bars = Observed taxa", "Diamonds = Pielou J′", "Dual axis"],
        ),
    }


def _explain_beta(
    chart_data: dict[str, Any],
    taxa: list[str],
    patients: list[str],
    rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    ex: dict[str, dict[str, Any]] = {}
    bray: dict[str, Any] = chart_data.get("pcoa_bray") or {}
    jacc: dict[str, Any] = chart_data.get("pcoa_jaccard") or {}
    nmds: dict[str, Any] = chart_data.get("nmds") or {}
    pm: dict[str, Any] = chart_data.get("permanova") or {}

    if bray.get("pct1") is not None:
        top_driver = str(pm.get("top_name") or "individual identity")
        top_R2 = float(pm.get("top_R2") or 0.0)
        ex["pcoa_bray"] = _entry(
            _WHAT_PCOA_BRAY,
            f"PC1 explains {bray['pct1']:.1f}% and PC2 {bray.get('pct2', 0):.1f}% of total variance. "
            f"The PERMANOVA confirms that {top_driver} is the primary driver of community structure "
            f"(R²={top_R2:.3f}). "
            "This is typical in microbiome studies — each person's gut is highly individual.",
            [f"PC1 {bray['pct1']:.1f}%", f"PC2 {bray.get('pct2', 0):.1f}%", "Bray-Curtis"],
        )

    if jacc.get("pct1") is not None:
        ex["pcoa_jaccard"] = _entry(
            _WHAT_PCOA_JACCARD,
            f"PC1 explains {jacc['pct1']:.1f}% and PC2 {jacc.get('pct2', 0):.1f}% of variance. "
            "A 5% presence threshold was applied (families < 5% reads treated as absent). "
            "Compare clustering patterns with the Bray-Curtis plot above.",
            [
                f"PC1 {jacc['pct1']:.1f}%",
                f"PC2 {jacc.get('pct2', 0):.1f}%",
                "Presence/absence",
                "5% threshold",
            ],
        )

    if nmds.get("pct1") is not None:
        ex["nmds"] = _entry(
            _WHAT_NMDS,
            f"NMDS1 explains {nmds['pct1']:.1f}% and NMDS2 {nmds.get('pct2', 0):.1f}% of variance. "
            "If NMDS and PCoA show similar groupings, the community structure is robust. "
            "Discrepancies suggest non-linear variation that PCoA misses.",
            [f"NMDS1 {nmds['pct1']:.1f}%", f"NMDS2 {nmds.get('pct2', 0):.1f}%", "Rank-order"],
        )

    ex["dendrogram"] = _entry(
        _WHAT_DENDROGRAM,
        "Look for same-patient pairs (e.g., EAA01_T0 and EAA01_T84) branching near each other — "
        "this confirms the personal microbiome fingerprint is stable over the study period. "
        "Cross-group clustering at T84 would suggest supplementation converged community composition.",
        ["Average linkage", "Bray-Curtis", "Leaf = 1 sample"],
    )

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
        ex["delta_heatmap"] = _entry(
            _WHAT_DELTA_HEATMAP,
            f"The family with the largest overall shift is {max_taxon} (Δ={max_delta:+.1f}%). "
            "Rows that are consistently red or blue across all patients suggest a treatment-wide effect. "
            "Patchy rows (mixed red/blue) indicate high inter-individual variability in response.",
            [f"Max Δ: {max_taxon}", f"{max_delta:+.1f}%", "Red=increase · Blue=decrease"],
        )

    return ex


def _explain_individual(
    chart_data: dict[str, Any],
    taxa: list[str],
    patients: list[str],
    timepoints: list[str],
    rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    ex: dict[str, dict[str, Any]] = {}

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
        majority = (
            "The majority gained diversity"
            if increases > decreases
            else "The majority lost diversity"
        )
        ex["paired_slope"] = _entry(
            _WHAT_PAIRED_SLOPE,
            f"{increases} of {total_p} patients showed increased Shannon diversity by T84; "
            f"{decreases} showed a decrease. {majority}, "
            "but check whether the group mean line (dashed) differs between EAA and Whey.",
            [
                f"{increases} increased",
                f"{decreases} decreased",
                f"{total_p} patients",
                "Toggle Shannon/Simpson",
            ],
        )

    stab = chart_data.get("stability_bar", [])
    if stab and stab[0].get("x") and stab[0].get("y"):
        bc_vals = [float(v) for v in stab[0]["x"]]
        pts = stab[0]["y"]
        median_bc = float(np.median(bc_vals))
        n_stable = sum(1 for v in bc_vals if v < 0.2)
        ex["stability"] = _entry(
            _WHAT_STABILITY,
            f"Median stability score: {median_bc:.3f}. "
            f"{pts[0]} was most stable (BC={bc_vals[0]:.3f}); "
            f"{pts[-1]} showed the most change (BC={bc_vals[-1]:.3f}). "
            f"{n_stable}/{len(bc_vals)} patients maintained stable composition (BC < 0.2).",
            [
                f"Median BC {median_bc:.3f}",
                f"{n_stable}/{len(bc_vals)} stable",
                "0 = identical · 1 = different",
            ],
        )

    ex["diversity_rank"] = _entry(
        _WHAT_DIVERSITY_RANK,
        "Look for whether T84 samples (diamonds) shift right relative to T0 (circles) within each group. "
        "A systematic rightward shift indicates supplementation increased diversity. "
        "Mixing of group colours along the rank axis suggests groups had similar baseline diversity.",
        ["○ = T0 · ◆ = T84", "Ranked low→high", "Toggle metric"],
    )
    ex["radar"] = _entry(
        _WHAT_RADAR,
        "The difference between the filled (T0) and dashed (T84) shapes shows how each patient's "
        "microbiome shifted during the study. "
        "The table on the right gives exact percentages and Δ values for each family — "
        "orange rows increased >2%, blue rows decreased >2%.",
        ["Filled = T0", "Dashed = T84", "Dotted = group mean", "Per-patient selector"],
    )
    ex["nmds_traj"] = _entry(
        _WHAT_NMDS_TRAJ,
        "Short arrows indicate the microbiome was resistant to change over the study period. "
        "Arrows pointing toward a group cluster at T84 would suggest supplementation converged "
        "community composition. Look for whether EAA and Whey arrows point in similar directions.",
        ["○ = T0 · ● = T84", "Line length = amount of change", "Per-patient trajectories"],
    )
    ex["faceted"] = _entry(
        _WHAT_FACETED,
        "Compare T0 and T84 bars within each patient panel. "
        "Patients whose bar colours shift noticeably are the responders. "
        "Patients whose bars look nearly identical were microbiome-stable "
        "(consistent with low BC scores).",
        ["T0 vs T84 per patient", "Colour = family", "Click legend to isolate family"],
    )
    return ex


def _explain_comparative(chart_data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    # Largest LFC hit per direction
    da_traces = chart_data.get("diff_abundance") or []
    lfc_up: tuple[str, float, str] = ("", 0.0, "")
    lfc_dn: tuple[str, float, str] = ("", 0.0, "")
    for tr in da_traces:
        bg = tr.get("name", "")
        for taxon, v in zip(tr.get("y", []), tr.get("x", []), strict=False):
            lfc_v = float(v)
            if lfc_v > lfc_up[1]:
                lfc_up = (str(taxon), lfc_v, bg)
            if lfc_v < lfc_dn[1]:
                lfc_dn = (str(taxon), lfc_v, bg)

    # Significant volcano hits (red marker = passed both thresholds)
    vol_traces = chart_data.get("volcano") or []
    sig_vol: list[str] = []
    for tr in vol_traces:
        colors = (tr.get("marker") or {}).get("color", [])
        for color, text in zip(colors, tr.get("text", []), strict=False):
            if color == "#D84E6A":
                sig_vol.append(f"{text} ({tr.get('name', '')})")
    n_sig_vol = len(sig_vol)
    vol_top = ", ".join(sig_vol[:3]) + (" …" if n_sig_vol > 3 else "")

    # Significant ANCOM hits
    ancom_traces = chart_data.get("ancom_style") or []
    n_sig_ancom = sum(
        1
        for tr in ancom_traces
        for c in (tr.get("marker") or {}).get("color", [])
        if c in ("#D84E6A", "#4A7ED4")
    )

    diff_finding = (
        f"Largest increase: {lfc_up[0]} (Log2FC={lfc_up[1]:+.2f} in {lfc_up[2]}); "
        f"largest decrease: {lfc_dn[0]} (Log2FC={lfc_dn[1]:+.2f} in {lfc_dn[2]}). "
        "A large difference between EAA and Whey bars for the same family indicates "
        "a group-specific response — cross-check in the ANCOM and Volcano charts."
        if lfc_up[0]
        else "Families on the right of zero increased during the study; those on the left decreased. "
        "A large difference in bar height between EAA and Whey for the same family indicates "
        "a group-specific response worth further investigation."
    )
    vol_finding = (
        f"{n_sig_vol} {'family' if n_sig_vol == 1 else 'families'} passed both thresholds "
        f"(|LFC| > 0.5 and FDR q < 0.1): {vol_top}. "
        "These are the most biologically meaningful hits — large fold-change AND statistically robust. "
        "Cross-check against the ANCOM chart for CLR-corrected confirmation."
        if n_sig_vol > 0
        else "No family passed both thresholds — changes were either small or high-variance. "
        "Broaden the LFC or FDR cut-offs, or inspect the ANCOM chart for CLR-corrected results."
    )

    return {
        "diff_abundance": _entry(
            _WHAT_DIFF_ABUNDANCE,
            diff_finding,
            [f"Top +: {lfc_up[0]}", f"Log2FC={lfc_up[1]:+.2f}", "EAA vs Whey grouped"]
            if lfc_up[0]
            else ["Log2 scale", "T84 vs T0", "EAA vs Whey grouped"],
        ),
        "volcano": _entry(
            _WHAT_VOLCANO,
            vol_finding,
            [
                f"{n_sig_vol} significant {'hit' if n_sig_vol == 1 else 'hits'}",
                "FDR q < 0.1 · |LFC| > 0.5",
                "Red = both thresholds",
            ],
        ),
        "ancom": _entry(
            _WHAT_ANCOM,
            f"{n_sig_ancom} {'family reached' if n_sig_ancom == 1 else 'families reached'} FDR q < 0.1 "
            "after CLR transformation. "
            f"{'Overlap with the Volcano hits is the strongest evidence of a real shift. ' if n_sig_ancom > 0 else ''}"
            "Discrepancies between CLR and LFC charts reveal composition-driven artefacts in the raw fold-changes.",
            [
                f"{n_sig_ancom} significant {'family' if n_sig_ancom == 1 else 'families'}",
                "CLR-transformed",
                "Paired Wilcoxon · BH-FDR",
                "More robust than LFC",
            ],
        ),
        "heatmap": _entry(
            _WHAT_HEATMAP,
            "Look for column-wise patterns (families abundant across a subgroup of patients) and "
            "row-wise patterns (patients with similar overall composition). "
            "A column that is dark only in T0 or only in T84 rows suggests a temporal shift.",
            ["Sample × taxon matrix", "Colour = relative abundance (%)", "Hover for exact value"],
        ),
        "corr_matrix": _entry(
            _WHAT_CORR_MATRIX,
            "Positively correlated families tend to grow together — they may share ecological niches or "
            "benefit from the same dietary substrates. "
            "Negatively correlated families compete. "
            "NOTE: correlations in compositional data can be spurious; interpret with caution.",
            [
                "Pearson r",
                "Red = co-occurring",
                "Blue = competing",
                "⚠ Compositional artefacts possible",
            ],
        ),
    }


def _explain_clinical(
    result: Any,
    rows: list[dict[str, Any]],
    chart_data: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    ex: dict[str, dict[str, Any]] = {}
    corr_mwt: dict[str, Any] = chart_data.get("corr_mwt") or {}
    corr_il18: dict[str, Any] = chart_data.get("corr_il18") or {}

    mwt_per_p: dict[str, dict[str, float]] = {}
    il18_per_p: dict[str, dict[str, float]] = {}
    for r in rows:
        p, tp = r["patient"], r.get("timepoint", "")
        if mwt_v := float(r.get("sixmwt") or 0):
            mwt_per_p.setdefault(p, {})[tp] = mwt_v
        if il18_v := float(r.get("il18") or 0):
            il18_per_p.setdefault(p, {})[tp] = il18_v

    def _paired(per_p: dict[str, dict[str, float]], improved: bool) -> tuple[int, int]:
        total = sum(1 for d in per_p.values() if "T0" in d and "T84" in d)
        n_imp = sum(
            1
            for d in per_p.values()
            if "T0" in d and "T84" in d and (d["T84"] > d["T0"] if improved else d["T84"] < d["T0"])
        )
        return n_imp, total

    n_mwt_imp, n_mwt_tot = _paired(mwt_per_p, improved=True)
    n_il18_imp, n_il18_tot = _paired(il18_per_p, improved=False)

    ex["clinical_mwt"] = _entry(
        _WHAT_CLINICAL_MWT,
        f"{n_mwt_imp} of {n_mwt_tot} patients improved their 6MWT distance by T84. "
        "Compare EAA vs Whey mean lines (dashed) — a larger upward shift in one group "
        "suggests that supplement may better support physical function."
        if n_mwt_tot > 0
        else "Upward-sloping lines indicate patients who improved their walking capacity during the study. "
        "Compare EAA vs Whey mean lines — a larger upward shift in one group suggests "
        "that supplement may better support physical function.",
        [
            f"{n_mwt_imp}/{n_mwt_tot} improved" if n_mwt_tot > 0 else "6MWT",
            "6MWT = 6-minute walk distance",
            "Higher = better function",
        ],
    )
    ex["clinical_il18"] = _entry(
        _WHAT_CLINICAL_IL18,
        f"{n_il18_imp} of {n_il18_tot} patients reduced their IL-18 level by T84. "
        "Compare EAA vs Whey mean lines (dashed) — a greater downward shift in one group "
        "suggests that supplement had a larger anti-inflammatory effect."
        if n_il18_tot > 0
        else "Downward-sloping lines indicate patients who reduced systemic inflammation during the study. "
        "Compare EAA vs Whey to see which supplement had a larger anti-inflammatory effect. "
        "IL-18 reductions may reflect changes in gut permeability or microbiome composition.",
        [
            f"{n_il18_imp}/{n_il18_tot} reduced" if n_il18_tot > 0 else "IL-18",
            "IL-18 pro-inflammatory marker",
            "Lower = less inflammation",
        ],
    )

    if (r_mwt := float(corr_mwt["r"])) is not None and corr_mwt.get("r") is not None:
        p_mwt = float(corr_mwt.get("p") or 1.0)
        sig_mwt = (
            "Statistically significant (p < 0.05)."
            if p_mwt < 0.05
            else "Not statistically significant (p ≥ 0.05)."
        )
        dir_mwt = (
            "A positive correlation suggests higher gut diversity associates with better physical function."
            if r_mwt > 0
            else "A negative correlation was unexpected — inspect for outlier samples."
        )
        ex["corr_mwt"] = _entry(
            _WHAT_CORR_MWT,
            f"Pearson r = {r_mwt:.2f}, p = {p_mwt:.3f}. {sig_mwt} {dir_mwt}",
            [
                f"r = {r_mwt:.2f}",
                f"p = {p_mwt:.3f}",
                "✱ Significant" if p_mwt < 0.05 else "Not significant",
                "Pearson",
            ],
        )

    if (r_il18 := float(corr_il18["r"])) is not None and corr_il18.get("r") is not None:
        p_il18 = float(corr_il18.get("p") or 1.0)
        sig_il18 = (
            "Statistically significant (p < 0.05)."
            if p_il18 < 0.05
            else "Not statistically significant (p ≥ 0.05)."
        )
        dir_il18 = (
            "A negative r suggests higher diversity associates with lower inflammation."
            if r_il18 < 0
            else "A positive r was unexpected — check for confounders or outliers."
        )
        ex["corr_il18"] = _entry(
            _WHAT_CORR_IL18,
            f"Pearson r = {r_il18:.2f}, p = {p_il18:.3f}. {sig_il18} {dir_il18}",
            [
                f"r = {r_il18:.2f}",
                f"p = {p_il18:.3f}",
                "✱ Significant" if p_il18 < 0.05 else "Not significant",
                "Pearson",
            ],
        )

    ex["taxa_clinical"] = _entry(
        _WHAT_TAXA_CLINICAL,
        "Families with starred red cells in the 6MWT column are the ones whose increase "
        "associates with improved physical function. "
        "Families with starred blue cells in the IL-18 column associate with reduced inflammation. "
        "These are the most actionable findings for future mechanistic studies.",
        ["Spearman ρ", "BH-FDR corrected", "* p<0.05 · ** p<0.01", "Δ taxon × Δ clinical"],
    )
    return ex


def _explain_longitudinal(chart_data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    long_traces = chart_data.get("longitudinal") or []
    long_changes: list[tuple[str, float]] = []
    for tr in long_traces:
        xs, ys = tr.get("x", []), tr.get("y", [])
        if "T0" in xs and "T84" in xs:
            long_changes.append(
                (str(tr.get("name", "")), float(ys[xs.index("T84")]) - float(ys[xs.index("T0")]))
            )
    long_str = "; ".join(f"{g} {d:+.3f}" for g, d in long_changes)

    return {
        "longitudinal": _entry(
            _WHAT_LONGITUDINAL,
            f"Shannon H′ change T0→T84: {long_str}. "
            "A positive value means average diversity increased; a diverging pattern between groups "
            "suggests a supplement-specific effect on the gut microbiome."
            if long_changes
            else "Compare the slope and endpoint height of EAA vs Whey lines. "
            "A diverging pattern (one group increasing while the other stays flat) would suggest "
            "the supplements have different effects on gut microbiome diversity.",
            ([f"{g}: {d:+.3f}" for g, d in long_changes[:2]] + ["Shannon H′"])
            if long_changes
            else ["Group mean per timepoint", "Shannon H′", "Summary view"],
        ),
        "lme": _entry(
            _WHAT_LME,
            "Wide confidence bands indicate high inter-individual variability — patients respond differently. "
            "Overlapping bands mean the groups cannot be reliably distinguished at this sample size. "
            "Use the metric toggle to compare Shannon, Simpson, and other alpha metrics.",
            ["Mean ± 95% CI", "Individual patient lines", "Wilcoxon p annotated"],
        ),
    }


def _explain_stats(chart_data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    ex: dict[str, dict[str, Any]] = {}
    pm: dict[str, Any] = chart_data.get("permanova") or {}

    if pm.get("rows"):
        top_name = str(pm.get("top_name") or "Unknown")
        perm_r2 = float(pm.get("top_R2") or 0.0)
        pm_rows: list[Any] = pm.get("rows") or []
        supp_row = next((r for r in pm_rows if "Group" in str(r[0]) or "Suppl" in str(r[0])), None)
        supp_p = float(supp_row[3]) if supp_row is not None else 1.0
        supp_R2 = float(supp_row[1]) if supp_row is not None else 0.0
        sig = "significant ✓" if supp_p < 0.05 else "not significant"
        ex["permanova"] = _entry(
            _WHAT_PERMANOVA,
            f"{top_name} explains the most variance (R²={perm_r2:.3f}). "
            f"Supplementation group accounts for R²={supp_R2:.3f} (p={supp_p:.3f} — {sig}). "
            "Individual identity dominating over group is expected and does not mean the intervention failed.",
            [
                f"{top_name} R²={perm_r2:.3f}",
                f"Group R²={supp_R2:.3f}",
                f"Group p={supp_p:.3f}",
                "99 permutations",
            ],
        )

    ex["stats_summary"] = _entry(
        _WHAT_STATS_SUMMARY,
        "Groups with similar T0 means are well-matched at baseline. "
        "A larger difference between T0 and T84 rows within the same group suggests "
        "the intervention had a measurable effect on that diversity metric. "
        "All five metrics are shown; consistent patterns across metrics are more reliable.",
        ["Mean ± SD per group", "All 5 alpha metrics", "Groups × timepoints"],
    )
    return ex


# ── Public API ────────────────────────────────────────────────────────────────


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
    ex.update(_explain_taxonomy(result, taxa, rows))
    ex.update(_explain_alpha(base_groups, rows))
    ex.update(_explain_beta(chart_data, taxa, patients, rows))
    ex.update(_explain_individual(chart_data, taxa, patients, timepoints, rows))
    ex.update(_explain_comparative(chart_data))
    if result.has_clinical:
        ex.update(_explain_clinical(result, rows, chart_data))
    ex.update(_explain_longitudinal(chart_data))
    ex.update(_explain_stats(chart_data))
    return ex
