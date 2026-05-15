import { useState } from 'react'
import { useAppStore }       from '@/store/appStore'
import { StatsBar }          from '@/components/StatsBar'
import { SectionNav, useSectionNav } from '@/components/SectionNav'
import { FlashCard, useFlashCard } from '@/components/FlashCard'
import { mean }              from '@/charts/shared/usePlotly'
import styles from './ResultsPage.module.css'

// Taxonomy charts
import { TaxonomyComposition } from '@/charts/taxonomy/TaxonomyComposition'
import { TopTaxaRanking }      from '@/charts/taxonomy/TopTaxaRanking'
import { DonutPerGroup }       from '@/charts/taxonomy/DonutPerGroup'
import { Sunburst }            from '@/charts/taxonomy/Sunburst'

// Alpha diversity charts
import { AlphaStrip }       from '@/charts/alpha/AlphaStrip'
import { BoxPlot }          from '@/charts/alpha/BoxPlot'
import { ViolinPlot }       from '@/charts/alpha/ViolinPlot'
import { RarefactionCurves } from '@/charts/alpha/RarefactionCurves'
import { MultiMetricAlpha } from '@/charts/alpha/MultiMetricAlpha'

// Beta charts
import { PCoAPlot }       from '@/charts/beta/PCoAPlot'
import { NMDSPlot }       from '@/charts/beta/NMDSPlot'
import { DendrogramPlot } from '@/charts/beta/DendrogramPlot'
import { DeltaHeatmap }   from '@/charts/beta/DeltaHeatmap'

// Individual charts
import { NMDSTrajectories } from '@/charts/individual/NMDSTrajectories'
import { PairedSlopegraph, StabilityBar, DiversityRank } from '@/charts/individual/IndividualCharts'
import { PatientRadar, FacetedComposition } from '@/charts/individual/PatientCharts'

// Comparative charts
import { DifferentialAbundance, VolcanoPlot, AbundanceHeatmap, CorrelationMatrix } from '@/charts/comparative/ComparativeCharts'

// Clinical charts
import { ClinicalSlopegraph, ClinicalCorrelation } from '@/charts/clinical/ClinicalCharts'

// Longitudinal + Stats
import { LongitudinalChart, PERMANOVATable, DiversitySummary } from '@/charts/longitudinal/LongitudinalStats'

// ── FlashCard info registry ──────────────────────────────────────────────────

const CHART_INFO: Record<string, { title: string; what: string; insight: string }> = {
  comp: {
    title:   'Taxonomic Composition',
    what:    'Stacked bars showing the relative abundance (%) of each bacterial family per sample. Each colour is one family; bar height is proportional to its share of the community.',
    insight: 'Look for families that dominate one group but not another. Samples that look very different from their group peers could be low-diversity outliers worth investigating.',
  },
  topn: {
    title:   'Top Taxa Ranking',
    what:    'Mean relative abundance of the most common bacterial families across all samples, ranked highest to lowest.',
    insight: 'In healthy older adults, Bacteroidaceae and Lachnospiraceae typically account for 50–60% combined. Large differences between groups may indicate a treatment effect.',
  },
  donut: {
    title:   'Composition per Group',
    what:    'Doughnut charts showing the average relative abundance of each family within each group.',
    insight: 'If the doughnuts look nearly identical across groups, the microbiome composition was not significantly changed by the intervention.',
  },
  sunburst: {
    title:   'Sunburst — Group → Taxa',
    what:    'Two-ring sunburst: the inner ring divides by group, the outer ring subdivides each group by family. Arc width equals relative abundance.',
    insight: 'Wedges of the same family at the same angular size across groups means balanced composition. Hover over an arc to see exact values.',
  },
  'pcoa-bray': {
    title:   'PCoA — Bray-Curtis',
    what:    'Principal Coordinates Analysis of Bray-Curtis distance between samples. Points close together have more similar microbial communities.',
    insight: 'Look for group separation along the first two axes. Tight clusters indicate consistent community structure within a group.',
  },
  'pcoa-jaccard': {
    title:   'PCoA — Jaccard',
    what:    'Principal Coordinates Analysis of Jaccard distance, which compares sample similarity based on shared taxa presence/absence.',
    insight: 'This view highlights composition changes independent of abundance. Differences here may reflect taxa gain/loss between groups.',
  },
  nmds: {
    title:   'NMDS Plot',
    what:    'Non-metric multidimensional scaling of beta diversity distances, optimized to preserve rank ordering of sample distances.',
    insight: 'Stress values low enough mean the 2D plot is a faithful representation. Closer points still indicate more similar communities.',
  },
  dendro: {
    title:   'Dendrogram',
    what:    'Hierarchical clustering tree of sample distances. Branch length reflects similarity between communities.',
    insight: 'Samples that join together early are more similar. Look for group-level branching patterns for treatment effects.',
  },
  delta: {
    title:   'Delta Heatmap',
    what:    'Heatmap showing pairwise differences between samples. Darker colors indicate larger beta diversity distance.',
    insight: 'This is useful for spotting strong outliers or subgroups with unusually large within-group differences.',
  },
  traj: {
    title:   'Trajectories',
    what:    'Individual patient trajectories across time in NMDS space, showing how each subject moves through microbiome composition space.',
    insight: 'Consistent movement in one direction across patients suggests a shared treatment impact on community structure.',
  },
  slope: {
    title:   'Paired Slopegraph',
    what:    'Each line connects a patient’s diversity score at baseline and follow-up. Dashed lines show group means.',
    insight: 'Parallel lines indicate consistent change; crossing lines suggest heterogeneous responses.',
  },
  stability: {
    title:   'Stability Bar',
    what:    'Bar plot quantifying within-patient stability over time, comparing community dissimilarity at baseline and follow-up.',
    insight: 'Smaller bars mean more stable microbiomes. Larger values may indicate strong changes from intervention or recovery.',
  },
  divrank: {
    title:   'Diversity Rank Plot',
    what:    'Samples are ranked by diversity score and colored by group. Shapes distinguish baseline vs follow-up samples.',
    insight: 'Group separation in rank space helps identify whether one group consistently has higher diversity than another.',
  },
  radar: {
    title:   'Patient Radar',
    what:    'Radar charts compare multiple taxa profiles for selected patients or groups, showing relative abundance patterns.',
    insight: 'Distinct radar shapes indicate different dominant taxa ecosystems between patients.',
  },
  faceted: {
    title:   'Faceted Composition',
    what:    'Multiple composition plots shown side-by-side for each patient or group, highlighting differences in relative abundance.',
    insight: 'Use this view to compare several individual profiles at once and spot consistent format changes.',
  },
  diff: {
    title:   'Differential Abundance',
    what:    'Scatter plot of taxa changes between groups, identifying those with the largest effect sizes and significance.',
    insight: 'Taxa far from zero are the most differentially abundant and may be responsible for group differences.',
  },
  volcano: {
    title:   'Volcano Plot',
    what:    'Displays log fold change versus significance for each taxon. Top hits are both large and statistically significant.',
    insight: 'Taxa in the top left/right corners are strong candidates for biomarkers or treatment-responsive species.',
  },
  heatmap: {
    title:   'Abundance Heatmap',
    what:    'Heatmap of taxa abundance across samples, clustered to reveal patterns and group structure.',
    insight: 'Clusters of similar rows or columns indicate taxa and samples with correlated abundance patterns.',
  },
  corrmat: {
    title:   'Correlation Matrix',
    what:    'Pairwise taxon correlations shown as a matrix. Strong positive or negative correlations highlight co-occurrence patterns.',
    insight: 'Highly correlated taxa may share ecological niches or respond similarly to the same conditions.',
  },
  mwt: {
    title:   'Clinical Slopegraph — 6MWT',
    what:    'Changes in 6-minute walk test performance for each patient, plotted as baseline-to-follow-up slopes.',
    insight: 'Consistent upward slopes mean functional improvement across patients.',
  },
  il18: {
    title:   'Clinical Slopegraph — IL-18',
    what:    'Changes in IL-18 cytokine concentration for each patient, plotted as baseline-to-follow-up slopes.',
    insight: 'Shifts in IL-18 may reflect inflammatory response differences between visits.',
  },
  'corr-mwt': {
    title:   'Clinical Correlation — 6MWT',
    what:    'Scatter plot correlating Shannon or Simpson diversity with 6-minute walk test performance.',
    insight: 'A strong correlation suggests diversity may be linked to functional mobility.',
  },
  'corr-il18': {
    title:   'Clinical Correlation — IL-18',
    what:    'Scatter plot correlating Shannon or Simpson diversity with IL-18 cytokine levels.',
    insight: 'This can reveal whether microbial diversity associates with inflammation markers.',
  },
  long: {
    title:   'Longitudinal Diversity',
    what:    'Tracks each patient’s diversity over time to show how microbiome diversity changes across visits.',
    insight: 'Diverging trajectories may indicate different patient responses to treatment or time.',
  },
  permanova: {
    title:   'PERMANOVA',
    what:    'Statistical test results for group, timepoint, and individual differences based on Bray-Curtis distance.',
    insight: 'Significant p-values indicate that sample grouping explains a meaningful amount of community variation.',
  },
  divsum: {
    title:   'Diversity Summary',
    what:    'Summary table with mean and standard deviation diversity scores per group.',
    insight: 'Compare the means and spread to see which groups have higher and more stable diversity.',
  },
  strip: {
    title:   'Alpha Diversity — Strip Chart',
    what:    'Each dot is one sample. The horizontal bar is the group mean. Shannon entropy captures both richness (how many taxa) and evenness (how equal their abundances are).',
    insight: 'A low Shannon score suggests the community is dominated by one or a few species. Compare groups to assess whether the intervention changed diversity.',
  },

  box: {
    title:   'Box Plot — Alpha Diversity',
    what:    'Box shows the interquartile range (25th–75th percentile). Whiskers extend to 1.5× IQR. Individual dots show all samples.',
    insight: 'Overlapping boxes between groups suggest no significant diversity difference. Non-overlapping boxes warrant a Kruskal–Wallis test.',
  },
  violin: {
    title:   'Violin Plot — Alpha Diversity',
    what:    'Width of the violin at each value shows how many samples have that diversity score. The box inside shows median and IQR.',
    insight: 'A wide violin at high Shannon values means many diverse samples. Flat, thin violins mean the group has consistent but low diversity.',
  },
  rar: {
    title:   'Rarefaction Curves',
    what:    'Each line shows how many bacterial species are detected as you sample more reads. Curves that plateau indicate the community has been adequately sequenced.',
    insight: 'Steep curves that have not plateaued suggest more sequencing depth is needed. Bold lines are group means.',
  },
  multimet: {
    title:   'Multi-Metric Alpha Diversity',
    what:    'Compares three complementary alpha diversity metrics: Observed OTUs (raw count), Chao1 (estimate accounting for rare species), and Pielou J′ (evenness from 0 to 1).',
    insight: 'Using multiple metrics guards against bias from any single index. High Pielou J′ (near 1) confirms even distribution with no single family dominating.',
  },
}

// ── Main page ────────────────────────────────────────────────────────────────

export function ResultsPage() {
  const { result, filteredRows, filteredTaxa } = useAppStore()
  const { activeSection, setActiveSection, isVisible }  = useSectionNav()
  const fc = useFlashCard()
  const [fcKey, setFcKey] = useState<string>('')

  if (!result) return null

  const rows = filteredRows()
  const taxa = filteredTaxa()

  const openFC = (key: string) => { setFcKey(key); fc.show() }

  // Dynamic pills for flashcard
  const fcPills = (): string[] => {
    if (!fcKey) return []
    if (['strip','box','violin','rar','multimet'].includes(fcKey)) {
      const groups = [...new Set(rows.map((r) => r.group))].sort()
      return groups.map((g) => {
        const v = mean(rows.filter((r) => r.group === g).map((r) => Number(r.shannon)))
        return `${g}: H′ ${v.toFixed(3)}`
      })
    }
    if (['comp','topn','donut','sunburst'].includes(fcKey)) {
      return [`${taxa.length} families shown`, `${rows.length} samples`]
    }
    return []
  }

  return (
    <div className={styles.page}>
      <StatsBar />
      <SectionNav activeSection={activeSection} onChange={setActiveSection} />

      <div className={styles.content}>

        {/* ── TAXONOMY ── */}
        {isVisible('taxonomy') && (
          <section id="sec-taxonomy" className={styles.section}>
            <div className={styles.secHeader}>🔬 Taxonomy</div>
            <div className={styles.grid2}>
              <TaxonomyComposition rows={rows} taxa={taxa} onExpand={() => openFC('comp')} />
              <TopTaxaRanking      rows={rows} taxa={taxa} onExpand={() => openFC('topn')} />
            </div>
            <div className={styles.grid2}>
              <DonutPerGroup rows={rows} taxa={taxa} onExpand={() => openFC('donut')} />
              <Sunburst      rows={rows} taxa={taxa} onExpand={() => openFC('sunburst')} />
            </div>
          </section>
        )}

        {/* ── ALPHA ── */}
        {isVisible('alpha') && (
          <section id="sec-alpha" className={styles.section}>
            <div className={styles.secHeader}>
              📊 Alpha Diversity
            </div>
            <div className={styles.grid2}>
              <AlphaStrip rows={rows} onExpand={() => openFC('strip')} />
              <BoxPlot    rows={rows} onExpand={() => openFC('box')} />
            </div>
            <div className={styles.grid1}>
              <ViolinPlot rows={rows} onExpand={() => openFC('violin')} />
            </div>
            <div className={styles.grid2}>
              <RarefactionCurves rows={rows} taxa={taxa} onExpand={() => openFC('rar')} />
              <MultiMetricAlpha  rows={rows} taxa={taxa} onExpand={() => openFC('multimet')} />
            </div>
          </section>
        )}

        {/* ── BETA ── */}
        {isVisible('beta') && (
          <section id="sec-beta" className={styles.section}>
            <div className={styles.secHeader}>🌐 Beta Diversity</div>
            <div className={styles.grid2}>
              <PCoAPlot rows={rows} taxa={taxa} distType="bray"    onExpand={() => openFC('pcoa-bray')} />
              <PCoAPlot rows={rows} taxa={taxa} distType="jaccard" onExpand={() => openFC('pcoa-jaccard')} />
            </div>
            <div className={styles.grid1}>
              <NMDSPlot rows={rows} taxa={taxa} onExpand={() => openFC('nmds')} />
            </div>
            <div className={styles.grid1}>
              <DendrogramPlot rows={rows} taxa={taxa} onExpand={() => openFC('dendro')} />
            </div>
            <div className={styles.grid1}>
              <DeltaHeatmap rows={rows} taxa={taxa} onExpand={() => openFC('delta')} />
            </div>
          </section>
        )}

        {/* ── INDIVIDUAL ── */}
        {isVisible('individual') && (
          <section id="sec-individual" className={styles.section}>
            <div className={styles.secHeader}>👤 Individual Analysis</div>
            <div className={styles.grid2}>
              <NMDSTrajectories  rows={rows} taxa={taxa} onExpand={() => openFC('traj')} />
              <PairedSlopegraph  rows={rows} onExpand={() => openFC('slope')} />
            </div>
            <div className={styles.grid2}>
              <StabilityBar  rows={rows} taxa={taxa} onExpand={() => openFC('stability')} />
              <DiversityRank rows={rows} onExpand={() => openFC('divrank')} />
            </div>
            <div className={styles.grid1}>
              <PatientRadar rows={rows} taxa={taxa} onExpand={() => openFC('radar')} />
            </div>
            <div className={styles.grid1}>
              <FacetedComposition rows={rows} taxa={taxa} onExpand={() => openFC('faceted')} />
            </div>
          </section>
        )}

        {/* ── COMPARATIVE ── */}
        {isVisible('compare') && (
          <section id="sec-compare" className={styles.section}>
            <div className={styles.secHeader}>⚖️ Comparative</div>
            <div className={styles.grid2}>
              <DifferentialAbundance rows={rows} taxa={taxa} onExpand={() => openFC('diff')} />
              <VolcanoPlot           rows={rows} taxa={taxa} onExpand={() => openFC('volcano')} />
            </div>
            <div className={styles.grid1}>
              <AbundanceHeatmap rows={rows} taxa={taxa} onExpand={() => openFC('heatmap')} />
            </div>
            <div className={styles.grid1}>
              <CorrelationMatrix rows={rows} taxa={taxa} onExpand={() => openFC('corrmat')} />
            </div>
          </section>
        )}

        {/* ── CLINICAL ── */}
        {isVisible('clinical') && (
          <section id="sec-clinical" className={styles.section}>
            <div className={styles.secHeader}>💊 Clinical Outcomes</div>
            {!result.has_clinical && (
              <div className={styles.clinNote}>
                No clinical data found. Add <code>sixmwt</code> and <code>il18</code> columns
                to metadata.tsv to enable this section.
              </div>
            )}
            <div className={styles.grid2}>
              <ClinicalSlopegraph rows={rows} field="sixmwt" label="6-Min Walk Test" unit="m"     onExpand={() => openFC('mwt')} />
              <ClinicalSlopegraph rows={rows} field="il18"   label="IL-18 Cytokine"  unit="pg/mL" onExpand={() => openFC('il18')} />
            </div>
            <div className={styles.grid2}>
              <ClinicalCorrelation rows={rows} yField="sixmwt" yLabel="6MWT"  yUnit="m"     onExpand={() => openFC('corr-mwt')} />
              <ClinicalCorrelation rows={rows} yField="il18"   yLabel="IL-18" yUnit="pg/mL" onExpand={() => openFC('corr-il18')} />
            </div>
          </section>
        )}

        {/* ── LONGITUDINAL ── */}
        {isVisible('longitudinal') && (
          <section id="sec-longitudinal" className={styles.section}>
            <div className={styles.secHeader}>📈 Longitudinal</div>
            <div className={styles.grid1}>
              <LongitudinalChart rows={rows} onExpand={() => openFC('long')} />
            </div>
          </section>
        )}

        {/* ── STATS ── */}
        {isVisible('stats') && (
          <section id="sec-stats" className={styles.section}>
            <div className={styles.secHeader}>📋 Statistics</div>
            <div className={styles.grid2}>
              <PERMANOVATable  rows={rows} taxa={taxa} onExpand={() => openFC('permanova')} />
              <DiversitySummary rows={rows} onExpand={() => openFC('divsum')} />
            </div>
          </section>
        )}

      </div>

      {/* ── FlashCard modal ── */}
      {fcKey && CHART_INFO[fcKey] && (
        <FlashCard
          isOpen={fc.open}
          onClose={fc.close}
          info={{ ...CHART_INFO[fcKey], pills: fcPills() }}
        >
          {fcKey === 'comp'     && <TaxonomyComposition rows={rows} taxa={taxa} />}
          {fcKey === 'topn'     && <TopTaxaRanking rows={rows} taxa={taxa} />}
          {fcKey === 'donut'    && <DonutPerGroup rows={rows} taxa={taxa} />}
          {fcKey === 'sunburst' && <Sunburst rows={rows} taxa={taxa} />}
          {fcKey === 'strip'    && <AlphaStrip rows={rows} />}
          {fcKey === 'box'      && <BoxPlot    rows={rows} />}
          {fcKey === 'violin'   && <ViolinPlot rows={rows} />}
          {fcKey === 'rar'      && <RarefactionCurves rows={rows} taxa={taxa} />}
          {fcKey === 'multimet' && <MultiMetricAlpha  rows={rows} taxa={taxa} />}
          {fcKey === 'pcoa-bray' && <PCoAPlot rows={rows} taxa={taxa} distType="bray" />}
          {fcKey === 'pcoa-jaccard' && <PCoAPlot rows={rows} taxa={taxa} distType="jaccard" />}
          {fcKey === 'nmds'     && <NMDSPlot rows={rows} taxa={taxa} />}
          {fcKey === 'dendro'   && <DendrogramPlot rows={rows} taxa={taxa} />}
          {fcKey === 'delta'    && <DeltaHeatmap rows={rows} taxa={taxa} />}
          {fcKey === 'traj'     && <NMDSTrajectories rows={rows} taxa={taxa} />}
          {fcKey === 'slope'    && <PairedSlopegraph rows={rows} />}
          {fcKey === 'stability' && <StabilityBar rows={rows} taxa={taxa} />}
          {fcKey === 'divrank'  && <DiversityRank rows={rows} />}
          {fcKey === 'radar'    && <PatientRadar rows={rows} taxa={taxa} />}
          {fcKey === 'faceted'  && <FacetedComposition rows={rows} taxa={taxa} />}
          {fcKey === 'diff'     && <DifferentialAbundance rows={rows} taxa={taxa} />}
          {fcKey === 'volcano'  && <VolcanoPlot rows={rows} taxa={taxa} />}
          {fcKey === 'heatmap'  && <AbundanceHeatmap rows={rows} taxa={taxa} />}
          {fcKey === 'corrmat'  && <CorrelationMatrix rows={rows} taxa={taxa} />}
          {fcKey === 'mwt'      && <ClinicalSlopegraph rows={rows} field="sixmwt" label="6-Min Walk Test" unit="m" />}
          {fcKey === 'il18'     && <ClinicalSlopegraph rows={rows} field="il18" label="IL-18 Cytokine" unit="pg/mL" />}
          {fcKey === 'corr-mwt' && <ClinicalCorrelation rows={rows} yField="sixmwt" yLabel="6MWT" yUnit="m" />}
          {fcKey === 'corr-il18' && <ClinicalCorrelation rows={rows} yField="il18" yLabel="IL-18" yUnit="pg/mL" />}
          {fcKey === 'long'     && <LongitudinalChart rows={rows} />}
          {fcKey === 'permanova' && <PERMANOVATable rows={rows} taxa={taxa} />}
          {fcKey === 'divsum'   && <DiversitySummary rows={rows} />}
        </FlashCard>
      )}
    </div>
  )
}
