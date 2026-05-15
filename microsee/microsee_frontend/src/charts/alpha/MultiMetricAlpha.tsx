/**
 * charts/alpha/MultiMetricAlpha.tsx
 *
 * Grouped bar chart comparing Observed OTUs, Chao1, Pielou J′ per group.
 * Replaces buildAlpha2() from MicroSee.html.
 */

import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseConfig, groupColor, mean, uniqueGroups,
  usePlotlyDownload, withAlpha, THEME,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  taxa:     string[]
  onExpand?: () => void
}

function observedTaxa(row: SampleRow, taxa: string[]): number {
  return taxa.filter((t) => Number(row[t] ?? 0) > 0).length
}

function chao1(row: SampleRow, taxa: string[]): number {
  const obs = observedTaxa(row, taxa)
  // F1 = singletons, F2 = doubletons (by rank within sample, since values are proportions)
  const sorted = taxa.map(t => Number(row[t] ?? 0)).filter(v => v > 0).sort((a, b) => a - b)
  const f1 = sorted.length >= 1 ? 1 : 0          // rarest taxon counts as F1
  const f2 = sorted.length >= 2 ? 1 : 0          // second rarest counts as F2
  // Chao1 = S_obs + F1² / (2·F2)
  return obs + (f1 * f1) / (2 * Math.max(f2, 1))
}

function pielou(row: SampleRow, taxa: string[]): number {
  const ps  = taxa.map((t) => Number(row[t] ?? 0))
  const tot = ps.reduce((a, b) => a + b, 0) || 1
  const p   = ps.map((v) => v / tot)
  const H   = -p.reduce((s, pi) => (pi > 0 ? s + pi * Math.log(pi) : s), 0)
  const S   = taxa.filter((t) => Number(row[t] ?? 0) > 0).length
  const maxH = Math.log(S || 1)
  return maxH > 0 ? H / maxH : 0
}

const METRICS = [
  { key: 'observed', label: 'Observed Taxa', fn: observedTaxa },
  { key: 'chao1',    label: 'Chao1 est.',    fn: chao1        },
  { key: 'pielou',   label: "Pielou J′",     fn: pielou       },
]

export function MultiMetricAlpha({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('multi-metric-alpha')
  const groups = uniqueGroups(rows)

  const traces = groups.map((g) => {
    const gd = rows.filter((r) => r.group === g)
    const c  = groupColor(g, groups)
    return {
      type:          'bar'       as const,
      name:          g,
      x:             METRICS.map((m) => m.label),
      y:             METRICS.map((m) => mean(gd.map((r) => m.fn(r, taxa)))),
      marker:        { color: withAlpha(c, 0.8), line: { color: c, width: 1.5 } },
      hovertemplate: `<b>${g}</b><br>%{x}: %{y:.2f}<extra></extra>`,
    }
  })

  return (
    <ChartCard
      title="Multi-Metric Alpha Diversity"
      subtitle="Observed Taxa · Chao1 richness estimate · Pielou J′ evenness"
      onExpand={onExpand}
      onDownload={download}
    >
      <Plot
        ref={plotRef}
        data={traces}
        layout={{
          barmode:        'group',
          paper_bgcolor:  THEME.paper,
          plot_bgcolor:   THEME.bg,
          font:           { family: THEME.font, color: THEME.text, size: 11 },
          height:         320,
          margin:         { l: 52, r: 16, t: 24, b: 60 },
          xaxis:          { tickfont: { size: 10 }, gridcolor: 'transparent' },
          yaxis:          { gridcolor: THEME.grid },
          legend:         { orientation: 'h', y: -0.3, x: 0, font: { size: 10 } },
          showlegend:     true,
        }}
        config={baseConfig}
        style={{ width: '100%', height: '320px' }}
        useResizeHandler
      />
    </ChartCard>
  )
}
