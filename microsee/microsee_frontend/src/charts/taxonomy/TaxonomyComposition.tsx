/**
 * charts/taxonomy/TaxonomyComposition.tsx
 *
 * Stacked bar chart: one bar per sample, one trace per taxon.
 * Replaces buildComp() from MicroSee.html.
 */

import Plot from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseLayout, baseConfig, taxonColor, normaliseRows,
  usePlotlyDownload, THEME,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  taxa:     string[]
  onExpand?: () => void
}

export function TaxonomyComposition({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('taxonomy-composition')
  const norm = normaliseRows(rows, taxa)

  const traces = taxa.map((t) => ({
    type:        'bar'  as const,
    name:        t,
    x:           norm.map((r) => r.sample_id as string),
    y:           norm.map((r) => Number(r[t] ?? 0)),
    marker:      { color: taxonColor(t) },
    hovertemplate: `<b>%{x}</b><br>${t}: %{y:.1f}%<extra></extra>`,
  }))

  const layout = baseLayout({
    barmode:     'stack',
    height:      380,
    margin:      { l: 52, r: 16, t: 28, b: 90 },
    xaxis:       { tickangle: -40, tickfont: { size: 9 }, title: undefined },
    yaxis:       { title: { text: 'Relative Abundance (%)', font: { size: 11 } } },
    showlegend:  true,
    legend:      { orientation: 'h', y: -0.35, x: 0, font: { size: 9 } },
  })

  return (
    <ChartCard
      title="Taxonomic Composition"
      subtitle="Relative abundance % per sample · toggle taxa in sidebar"
      onExpand={onExpand}
      onDownload={download}
    >
      <Plot
        ref={plotRef}
        data={traces}
        layout={layout}
        config={baseConfig}
        style={{ width: '100%', height: '380px' }}
        useResizeHandler
      />
    </ChartCard>
  )
}
