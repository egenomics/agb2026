/**
 * charts/taxonomy/DonutPerGroup.tsx
 *
 * 2×2 grid of donut charts, one per group.
 * Replaces buildDonut() from MicroSee.html.
 */

import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import { baseConfig, taxonColor, uniqueGroups, usePlotlyDownload, THEME } from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'
import type { Layout }    from 'plotly.js'

interface Props {
  rows:     SampleRow[]
  taxa:     string[]
  onExpand?: () => void
}

export function DonutPerGroup({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('donut-per-group')
  const groups = uniqueGroups(rows)

  // Grid: up to 2 cols
  const cols = Math.min(2, groups.length)
  const gridRows = Math.ceil(groups.length / cols)

  const traces = groups.map((g, i) => {
    const gRows = rows.filter((r) => r.group === g)
    const rawMeans = taxa.map((t) => {
      const vals = gRows.map((r) => Number(r[t] ?? 0))
      return vals.length ? vals.reduce((a, b) => a + b, 0) / vals.length : 0
    })
    const meanSum = rawMeans.reduce((a, b) => a + b, 0) || 1
    const means = rawMeans.map(v => +(v / meanSum * 100).toFixed(2))
    const row = Math.floor(i / cols) + 1
    const col = (i % cols) + 1
    return {
      type:     'pie' as const,
      name:     g,
      labels:   taxa,
      values:   means,
      hole:     0.45,
      domain:   { row: row - 1, column: col - 1 },
      marker:   { colors: taxa.map(taxonColor) },
      textinfo: 'none' as const,
      title:    {
        text:     g,
        position: 'bottom center' as const,
        font:     { size: 11, color: THEME.text, family: THEME.font },
      },
      hovertemplate: '<b>%{label}</b>: %{percent} (%{value:.1f}%)<extra></extra>',
    }
  })

  const layout: Partial<Layout> = {
    paper_bgcolor: THEME.paper,
    plot_bgcolor:  THEME.bg,
    font:          { family: THEME.font, color: THEME.text, size: 11 },
    height:        gridRows * 240,
    margin:        { l: 16, r: 16, t: 16, b: 30 },
    showlegend:    true,
    legend:        { orientation: 'h', y: -0.08, x: 0.5, xanchor: 'center', font: { size: 9 } },
    grid:          { rows: gridRows, columns: cols, pattern: 'independent' as const },
  }

  return (
    <ChartCard
      title="Composition per Group"
      subtitle="Average taxon proportions by group · hover for values"
      onExpand={onExpand}
      onDownload={download}
    >
      <Plot
        ref={plotRef}
        data={traces}
        layout={layout}
        config={baseConfig}
        style={{ width: '100%', height: `${gridRows * 240}px` }}
        useResizeHandler
      />
    </ChartCard>
  )
}
