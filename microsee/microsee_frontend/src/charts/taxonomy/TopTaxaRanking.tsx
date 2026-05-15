/**
 * charts/taxonomy/TopTaxaRanking.tsx
 *
 * Horizontal bar chart: taxa ranked by mean abundance.
 * Replaces buildTopN() from MicroSee.html.
 */

import { useState }   from 'react'
import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseLayout, baseConfig, taxonColor,
  mean, taxonValues, usePlotlyDownload,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  taxa:     string[]
  onExpand?: () => void
}

export function TopTaxaRanking({ rows, taxa, onExpand }: Props) {
  const [topN, setTopN] = useState(9)
  const { plotRef, download } = usePlotlyDownload('top-taxa')

  const ranked = [...taxa]
    .map((t) => ({ name: t, val: mean(taxonValues(rows, t)) }))
    .sort((a, b) => b.val - a.val)
    .slice(0, topN)

  const trace = {
    type:            'bar' as const,
    orientation:     'h' as const,
    x:               ranked.map((d) => d.val),
    y:               ranked.map((d) => d.name),
    marker:          { color: ranked.map((d) => taxonColor(d.name)) },
    text:            ranked.map((d) => `${d.val.toFixed(1)}%`),
    textposition:    'outside' as const,
    hovertemplate:   '<b>%{y}</b>: %{x:.2f}%<extra></extra>',
  }

  const layout = baseLayout({
    height:    340,
    margin:    { l: 140, r: 60, t: 28, b: 40 },
    yaxis:     { autorange: 'reversed', tickfont: { size: 10 } },
    xaxis:     { title: { text: 'Mean Relative Abundance (%)', font: { size: 11 } } },
    showlegend: false,
  })

  const controls = (
    <select
      value={topN}
      onChange={(e) => setTopN(+e.target.value)}
      style={{ fontSize: '11px', padding: '3px 6px', borderRadius: '6px',
               border: '1px solid #E8C0AD', background: '#FFF8F4',
               fontFamily: 'Nunito, system-ui' }}
    >
      <option value={5}>Top 5</option>
      <option value={7}>Top 7</option>
      <option value={9}>Top 9</option>
    </select>
  )

  return (
    <ChartCard
      title="Top Taxa Ranking"
      subtitle="Mean relative abundance across all samples"
      onExpand={onExpand}
      onDownload={download}
      controls={controls}
    >
      <Plot
        ref={plotRef}
        data={[trace]}
        layout={layout}
        config={baseConfig}
        style={{ width: '100%', height: '340px' }}
        useResizeHandler
      />
    </ChartCard>
  )
}
