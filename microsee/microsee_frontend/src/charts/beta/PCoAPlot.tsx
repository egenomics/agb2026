/**
 * charts/beta/PCoAPlot.tsx
 * Generic PCoA scatter — used for both Bray-Curtis and Jaccard.
 */

import { useMemo } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import { baseLayout, baseConfig, groupColor, uniqueGroups, usePlotlyDownload, withAlpha } from '@/charts/shared/usePlotly'
import { computePCoA, brayCurtis, jaccard } from './distances'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:      SampleRow[]
  taxa:      string[]
  distType:  'bray' | 'jaccard'
  onExpand?: () => void
}

export function PCoAPlot({ rows, taxa, distType, onExpand }: Props) {
  const label    = distType === 'bray' ? 'Bray-Curtis' : 'Jaccard (5% threshold)'
  const filename = `pcoa-${distType}`
  const { plotRef, download } = usePlotlyDownload(filename)
  const groups = uniqueGroups(rows)

  const coords = useMemo(
    () => computePCoA(rows, taxa, distType === 'bray' ? brayCurtis : jaccard),
    [rows, taxa, distType],
  )

  const pct1 = coords[0]?.pct1 ?? 0
  const pct2 = coords[0]?.pct2 ?? 0

  const traces = groups.map((g) => {
    const pts = coords.filter((_, i) => rows[i].group === g)
    const c   = groupColor(g, groups)
    return {
      type:          'scatter' as const,
      mode:          'markers' as const,
      name:          g,
      x:             pts.map(p => p.x),
      y:             pts.map(p => p.y),
      text:          pts.map(p => p.sample_id),
      marker:        { color: withAlpha(c, 0.85), size: 10, line: { width: 1, color: 'white' } },
      hovertemplate: '<b>%{text}</b><extra></extra>',
    }
  })

  const layout = baseLayout({
    height: 400,
    xaxis: { ...baseLayout().xaxis, title: { text: `PC1 (${pct1}%)`, font: { size: 11 } } },
    yaxis: { ...baseLayout().yaxis, title: { text: `PC2 (${pct2}%)`, font: { size: 11 } } },
    showlegend: true,
  })

  return (
    <ChartCard
      title={`PCoA — ${label}`}
      subtitle={`Distance = compositional dissimilarity · each dot = one sample`}
      onExpand={onExpand} onDownload={download}
    >
      <Plot ref={plotRef} data={traces} layout={layout} config={baseConfig}
        style={{ width: '100%', height: '400px' }} useResizeHandler />
    </ChartCard>
  )
}
