/**
 * charts/beta/NMDSPlot.tsx
 * Uses pre-computed nmds1/nmds2 from demo data, falls back to PCoA coords.
 */

import { useMemo } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import { baseLayout, baseConfig, groupColor, uniqueGroups, usePlotlyDownload, withAlpha } from '@/charts/shared/usePlotly'
import { computePCoA, brayCurtis } from './distances'
import type { SampleRow } from '@/types/sample'

interface Props { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function NMDSPlot({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('nmds')
  const groups = uniqueGroups(rows)

  const coords = useMemo(() => {
    const hasNmds = rows.some(r => r.nmds1 !== undefined && Number(r.nmds1) !== 0)
    if (hasNmds) {
      return rows.map(r => ({ sample_id: String(r.sample_id), x: Number(r.nmds1 ?? 0), y: Number(r.nmds2 ?? 0) }))
    }
    return computePCoA(rows, taxa, brayCurtis).map(p => ({ ...p }))
  }, [rows, taxa])

  const traces = groups.map(g => {
    const idx = rows.map((r, i) => r.group === g ? i : -1).filter(i => i >= 0)
    const c   = groupColor(g, groups)
    return {
      type:          'scatter' as const,
      mode:          'markers' as const,
      name:          g,
      x:             idx.map(i => coords[i].x),
      y:             idx.map(i => coords[i].y),
      text:          idx.map(i => coords[i].sample_id),
      marker:        { color: withAlpha(c, 0.85), size: 10, symbol: 'diamond', line: { width: 1, color: 'white' } },
      hovertemplate: '<b>%{text}</b><extra></extra>',
    }
  })

  const layout = baseLayout({
    height: 400,
    xaxis: { ...baseLayout().xaxis, title: { text: 'NMDS1', font: { size: 11 } } },
    yaxis: { ...baseLayout().yaxis, title: { text: 'NMDS2', font: { size: 11 } } },
    showlegend: true,
  })

  return (
    <ChartCard title="NMDS — Bray-Curtis"
      subtitle="Rank-order distances · stress < 0.2 is acceptable"
      onExpand={onExpand} onDownload={download}>
      <Plot ref={plotRef} data={traces} layout={layout} config={baseConfig}
        style={{ width: '100%', height: '400px' }} useResizeHandler />
    </ChartCard>
  )
}
