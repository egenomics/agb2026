/**
 * charts/individual/NMDSTrajectories.tsx
 * Arrows connecting T0 → T84 per patient on NMDS space.
 */

import { useMemo } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import { baseLayout, baseConfig, groupColor, usePlotlyDownload, withAlpha } from '@/charts/shared/usePlotly'
import { computePCoA, brayCurtis } from '@/charts/beta/distances'
import type { SampleRow } from '@/types/sample'

interface Props { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function NMDSTrajectories({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('nmds-trajectories')

  const { traces } = useMemo(() => {
    const hasNmds = rows.some(r => Number.isFinite(Number(r.nmds1)) && Number(r.nmds1) !== 0)
    let coords: { sample_id: string; x: number; y: number }[]
    if (hasNmds) {
      coords = rows.map(r => ({ sample_id: String(r.sample_id), x: Number(r.nmds1 ?? 0), y: Number(r.nmds2 ?? 0) }))
    } else {
      coords = computePCoA(rows, taxa, brayCurtis)
    }
    const coordMap = Object.fromEntries(coords.map(c => [c.sample_id, c]))
    const pid = (r: SampleRow) => String(r.patient ?? r.sample_id).replace(/_T\d+$/, '')
    const baseG = (r: SampleRow) => String(r.base_group ?? r.group).replace(/_T\d+$/, '')
    const allBase = [...new Set(rows.map(baseG))].sort()
    const patients = [...new Set(rows.map(pid))]

    const traces: object[] = []

    patients.forEach(p => {
      const r0  = rows.find(r => pid(r) === p && Number(r.time) === 0)
      const r84 = rows.filter(r => pid(r) === p && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
      if (!r0 || !r84) return
      const bg = baseG(r0)
      const c  = groupColor(bg, allBase)
      const c0 = coordMap[String(r0.sample_id)]
      const c1 = coordMap[String(r84.sample_id)]
      if (!c0 || !c1) return

      // Line
      traces.push({
        type: 'scatter', mode: 'lines',
        x: [c0.x, c1.x], y: [c0.y, c1.y],
        line: { color: withAlpha(c, 0.6), width: 1.8 },
        showlegend: false, hoverinfo: 'skip',
      })
      // T0 circle
      traces.push({
        type: 'scatter', mode: 'markers',
        x: [c0.x], y: [c0.y],
        marker: { color: withAlpha(c, 0.5), size: 8, symbol: 'circle', line: { color: c, width: 1 } },
        name: `${p} T0`, showlegend: false,
        hovertemplate: `<b>${p} T0</b><extra></extra>`,
      })
      // T84 filled dot
      traces.push({
        type: 'scatter', mode: 'markers',
        x: [c1.x], y: [c1.y],
        marker: { color: c, size: 9, symbol: 'circle', line: { color: 'white', width: 1 } },
        name: `${p} T84`, showlegend: false,
        hovertemplate: `<b>${p} T84</b><extra></extra>`,
      })
    })

    // Legend
    allBase.forEach(bg => {
      traces.push({
        type: 'scatter', mode: 'markers+lines',
        x: [null], y: [null],
        marker: { color: groupColor(bg, allBase), size: 8 },
        line: { color: groupColor(bg, allBase), width: 2 },
        name: bg, showlegend: true,
      })
    })

    return { traces }
  }, [rows, taxa])

  const layout = baseLayout({
    height: 420, showlegend: true,
    xaxis: { ...baseLayout().xaxis, title: { text: 'NMDS1 / PC1', font: { size: 11 } } },
    yaxis: { ...baseLayout().yaxis, title: { text: 'NMDS2 / PC2', font: { size: 11 } } },
  })

  return (
    <ChartCard title="NMDS — Individual Trajectories"
      subtitle="○ = T0 start · ● = T84 end · short arrow = stable microbiome"
      onExpand={onExpand} onDownload={download}>
      <Plot ref={plotRef} data={traces as never[]} layout={layout} config={baseConfig}
        style={{ width: '100%', height: '420px' }} useResizeHandler />
    </ChartCard>
  )
}
