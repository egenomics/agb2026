/**
 * charts/beta/DendrogramPlot.tsx
 * Hierarchical clustering dendrogram using average-linkage on Bray-Curtis.
 */

import { useMemo } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import { baseConfig, groupColor, uniqueGroups, usePlotlyDownload, THEME } from '@/charts/shared/usePlotly'
import { buildDistMatrix, brayCurtis } from './distances'
import type { SampleRow } from '@/types/sample'

interface Props { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

/** Average-linkage hierarchical clustering — returns dendrogram line segments (right-to-left) */
function averageLinkage(mat: number[][], labels: string[]) {
  const n = labels.length
  const clusters: number[][] = labels.map((_, i) => [i])
  const dist = mat.map(r => [...r])
  const active = new Set(Array.from({ length: n }, (_, i) => i))
  const segments: { x: [number, number, number, number]; y: [number, number, number, number] }[] = []
  const heights: number[] = new Array(n).fill(0)

  let nextId = n
  const clusterPos: Record<number, number> = {}
  labels.forEach((_, i) => { clusterPos[i] = i })

  for (let step = 0; step < n - 1; step++) {
    let minD = Infinity, ci = -1, cj = -1
    const arr = [...active]
    for (let a = 0; a < arr.length; a++) {
      for (let b = a + 1; b < arr.length; b++) {
        if (dist[arr[a]][arr[b]] < minD) { minD = dist[arr[a]][arr[b]]; ci = arr[a]; cj = arr[b] }
      }
    }
    if (ci < 0) break
    const posI = clusterPos[ci], posJ = clusterPos[cj]
    const posNew = (posI + posJ) / 2
    clusterPos[nextId] = posNew
    const heightNew = minD
    // Right-to-left: x = distance (axis reversed), y = position
    segments.push({ x: [heights[ci] ?? 0, heightNew, heightNew, heights[cj] ?? 0], y: [posI, posI, posJ, posJ] })
    heights[nextId] = heightNew

    active.delete(ci); active.delete(cj); active.add(nextId)
    dist[nextId] = []
    dist.forEach((_, k) => { dist[nextId][k] = dist[k][nextId] = 0 })
    for (const k of active) {
      if (k === nextId) continue
      const ni = clusters[ci]?.length ?? 1, nj = clusters[cj]?.length ?? 1
      dist[nextId][k] = dist[k][nextId] = (dist[ci][k] * ni + dist[cj][k] * nj) / (ni + nj)
    }
    clusters[nextId] = [...(clusters[ci] ?? []), ...(clusters[cj] ?? [])]
    nextId++
  }

  const leafOrder = Object.entries(clusterPos)
    .filter(([id]) => +id < n)
    .sort((a, b) => a[1] - b[1])
    .map(([id]) => labels[+id])

  return { segments, leafOrder }
}

export function DendrogramPlot({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('dendrogram')
  const groups = uniqueGroups(rows)
  const labels = rows.map(r => String(r.sample_id))
  const sid2grp = Object.fromEntries(rows.map(r => [String(r.sample_id), r.group]))

  const { segments, leafOrder } = useMemo(() => {
    const mat = buildDistMatrix(rows, taxa, brayCurtis)
    return averageLinkage(mat, labels)
  }, [rows, taxa])

  const n = leafOrder.length
  const dynamicHeight = Math.max(380, n * 22 + 80)

  const traces = [
    // Dendrogram lines (x = distance, y = leaf position)
    ...segments.map(s => ({
      type:       'scatter' as const,
      mode:       'lines' as const,
      x:          s.x, y: s.y,
      line:       { color: '#C4A0A8', width: 1.5 },
      showlegend: false,
      hoverinfo:  'skip' as const,
    })),
    // Leaf markers at x=0 (rightmost due to reversed axis) — labels shown via yticks on right
    {
      type:          'scatter' as const,
      mode:          'markers' as const,
      x:             leafOrder.map(() => 0),
      y:             leafOrder.map((_, i) => i),
      marker:        { color: leafOrder.map(s => groupColor(sid2grp[s] ?? '', groups)), size: 8 },
      hovertemplate: leafOrder.map(s => `${s}<extra></extra>`),
      showlegend:    false,
    },
    // Legend-only traces per group
    ...groups.map(g => ({
      type:       'scatter' as const,
      mode:       'markers' as const,
      x:          [null as unknown as number], y: [null as unknown as number],
      name:       g, marker: { color: groupColor(g, groups), size: 8 },
      showlegend: true,
    })),
  ]

  return (
    <ChartCard title="Hierarchical Clustering Dendrogram"
      subtitle="Average-linkage on Bray-Curtis · leaves on right, root on left · most similar merge last"
      onExpand={onExpand} onDownload={download}>
      <Plot ref={plotRef} data={traces}
        layout={{
          autosize: true,
          paper_bgcolor: THEME.paper, plot_bgcolor: THEME.bg,
          font: { family: THEME.font, color: THEME.text, size: 10 },
          height: dynamicHeight,
          margin: { l: 24, r: 140, t: 28, b: 48 },
          xaxis: {
            title:     { text: 'Distance', font: { size: 11 } },
            zeroline:  false,
            gridcolor: THEME.grid,
            autorange: 'reversed' as const,
          },
          yaxis: {
            side:      'right' as const,
            tickmode:  'array' as const,
            tickvals:  leafOrder.map((_, i) => i),
            ticktext:  leafOrder.map(s => s.replace('_', '·')),
            tickfont:  { size: 9, color: leafOrder.map(s => groupColor(sid2grp[s] ?? '', groups)) },
            zeroline:  false,
            gridcolor: 'transparent',
            showgrid:  false,
            autorange: 'reversed' as const,
          },
          showlegend: true,
          legend: { orientation: 'h', y: -0.08, x: 0, font: { size: 10 } },
        }}
        config={baseConfig} style={{ width: '100%', height: `${dynamicHeight}px` }} useResizeHandler />
    </ChartCard>
  )
}
