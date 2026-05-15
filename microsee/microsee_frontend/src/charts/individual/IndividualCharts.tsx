/**
 * charts/individual/IndividualCharts.tsx
 * PairedSlopegraph, StabilityBar, DiversityRank — all in one file.
 */

import { useMemo, useState } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import {
  baseLayout, baseConfig, groupColor, mean,
  usePlotlyDownload, withAlpha, THEME,
} from '@/charts/shared/usePlotly'
import { brayCurtis } from '@/charts/beta/distances'
import type { SampleRow } from '@/types/sample'

const pid = (r: SampleRow) => String(r.patient ?? r.sample_id).replace(/_T\d+$/, '')
const baseG = (r: SampleRow) => String(r.base_group ?? r.group).replace(/_T\d+$/, '')

// ── Paired Slopegraph ──────────────────────────────────────────────────────

interface SlopeProps { rows: SampleRow[]; onExpand?: () => void }

export function PairedSlopegraph({ rows, onExpand }: SlopeProps) {
  const [metric, setMetric] = useState<'shannon' | 'simpson'>('shannon')
  const { plotRef, download } = usePlotlyDownload(`slopegraph-${metric}`)
  const allBase = [...new Set(rows.map(baseG))].sort()
  const patients = [...new Set(rows.map(pid))]
  const label = metric === 'shannon' ? "Shannon H′" : 'Simpson'
  const shown = new Set<string>()

  const traces = useMemo(() => {
    const out: object[] = []
    patients.forEach(p => {
      const r0  = rows.find(r => pid(r) === p && Number(r.time) === 0)
      const r84 = rows.filter(r => pid(r) === p && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
      if (!r0 || !r84) return
      const bg = baseG(r0); const c = groupColor(bg, allBase)
      out.push({
        type: 'scatter', mode: 'lines+markers',
        x: ['T0', 'T84'], y: [Number(r0[metric]), Number(r84[metric])],
        line: { color: withAlpha(c, 0.55), width: 1.8 },
        marker: { color: c, size: 7 },
        name: bg, legendgroup: bg, showlegend: !shown.has(bg),
        hovertemplate: `<b>${p}</b>: %{y:.3f}<extra></extra>`,
      })
      shown.add(bg)
    })
    // Group mean lines
    allBase.forEach(bg => {
      const c = groupColor(bg, allBase)
      const m0  = mean(rows.filter(r => baseG(r) === bg && Number(r.time) === 0).map(r => Number(r[metric])))
      const m84 = mean(rows.filter(r => baseG(r) === bg && Number(r.time) > 0).map(r => Number(r[metric])))
      out.push({
        type: 'scatter', mode: 'text+lines+markers',
        x: ['T0', 'T84'], y: [m0, m84],
        text: [m0.toFixed(2), m84.toFixed(2)],
        textposition: ['top left', 'top right'],
        textfont: { size: 11, color: c, family: THEME.font },
        line: { color: c, width: 3.5, dash: 'dot' },
        marker: { color: c, size: 11, symbol: 'diamond' },
        name: `${bg} mean`, showlegend: true,
        hovertemplate: `<b>${bg} mean</b>: %{y:.3f}<extra></extra>`,
      })
    })
    return out
  }, [rows, metric])

  return (
    <ChartCard title={`Paired Slopegraph — ${label}`}
      subtitle="Each line = one patient · dashed = group mean · flat = no diversity change"
      onExpand={onExpand} onDownload={download}
      controls={
        <select
          value={metric}
          onChange={(e) => setMetric(e.target.value as 'shannon' | 'simpson')}
          style={{ fontSize: '12px', padding: '2px 4px' }}
        >
          <option value="shannon">Shannon H′</option>
          <option value="simpson">Simpson</option>
        </select>
      }>
      <Plot ref={plotRef} data={traces as never[]}
        layout={baseLayout({ height: 400, showlegend: true, yaxis: { title: { text: label, font: { size: 11 } } } })}
        config={baseConfig} style={{ width: '100%', height: '400px' }} useResizeHandler />
    </ChartCard>
  )
}

// ── Stability Bar ──────────────────────────────────────────────────────────

interface StabProps { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function StabilityBar({ rows, taxa, onExpand }: StabProps) {
  const { plotRef, download } = usePlotlyDownload('stability')
  const allBase = [...new Set(rows.map(baseG))].sort()

  const scores = useMemo(() => {
    return [...new Set(rows.map(pid))].map(p => {
      const r0  = rows.find(r => pid(r) === p && Number(r.time) === 0)
      const r84 = rows.filter(r => pid(r) === p && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
      if (!r0 || !r84) return null
      return { pid: p, bc: +brayCurtis(r0, r84, taxa).toFixed(3), group: baseG(r0) }
    }).filter(Boolean).sort((a, b) => a!.bc - b!.bc) as { pid: string; bc: number; group: string }[]
  }, [rows, taxa])

  return (
    <ChartCard title="Microbiome Stability Score (T0 vs T84)"
      subtitle="Bray-Curtis dissimilarity · 0 = identical · 1 = completely different"
      onExpand={onExpand} onDownload={download}>
      <Plot ref={plotRef}
        data={[{
          type: 'bar', orientation: 'h',
          x: scores.map(s => s.bc), y: scores.map(s => s.pid),
          marker: { color: scores.map(s => withAlpha(groupColor(s.group, allBase), 0.8)) },
          text: scores.map(s => String(s.bc)), textposition: 'outside',
          hovertemplate: '<b>%{y}</b> BC = %{x:.3f}<extra></extra>',
        }]}
        layout={{
          paper_bgcolor: THEME.paper, plot_bgcolor: THEME.bg,
          font: { family: THEME.font, color: THEME.text, size: 11 },
          height: 400, margin: { l: 80, r: 60, t: 28, b: 40 },
          xaxis: { range: [0, Math.max(...scores.map(s => s.bc)) * 1.25 || 0.5],
                   title: { text: 'Bray-Curtis dissimilarity', font: { size: 11 } }, gridcolor: THEME.grid },
          yaxis: { gridcolor: 'transparent' }, showlegend: false,
        }}
        config={baseConfig} style={{ width: '100%', height: '400px' }} useResizeHandler />
    </ChartCard>
  )
}

// ── Diversity Rank Plot ────────────────────────────────────────────────────

interface RankProps { rows: SampleRow[]; onExpand?: () => void }

export function DiversityRank({ rows, onExpand }: RankProps) {
  const [metric, setMetric] = useState<'shannon' | 'simpson'>('shannon')
  const { plotRef, download } = usePlotlyDownload(`diversity-rank-${metric}`)
  const allBase = [...new Set(rows.map(baseG))].sort()
  const label   = metric === 'shannon' ? "Shannon H′" : 'Simpson'
  const sorted  = [...rows].sort((a, b) => Number(a[metric]) - Number(b[metric]))

  const traces = allBase.map(bg => {
    const pts = sorted.map((r, i) => ({ i, r })).filter(({ r }) => baseG(r) === bg)
    const c   = groupColor(bg, allBase)
    const sym = (r: SampleRow) => Number(r.time) === 0 ? 'circle' : 'diamond'
    return {
      type: 'scatter' as const, mode: 'markers' as const, name: bg,
      x: pts.map(({ i }) => i),
      y: pts.map(({ r }) => Number(r[metric])),
      marker: { color: pts.map(({ r }) => Number(r.time) === 0 ? withAlpha(c, 0.6) : c),
                size: pts.map(({ r }) => Number(r.time) === 0 ? 7 : 9),
                symbol: pts.map(({ r }) => sym(r)), line: { width: 1, color: 'white' } },
      text: pts.map(({ r }) => String(r.sample_id)),
      hovertemplate: '<b>%{text}</b>: %{y:.3f}<extra></extra>',
    }
  })

  return (
    <ChartCard title={`Diversity Rank Plot — ${label}`}
      subtitle="Samples ranked low→high · circles = T0 · diamonds = T84"
      onExpand={onExpand} onDownload={download}
      controls={
        <select
          value={metric}
          onChange={(e) => setMetric(e.target.value as 'shannon' | 'simpson')}
          style={{ fontSize: '12px', padding: '2px 4px' }}
        >
          <option value="shannon">Shannon H′</option>
          <option value="simpson">Simpson</option>
        </select>
      }>
      <Plot ref={plotRef} data={traces}
        layout={baseLayout({
          height: 360, showlegend: true,
          xaxis: { ...baseLayout().xaxis, showticklabels: false, title: { text: 'Rank', font: { size: 11 } }, gridcolor: 'transparent' },
          yaxis: { ...baseLayout().yaxis, title: { text: label, font: { size: 11 } } },
        })}
        config={baseConfig} style={{ width: '100%', height: '360px' }} useResizeHandler />
    </ChartCard>
  )
}
