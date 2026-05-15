/**
 * charts/clinical/ClinicalCharts.tsx
 * ClinicalSlopegraph, ClinicalCorrelation
 */

import { useMemo, useState } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import {
  baseLayout, baseConfig, groupColor, mean,
  usePlotlyDownload, withAlpha, THEME, twoTailP,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

const pid   = (r: SampleRow) => String(r.patient ?? r.sample_id).replace(/_T\d+$/, '')
const baseG = (r: SampleRow) => String(r.base_group ?? r.group).replace(/_T\d+$/, '')

// ── Clinical Slopegraph ────────────────────────────────────────────────────

interface SlopeProps {
  rows: SampleRow[]
  field: 'sixmwt' | 'il18'
  label: string
  unit: string
  onExpand?: () => void
}

export function ClinicalSlopegraph({ rows, field, label, unit, onExpand }: SlopeProps) {
  const { plotRef, download } = usePlotlyDownload(`clinical-${field}`)
  const allBase = [...new Set(rows.map(baseG))].sort()
  const patients = [...new Set(rows.map(pid))]
  const hasData = rows.some(r => Number(r[field]) > 0)
  const shown = new Set<string>()

  const traces = useMemo(() => {
    if (!hasData) return []
    const out: object[] = []
    patients.forEach(p => {
      const r0  = rows.find(r => pid(r) === p && Number(r.time) === 0)
      const r84 = rows.filter(r => pid(r) === p && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
      if (!r0 || !r84 || !Number(r0[field]) || !Number(r84[field])) return
      const bg = baseG(r0); const c = groupColor(bg, allBase)
      const v0 = Number(r0[field]), v84 = Number(r84[field])
      const delta = v84 - v0
      out.push({
        type: 'scatter', mode: 'lines+markers',
        x: ['T0', 'T84'], y: [v0, v84],
        line: { color: withAlpha(c, 0.55), width: 2 },
        marker: { color: c, size: 8 },
        name: bg, legendgroup: bg, showlegend: !shown.has(bg),
        hovertemplate: `<b>${p}</b>: %{y:.0f} ${unit} (${delta > 0 ? '+' : ''}${delta.toFixed(0)})<extra></extra>`,
      })
      shown.add(bg)
    })
    // Group means
    allBase.forEach(bg => {
      const c = groupColor(bg, allBase)
      const m0  = mean(rows.filter(r => baseG(r) === bg && Number(r.time) === 0).map(r => Number(r[field])).filter(v => v > 0))
      const m84 = mean(rows.filter(r => baseG(r) === bg && Number(r.time) > 0).map(r => Number(r[field])).filter(v => v > 0))
      if (!m0 || !m84) return
      out.push({
        type: 'scatter', mode: 'text+lines+markers',
        x: ['T0', 'T84'], y: [m0, m84],
        text: [m0.toFixed(0), m84.toFixed(0)],
        textposition: ['top left', 'top right'],
        textfont: { size: 11, color: c, family: THEME.font },
        line: { color: c, width: 3.5, dash: 'dot' },
        marker: { color: c, size: 11, symbol: 'diamond' },
        name: `${bg} mean`,
        hovertemplate: `<b>${bg} mean</b>: %{y:.0f} ${unit}<extra></extra>`,
      })
    })
    return out
  }, [rows, field, hasData])

  return (
    <ChartCard title={`${label} — Individual Trajectories`}
      subtitle={`${unit} per patient · T0 → T84 · dashed = group mean`}
      onExpand={onExpand} onDownload={download}>
      {!hasData ? (
        <div style={{ height: '200px', display: 'flex', alignItems: 'center', justifyContent: 'center',
          color: THEME.text3, fontSize: '12px', border: '1.5px dashed #F0D5C8', borderRadius: '8px', fontStyle: 'italic' }}>
          No {field} data. Add column to metadata.tsv.
        </div>
      ) : (
        <Plot ref={plotRef} data={traces as never[]}
          layout={baseLayout({ height: 400, showlegend: true, yaxis: { title: { text: `${label} (${unit})`, font: { size: 11 } } } })}
          config={baseConfig} style={{ width: '100%', height: '400px' }} useResizeHandler />
      )}
    </ChartCard>
  )
}

// ── Clinical Correlation ───────────────────────────────────────────────────

interface CorrProps {
  rows: SampleRow[]
  yField: 'sixmwt' | 'il18'
  yLabel: string
  yUnit: string
  onExpand?: () => void
}

export function ClinicalCorrelation({ rows, yField, yLabel, yUnit, onExpand }: CorrProps) {
  const [metric, setMetric] = useState<'shannon' | 'simpson'>('shannon')
  const { plotRef, download } = usePlotlyDownload(`corr-${metric}-${yField}`)
  const allBase = [...new Set(rows.map(baseG))].sort()
  const xLabel  = metric === 'shannon' ? "Shannon H′" : 'Simpson'
  const valid   = rows.filter(r => Number(r[yField]) > 0 && Number(r[metric]) > 0)

  // Pearson r
  const { r, p } = useMemo(() => {
    if (valid.length < 3) return { r: 0, p: 1 }
    const xs = valid.map(r => Number(r[metric])), ys = valid.map(r => Number(r[yField]))
    const mx = mean(xs), my = mean(ys)
    const num = xs.reduce((s, x, i) => s + (x - mx) * (ys[i] - my), 0)
    const dx  = Math.sqrt(xs.reduce((s, x) => s + (x - mx) ** 2, 0))
    const dy  = Math.sqrt(ys.reduce((s, y) => s + (y - my) ** 2, 0))
    const rv  = dx && dy ? num / (dx * dy) : 0
    const n   = valid.length
    const t   = rv * Math.sqrt((n - 2) / (1 - rv * rv + 1e-10))
    const pv  = n > 2 ? twoTailP(t, n - 2) : 1
    return { r: +rv.toFixed(3), p: +pv.toFixed(3) }
  }, [valid, metric, yField])

  const xs = valid.map(r => Number(r[metric])), ys = valid.map(r => Number(r[yField]))
  const xMin = Math.min(...xs), xMax = Math.max(...xs)
  const slope = xs.length > 1 ? (ys.reduce((s, y, i) => s + (xs[i] - mean(xs)) * (y - mean(ys)), 0) /
    Math.max(xs.reduce((s, x) => s + (x - mean(xs)) ** 2, 0), 1e-10)) : 0
  const intercept = mean(ys) - slope * mean(xs)

  const traces = [
    ...allBase.map(bg => {
      const pts = valid.filter(r => baseG(r) === bg)
      const c   = groupColor(bg, allBase)
      return {
        type: 'scatter', mode: 'markers', name: bg,
        x: pts.map(r => Number(r[metric])), y: pts.map(r => Number(r[yField])),
        text: pts.map(r => String(r.sample_id)),
        marker: { color: withAlpha(c, 0.85), size: 9, line: { width: 1, color: 'white' } },
        hovertemplate: '<b>%{text}</b><br>%{x:.3f} / %{y:.0f}<extra></extra>',
      }
    }),
    {
      type: 'scatter', mode: 'lines', name: 'Regression',
      x: [xMin, xMax], y: [slope * xMin + intercept, slope * xMax + intercept],
      line: { color: withAlpha('#C4A0A8', 0.7), width: 1.5, dash: 'dash' },
      showlegend: false, hoverinfo: 'skip',
    },
  ]

  const sig = p < 0.05

  return (
    <ChartCard title={`${xLabel} vs ${yLabel}`}
      subtitle={`r = ${r}, p = ${p} · ${sig ? '✱ significant' : 'not significant'}`}
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
      {valid.length < 3 ? (
        <div style={{ height: '200px', display: 'flex', alignItems: 'center', justifyContent: 'center',
          color: THEME.text3, fontSize: '12px', border: '1.5px dashed #F0D5C8', borderRadius: '8px', fontStyle: 'italic' }}>
          No {yField} data available.
        </div>
      ) : (
        <Plot ref={plotRef} data={traces as never[]}
          layout={baseLayout({
            height: 380, showlegend: true,
            xaxis: { ...baseLayout().xaxis, title: { text: xLabel, font: { size: 11 } } },
            yaxis: { ...baseLayout().yaxis, title: { text: `${yLabel} (${yUnit})`, font: { size: 11 } } },
            annotations: [{
              x: 0.97, y: 0.03, xref: 'paper', yref: 'paper',
              text: `r = ${r}<br>p = ${p}<br>${sig ? '✱ significant' : 'not significant'}`,
              showarrow: false, align: 'right',
              bgcolor: 'rgba(255,248,244,0.85)', bordercolor: '#C4A0A8',
              font: { size: 10, color: THEME.text },
            }],
          })}
          config={baseConfig} style={{ width: '100%', height: '380px' }} useResizeHandler />
      )}
    </ChartCard>
  )
}
