/**
 * charts/individual/PatientCharts.tsx
 * PatientRadar, FacetedComposition
 */

import { useMemo, useState } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import { baseConfig, groupColor, mean, usePlotlyDownload, withAlpha, THEME, taxonColor } from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

const pid   = (r: SampleRow) => String(r.patient ?? r.sample_id).replace(/_T\d+$/, '')
const baseG = (r: SampleRow) => String(r.base_group ?? r.group).replace(/_T\d+$/, '')

// ── Patient Radar ──────────────────────────────────────────────────────────

interface RadarProps { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function PatientRadar({ rows, taxa, onExpand }: RadarProps) {
  const patients = [...new Set(rows.map(pid))].sort()
  const [selPat, setSelPat] = useState(patients[0] ?? '')
  const { plotRef, download } = usePlotlyDownload('patient-radar')
  const allBase = [...new Set(rows.map(baseG))].sort()

  const { traces, tableRows } = useMemo(() => {
    const r0  = rows.find(r => pid(r) === selPat && Number(r.time) === 0)
    const r84 = rows.filter(r => pid(r) === selPat && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
    const bg  = r0 ? baseG(r0) : allBase[0]
    const c   = groupColor(bg, allBase)
    const short = taxa.map(t => t.replace('aceae', '').slice(0, 10))
    const normRow = (r: SampleRow) => {
      const tot = taxa.reduce((s, t) => s + Number(r[t] ?? 0), 0) || 1
      return taxa.map(t => +(Number(r[t] ?? 0) / tot * 100).toFixed(1))
    }
    const gd = rows.filter(r => baseG(r) === bg && Number(r.time) === 0)
    const gmean = taxa.map(t => mean(gd.map(r => Number(r[t] ?? 0) / (taxa.reduce((s, tt) => s + Number(r[tt] ?? 0), 0) || 1) * 100)))
    const v0  = r0  ? normRow(r0)  : null
    const v84 = r84 ? normRow(r84) : null

    const traces: object[] = []
    if (v0) traces.push({
      type: 'scatterpolar', fill: 'toself',
      r: [...v0, v0[0]], theta: [...short, short[0]],
      fillcolor: withAlpha(c, 0.2), line: { color: withAlpha(c, 0.8), width: 2 },
      name: `${selPat} T0`, hovertemplate: '<b>%{theta}</b>: %{r:.1f}%<extra></extra>',
    })
    if (v84) traces.push({
      type: 'scatterpolar', fill: 'none',
      r: [...v84, v84[0]], theta: [...short, short[0]],
      line: { color: c, width: 2.5, dash: 'dash' },
      name: `${selPat} T84`, hovertemplate: '<b>%{theta}</b>: %{r:.1f}%<extra></extra>',
    })
    traces.push({
      type: 'scatterpolar', fill: 'none',
      r: [...gmean, gmean[0]], theta: [...short, short[0]],
      line: { color: '#C4A0A8', width: 1.5, dash: 'dot' },
      name: `${bg} mean T0`, hovertemplate: '<b>%{theta}</b>: %{r:.1f}%<extra></extra>',
    })

    const tableRows = taxa.map((t, i) => ({
      family: t.replace('aceae', ''),
      v0: v0 ? v0[i].toFixed(1) : '—',
      v84: v84 ? v84[i].toFixed(1) : '—',
      delta: v0 && v84 ? +(v84[i] - v0[i]).toFixed(1) : null,
    }))

    return { traces, tableRows }
  }, [selPat, rows, taxa])

  const controls = (
    <select value={selPat} onChange={e => setSelPat(e.target.value)}
      style={{ fontSize: '11px', padding: '3px 6px', borderRadius: '6px', border: '1px solid #E8C0AD', background: '#FFF8F4', fontFamily: THEME.font }}>
      {patients.map(p => <option key={p} value={p}>{p}</option>)}
    </select>
  )

  return (
    <ChartCard title="Patient Radar Profile"
      subtitle="T0 filled · T84 dashed · group mean dotted · hover for values"
      onExpand={onExpand} onDownload={download} controls={controls}>
      <div style={{ display: 'flex', gap: '24px', alignItems: 'flex-start' }}>
        <div style={{ flex: '0 0 420px' }}>
          <Plot ref={plotRef} data={traces as never[]}
            layout={{
              paper_bgcolor: THEME.paper, font: { family: THEME.font, color: THEME.text, size: 10 },
              height: 420, margin: { l: 60, r: 60, t: 30, b: 60 },
              polar: { radialaxis: { visible: true, range: [0, 45], ticksuffix: '%', tickfont: { size: 9 } },
                       gridshape: 'circular' },
              showlegend: true,
              legend: { orientation: 'h', y: -0.12, x: 0, font: { size: 10 } },
            }}
            config={baseConfig} style={{ width: '420px', height: '420px' }} useResizeHandler={false} />
        </div>
        <div style={{ flex: 1, overflowX: 'auto', fontSize: '12px' }}>
          <div style={{ fontWeight: 700, color: THEME.text, marginBottom: '6px', paddingBottom: '4px', borderBottom: '1px solid #F0D5C8' }}>
            {selPat}
          </div>
          <table style={{ borderCollapse: 'collapse', width: '100%' }}>
            <thead>
              <tr style={{ borderBottom: '1px solid #F0D5C8', fontSize: '9px', color: THEME.text3 }}>
                <th style={{ textAlign: 'left', padding: '2px 6px' }}>Family</th>
                <th style={{ textAlign: 'right', padding: '2px 6px' }}>T0 %</th>
                <th style={{ textAlign: 'right', padding: '2px 6px' }}>T84 %</th>
                <th style={{ textAlign: 'right', padding: '2px 6px' }}>Δ</th>
              </tr>
            </thead>
            <tbody>
              {tableRows.map(r => (
                <tr key={r.family} style={{ borderBottom: '1px solid #F8EDE8',
                  background: Math.abs(r.delta ?? 0) > 2 ? (r.delta! > 0 ? 'rgba(217,122,58,0.06)' : 'rgba(74,126,212,0.06)') : 'transparent' }}>
                  <td style={{ padding: '3px 6px', fontWeight: 600, color: THEME.text }}>{r.family}</td>
                  <td style={{ padding: '3px 6px', textAlign: 'right', color: THEME.text2 }}>{r.v0}</td>
                  <td style={{ padding: '3px 6px', textAlign: 'right', color: THEME.text2 }}>{r.v84}</td>
                  <td style={{ padding: '3px 6px', textAlign: 'right', fontWeight: 700,
                    color: r.delta === null ? THEME.text3 : r.delta > 2 ? '#D97A3A' : r.delta < -2 ? '#4A7ED4' : THEME.text3 }}>
                    {r.delta !== null ? (r.delta > 0 ? '+' : '') + r.delta : '—'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          <div style={{ fontSize: '8px', color: THEME.text3, marginTop: '6px' }}>
            Orange Δ = increase &gt;2% · Blue Δ = decrease &gt;2%
          </div>
        </div>
      </div>
    </ChartCard>
  )
}

// ── Faceted Small Multiples ────────────────────────────────────────────────

interface FacetProps { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function FacetedComposition({ rows, taxa, onExpand }: FacetProps) {
  const patients = [...new Set(rows.map(pid))].sort()
  const [selPat, setSelPat] = useState<string>('all')
  const { plotRef, download } = usePlotlyDownload('faceted-composition')

  const displayPatients = selPat === 'all' ? patients : [selPat]

  const traces = useMemo(() => {
    const out: object[] = []
    const shown = new Set<string>()
    displayPatients.forEach((p, pi) => {
      const r0  = rows.find(r => pid(r) === p && Number(r.time) === 0)
      const r84 = rows.filter(r => pid(r) === p && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
      ;[r0, r84].forEach((r, ti) => {
        if (!r) return
        const tp = ti === 0 ? 'T0' : 'T84'
        const tot = taxa.reduce((s, t) => s + Number(r[t] ?? 0), 0) || 1
        taxa.forEach(t => {
          out.push({
            type: 'bar', name: t,
            x: [`${p}<br>${tp}`],
            y: [+(Number(r[t] ?? 0) / tot * 100).toFixed(1)],
            marker: { color: taxonColor(t) },
            legendgroup: t, showlegend: !shown.has(t),
            hovertemplate: `<b>${p} ${tp}</b><br>${t}: %{y:.1f}%<extra></extra>`,
            offsetgroup: String(pi * 2 + ti),
          })
          shown.add(t)
        })
      })
    })
    return out
  }, [displayPatients, rows, taxa])

  const controls = (
    <select value={selPat} onChange={e => setSelPat(e.target.value)}
      style={{ fontSize: '11px', padding: '3px 6px', borderRadius: '6px', border: '1px solid #E8C0AD', background: '#FFF8F4', fontFamily: THEME.font }}>
      <option value="all">All patients</option>
      {patients.map(p => <option key={p} value={p}>{p}</option>)}
    </select>
  )

  const h = selPat === 'all' ? 380 : 340

  return (
    <ChartCard title="Individual Composition — Small Multiples"
      subtitle="T0 vs T84 per patient · each colour = one family"
      onExpand={onExpand} onDownload={download} controls={controls}>
      <Plot ref={plotRef} data={traces as never[]}
        layout={{
          autosize: true,
          barmode: 'stack', paper_bgcolor: THEME.paper, plot_bgcolor: THEME.bg,
          font: { family: THEME.font, color: THEME.text, size: 9 },
          height: h, margin: { l: 40, r: 16, t: 16, b: 80 },
          xaxis: { tickangle: -30, tickfont: { size: 8 }, gridcolor: 'transparent' },
          yaxis: { range: [0, 100], ticksuffix: '%', gridcolor: THEME.grid },
          showlegend: true, legend: { orientation: 'h', y: -0.35, x: 0, font: { size: 9 } },
        }}
        config={baseConfig} style={{ width: '100%', height: `${h}px` }} useResizeHandler />
    </ChartCard>
  )
}
