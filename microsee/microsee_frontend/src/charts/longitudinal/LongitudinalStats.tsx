/**
 * charts/longitudinal/LongitudinalChart.tsx
 * charts/stats/StatsCharts.tsx — combined file
 */

import { useMemo, useState } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import {
  baseLayout, baseConfig, groupColor, mean,
  usePlotlyDownload, withAlpha, THEME,
} from '@/charts/shared/usePlotly'
import { buildDistMatrix, brayCurtis } from '@/charts/beta/distances'
import type { SampleRow } from '@/types/sample'

const pid   = (r: SampleRow) => String(r.patient ?? r.sample_id).replace(/_T\d+$/, '')
const baseG = (r: SampleRow) => String(r.base_group ?? r.group).replace(/_T\d+$/, '')

// ── Longitudinal Diversity ─────────────────────────────────────────────────

interface LongProps { rows: SampleRow[]; onExpand?: () => void }

export function LongitudinalChart({ rows, onExpand }: LongProps) {
  const [metric, setMetric] = useState<'shannon' | 'simpson'>('shannon')
  const { plotRef, download } = usePlotlyDownload(`longitudinal-${metric}`)
  const allBase = [...new Set(rows.map(baseG))].sort()
  const patients = [...new Set(rows.map(pid))]
  const label    = metric === 'shannon' ? "Shannon H′" : 'Simpson'
  const shown    = new Set<string>()

  const hasTime = rows.some(r => r.time !== null && r.time !== undefined)

  const traces = useMemo(() => {
    if (!hasTime) return []
    const out: object[] = []
    patients.forEach(p => {
      const pr = rows.filter(r => pid(r) === p).sort((a, b) => Number(a.time) - Number(b.time))
      if (pr.length < 2) return
      const bg = baseG(pr[0]); const c = groupColor(bg, allBase)
      out.push({
        type: 'scatter', mode: 'lines+markers',
        x: pr.map(r => Number(r.time)), y: pr.map(r => Number(r[metric])),
        line: { color: withAlpha(c, 0.45), width: 1.5 },
        marker: { color: c, size: 7 },
        name: bg, legendgroup: bg, showlegend: !shown.has(bg),
        text: pr.map(r => String(r.sample_id)),
        hovertemplate: '<b>%{text}</b> (day %{x}): %{y:.3f}<extra></extra>',
      })
      shown.add(bg)
    })
    // Group means
    allBase.forEach(bg => {
      const c = groupColor(bg, allBase)
      const times = [...new Set(rows.filter(r => baseG(r) === bg).map(r => Number(r.time)))].sort((a, b) => a - b)
      const means = times.map(t => mean(rows.filter(r => baseG(r) === bg && Number(r.time) === t).map(r => Number(r[metric]))))
      out.push({
        type: 'scatter', mode: 'text+lines+markers',
        x: times, y: means,
        text: means.map(m => m.toFixed(2)),
        textposition: 'top center',
        textfont: { size: 10, color: c, family: THEME.font },
        line: { color: c, width: 3.5 }, marker: { color: c, size: 10, symbol: 'diamond' },
        name: `${bg} mean`, showlegend: true,
        hovertemplate: `<b>${bg} mean</b> day %{x}: %{y:.3f}<extra></extra>`,
      })
    })
    return out
  }, [rows, metric, hasTime])

  return (
    <ChartCard title={`Longitudinal Diversity — ${label}`}
      subtitle="Individual patient trajectories · bold lines = group means"
      onExpand={onExpand} onDownload={download} fullWidth
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
      {!hasTime ? (
        <div style={{ height: '200px', display: 'flex', alignItems: 'center', justifyContent: 'center',
          color: THEME.text3, fontSize: '12px', border: '1.5px dashed #F0D5C8', borderRadius: '8px', fontStyle: 'italic' }}>
          No timepoint data found in metadata.
        </div>
      ) : (
        <Plot ref={plotRef} data={traces as never[]}
          layout={baseLayout({
            height: 420, showlegend: true,
            xaxis: { ...baseLayout().xaxis, title: { text: 'Time (days)', font: { size: 11 } } },
            yaxis: { ...baseLayout().yaxis, title: { text: label, font: { size: 11 } } },
          })}
          config={baseConfig} style={{ width: '100%', height: '420px' }} useResizeHandler />
      )}
    </ChartCard>
  )
}

// ── PERMANOVA Table ────────────────────────────────────────────────────────

interface PermProps { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function PERMANOVATable({ rows, taxa, onExpand }: PermProps) {
  const results = useMemo(() => {
    if (rows.length < 4 || !taxa.length) return []
    const mat   = buildDistMatrix(rows, taxa, brayCurtis)
    const flat  = mat.flat()
    const n     = rows.length
    const dataSeed = rows.reduce((s, r, i) => (s ^ (String(r.sample_id).charCodeAt(0) * (i + 1))) & 0x7fffffff, 0x12345678)
    const rng   = { v: dataSeed || 42 }
    const rand  = () => { rng.v = (rng.v * 1664525 + 1013904223) & 0xffffffff; return Math.abs(rng.v) / 0xffffffff }

    function pseudoF(labels: string[]) {
      const gs = [...new Set(labels)]
      let ssTot = 0
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) ssTot += flat[i * n + j] ** 2
      ssTot /= n
      let ssW = 0
      gs.forEach(g => {
        const idx = labels.map((l, i) => l === g ? i : -1).filter(i => i >= 0)
        const ng  = idx.length
        if (ng < 2) return
        let ss = 0
        for (let a = 0; a < idx.length; a++) for (let b = a + 1; b < idx.length; b++) ss += flat[idx[a] * n + idx[b]] ** 2
        ssW += ss / ng
      })
      const ssB = ssTot - ssW
      const r2  = ssB / (ssTot || 1)
      const k   = gs.length
      const F   = (ssB / (k - 1 || 1)) / (ssW / (n - k || 1))
      return { F, R2: r2 }
    }

    function permute(labels: string[], nperm = 99) {
      const obs  = pseudoF(labels)
      let cnt    = 1
      const lab  = [...labels]
      for (let i = 0; i < nperm; i++) {
        // Fisher-Yates shuffle
        for (let j = lab.length - 1; j > 0; j--) { const k = Math.floor(rand() * (j + 1)); [lab[j], lab[k]] = [lab[k], lab[j]] }
        if (pseudoF(lab).F >= obs.F) cnt++
      }
      return { ...obs, p: cnt / (nperm + 1) }
    }

    const tests = [
      { name: 'Supplementation Group', labels: rows.map(r => baseG(r)) },
      { name: 'Timepoint',             labels: rows.map(r => String(r.time ?? r.timepoint ?? '?')) },
      { name: 'Full Subgroup (4)',      labels: rows.map(r => String(r.group)) },
      { name: 'Individual',            labels: rows.map(r => pid(r)) },
    ]

    return tests.map(t => ({ name: t.name, ...permute(t.labels) }))
      .sort((a, b) => b.R2 - a.R2)
  }, [rows, taxa])

  const topFactor = results[0]

  return (
    <ChartCard title="PERMANOVA" subtitle="Bray-Curtis · 99 permutations · tests group / timepoint / individual"
      onExpand={onExpand}>
      <div style={{ fontSize: '10px', color: THEME.text3, marginBottom: '8px' }}>
        Bray–Curtis distances · 99 permutations
      </div>
      <table style={{ borderCollapse: 'collapse', width: '100%', fontSize: '11px' }}>
        <thead>
          <tr style={{ borderBottom: '1px solid #F0D5C8', background: '#FFF4EE' }}>
            {['Variable', 'R²', 'Pseudo-F', 'p-value', 'Result'].map(h => (
              <th key={h} style={{ padding: '6px 10px', textAlign: 'left', fontWeight: 700, color: THEME.text }}>{h}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {results.map(r => (
            <tr key={r.name} style={{ borderBottom: '1px solid #F8EDE8' }}>
              <td style={{ padding: '5px 10px', fontWeight: 600 }}>{r.name}</td>
              <td style={{ padding: '5px 10px' }}>{r.R2.toFixed(3)}</td>
              <td style={{ padding: '5px 10px' }}>{r.F.toFixed(2)}</td>
              <td style={{ padding: '5px 10px', fontWeight: 700, color: r.p < 0.05 ? '#2FA896' : THEME.text2 }}>{r.p.toFixed(3)}</td>
              <td style={{ padding: '5px 10px', fontSize: '10px', color: r.p < 0.05 ? '#2FA896' : THEME.text3 }}>
                {r.p < 0.05 ? '✓ Significant' : 'Not significant'}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
      {topFactor && (
        <div style={{ marginTop: '10px', padding: '9px 13px', background: '#FFE6DF', borderLeft: '3px solid #C05070',
          borderRadius: '0 8px 8px 0', fontSize: '11px', color: '#8B3A52', lineHeight: 1.5 }}>
          💡 <strong>{topFactor.name}</strong> explains the most variance (R² = {topFactor.R2.toFixed(3)}).
          {topFactor.name === 'Individual' && ' This confirms each person has a unique stable microbiome that dominates over any group-level effect.'}
        </div>
      )}
    </ChartCard>
  )
}

// ── Diversity Summary Table ────────────────────────────────────────────────

interface SummProps { rows: SampleRow[]; onExpand?: () => void }

export function DiversitySummary({ rows }: SummProps) {
  const groups = [...new Set(rows.map(r => String(r.group)))].sort()
  const sd = (arr: number[]) => {
    if (arr.length < 2) return 0
    const m = mean(arr)
    return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / (arr.length - 1))
  }
  const stats = groups.map(g => {
    const gd = rows.filter(r => r.group === g)
    const sh = gd.map(r => Number(r.shannon)), si = gd.map(r => Number(r.simpson))
    return {
      group: g, n: gd.length,
      shannon: `${mean(sh).toFixed(3)} ± ${sd(sh).toFixed(3)}`,
      simpson: `${mean(si).toFixed(3)} ± ${sd(si).toFixed(3)}`,
    }
  })

  return (
    <ChartCard title="Diversity Summary" subtitle="Mean ± SD per group">
      <table style={{ borderCollapse: 'collapse', width: '100%', fontSize: '11px' }}>
        <thead>
          <tr style={{ borderBottom: '1px solid #F0D5C8', background: '#FFF4EE' }}>
            {['Group', 'n', 'Shannon H′ (mean ± SD)', 'Simpson (mean ± SD)'].map(h => (
              <th key={h} style={{ padding: '6px 10px', textAlign: 'left', fontWeight: 700, color: THEME.text }}>{h}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {stats.map(s => (
            <tr key={s.group} style={{ borderBottom: '1px solid #F8EDE8' }}>
              <td style={{ padding: '5px 10px', fontWeight: 700, color: groupColor(s.group, groups) }}>{s.group}</td>
              <td style={{ padding: '5px 10px' }}>{s.n}</td>
              <td style={{ padding: '5px 10px' }}>{s.shannon}</td>
              <td style={{ padding: '5px 10px' }}>{s.simpson}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </ChartCard>
  )
}
