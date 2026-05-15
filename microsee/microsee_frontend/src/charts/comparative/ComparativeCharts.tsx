/**
 * charts/comparative/ComparativeCharts.tsx
 * DifferentialAbundance, VolcanoPlot, AbundanceHeatmap, CorrelationMatrix
 */

import { useMemo, useState } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import {
  baseConfig, groupColor, mean, uniqueGroups,
  usePlotlyDownload, withAlpha, THEME, twoTailP,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

// ── Differential Abundance ─────────────────────────────────────────────────

interface DiffProps { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function DifferentialAbundance({ rows, taxa, onExpand }: DiffProps) {
  const { plotRef, download } = usePlotlyDownload('differential-abundance')
  const groups = uniqueGroups(rows)
  const [grpA, setGrpA] = useState(groups[0] ?? '')
  const [grpB, setGrpB] = useState(groups[1] ?? '')

  const diffs = useMemo(() => {
    return taxa.map(t => {
      const a = mean(rows.filter(r => r.group === grpA).map(r => Number(r[t] ?? 0)))
      const b = mean(rows.filter(r => r.group === grpB).map(r => Number(r[t] ?? 0)))
      const lfc = Math.log2((a + 1e-6) / (b + 1e-6))
      return { taxon: t, lfc: +lfc.toFixed(3) }
    }).sort((a, b) => b.lfc - a.lfc)
  }, [rows, taxa, grpA, grpB])

  const colors = diffs.map(d => d.lfc > 0
    ? withAlpha(groupColor(grpA, groups), 0.8)
    : withAlpha(groupColor(grpB, groups), 0.8))

  const controls = (
    <div style={{ display: 'flex', gap: '4px', alignItems: 'center', fontSize: '10px' }}>
      <select value={grpA} onChange={e => setGrpA(e.target.value)}
        style={{ fontSize: '10px', padding: '2px 5px', borderRadius: '5px', border: '1px solid #E8C0AD', background: '#FFF8F4' }}>
        {groups.map(g => <option key={g}>{g}</option>)}
      </select>
      <span style={{ color: THEME.text3 }}>vs</span>
      <select value={grpB} onChange={e => setGrpB(e.target.value)}
        style={{ fontSize: '10px', padding: '2px 5px', borderRadius: '5px', border: '1px solid #E8C0AD', background: '#FFF8F4' }}>
        {groups.map(g => <option key={g}>{g}</option>)}
      </select>
    </div>
  )

  return (
    <ChartCard title="Differential Abundance" subtitle="Log₂ fold change · positive = higher in Group A"
      onExpand={onExpand} onDownload={download} controls={controls}>
      <Plot ref={plotRef}
        data={[{
          type: 'bar', orientation: 'h',
          x: diffs.map(d => d.lfc), y: diffs.map(d => d.taxon.replace('aceae', '')),
          marker: { color: colors }, text: diffs.map(d => d.lfc > 0 ? `+${d.lfc}` : String(d.lfc)),
          textposition: 'outside', hovertemplate: '<b>%{y}</b>: log₂FC = %{x:.3f}<extra></extra>',
        }]}
        layout={{
          paper_bgcolor: THEME.paper, plot_bgcolor: THEME.bg,
          font: { family: THEME.font, color: THEME.text, size: 10 },
          height: 360, margin: { l: 110, r: 60, t: 28, b: 40 },
          xaxis: { title: { text: 'Log₂ Fold Change', font: { size: 11 } }, zeroline: true, zerolinecolor: THEME.text3, gridcolor: THEME.grid },
          yaxis: { autorange: 'reversed', gridcolor: 'transparent' },
          showlegend: false,
        }}
        config={baseConfig} style={{ width: '100%', height: '360px' }} useResizeHandler />
    </ChartCard>
  )
}

// ── Volcano Plot ───────────────────────────────────────────────────────────


export function VolcanoPlot({ rows, taxa, onExpand }: DiffProps) {
  const { plotRef, download } = usePlotlyDownload('volcano')
  const groups = uniqueGroups(rows)
  const [grpA, setGrpA] = useState(groups[0] ?? '')
  const [grpB, setGrpB] = useState(groups[1] ?? '')

  const pts = useMemo(() => taxa.map(taxon => {
    const aVals = rows.filter(r => r.group === grpA).map(r => Number(r[taxon] ?? 0))
    const bVals = rows.filter(r => r.group === grpB).map(r => Number(r[taxon] ?? 0))
    const a = mean(aVals), b = mean(bVals)
    // Pseudo-count prevents log(0); LFC positive = higher in A
    const lfc = Math.log2((a + 1e-6) / (b + 1e-6))
    const na = aVals.length, nb = bVals.length
    const va = aVals.reduce((s, v) => s + (v - a) ** 2, 0) / Math.max(na - 1, 1)
    const vb = bVals.reduce((s, v) => s + (v - b) ** 2, 0) / Math.max(nb - 1, 1)
    const sa = va / na, sb = vb / nb
    const se = Math.sqrt(sa + sb)
    const t  = se > 0 ? (a - b) / se : 0
    // Welch-Satterthwaite df — consistent with Welch SE above
    const df = Math.max(1, Math.round((sa + sb) ** 2 / (sa ** 2 / Math.max(na - 1, 1) + sb ** 2 / Math.max(nb - 1, 1))))
    const p  = twoTailP(t, df)
    return { taxon, lfc: +lfc.toFixed(3), negLogP: +(-Math.log10(p)).toFixed(2), p }
  }), [rows, taxa, grpA, grpB])

  const sig    = pts.filter(p => p.p < 0.05 && Math.abs(p.lfc) > 1)
  const colors = pts.map(p =>
    p.p < 0.05 && Math.abs(p.lfc) > 1
      ? (p.lfc > 0 ? withAlpha(groupColor(grpA, groups), 0.9) : withAlpha(groupColor(grpB, groups), 0.9))
      : withAlpha('#C4A0A8', 0.45))

  const controls = (
    <div style={{ display: 'flex', gap: '4px', alignItems: 'center', fontSize: '10px' }}>
      <select value={grpA} onChange={e => setGrpA(e.target.value)}
        style={{ fontSize: '10px', padding: '2px 5px', borderRadius: '5px', border: '1px solid #E8C0AD', background: '#FFF8F4' }}>
        {groups.map(g => <option key={g}>{g}</option>)}
      </select>
      <span style={{ color: THEME.text3 }}>vs</span>
      <select value={grpB} onChange={e => setGrpB(e.target.value)}
        style={{ fontSize: '10px', padding: '2px 5px', borderRadius: '5px', border: '1px solid #E8C0AD', background: '#FFF8F4' }}>
        {groups.map(g => <option key={g}>{g}</option>)}
      </select>
    </div>
  )

  const pThresh = +(-Math.log10(0.05)).toFixed(4)

  return (
    <ChartCard title="Volcano Plot" subtitle="Coloured = p < 0.05 & |LFC| > 1 · right = higher in A · left = higher in B"
      onExpand={onExpand} onDownload={download} controls={controls}>
      <Plot ref={plotRef}
        data={[{
          type:          'scatter' as const,
          mode:          'text+markers' as const,
          x:             pts.map(p => p.lfc),
          y:             pts.map(p => p.negLogP),
          text:          pts.map(p => sig.includes(p) ? p.taxon.replace('aceae', '') : ''),
          textposition:  'top center' as const,
          textfont:      { size: 8 },
          customdata:    pts.map(p => p.taxon),
          marker:        { color: colors, size: 10, line: { width: 1, color: 'white' } },
          hovertemplate: '<b>%{customdata}</b><br>LFC = %{x:.2f}<br>-log₁₀(p) = %{y:.2f}<extra></extra>',
        }]}
        layout={{
          paper_bgcolor: THEME.paper, plot_bgcolor: THEME.bg,
          font: { family: THEME.font, color: THEME.text, size: 10 },
          height: 360, margin: { l: 52, r: 16, t: 28, b: 40 },
          xaxis: { title: { text: `Log₂ Fold Change  (${grpA} / ${grpB})`, font: { size: 11 } },
            zeroline: true, zerolinecolor: THEME.text3, gridcolor: THEME.grid },
          yaxis: { title: { text: '-log₁₀(p)', font: { size: 11 } }, gridcolor: THEME.grid },
          shapes: [
            { type: 'line', x0:  1, x1:  1, y0: 0, y1: 1, xref: 'x', yref: 'paper', line: { color: THEME.text3, dash: 'dot', width: 1 } },
            { type: 'line', x0: -1, x1: -1, y0: 0, y1: 1, xref: 'x', yref: 'paper', line: { color: THEME.text3, dash: 'dot', width: 1 } },
            { type: 'line', x0: 0,  x1: 1,  y0: pThresh, y1: pThresh, xref: 'paper', yref: 'y', line: { color: THEME.text3, dash: 'dot', width: 1 } },
          ],
          showlegend: false,
        }}
        config={baseConfig} style={{ width: '100%', height: '360px' }} useResizeHandler />
    </ChartCard>
  )
}

// ── Abundance Heatmap ──────────────────────────────────────────────────────

export function AbundanceHeatmap({ rows, taxa, onExpand }: DiffProps) {
  const { plotRef, download } = usePlotlyDownload('abundance-heatmap')
  const groups = uniqueGroups(rows)
  const sortedRows = [...rows].sort((a, b) => groups.indexOf(a.group) - groups.indexOf(b.group))
  const sortedTaxa = [...taxa].sort((a, b) => mean(rows.map(r => Number(r[b] ?? 0))) - mean(rows.map(r => Number(r[a] ?? 0))))
  const h = Math.max(420, Math.min(650, sortedTaxa.length * 34 + 150))

  return (
    <ChartCard title="Abundance Heatmap" subtitle="Samples × families · colour = relative abundance % · samples sorted by group"
      onExpand={onExpand} onDownload={download}>
      <Plot ref={plotRef}
        data={[{
          type: 'heatmap',
          z: sortedTaxa.map(t => sortedRows.map(r => Number(r[t] ?? 0))),
          x: sortedRows.map(r => String(r.sample_id)),
          y: sortedTaxa,
          colorscale: 'Oranges', colorbar: { title: { text: '%' }, len: 0.8, thickness: 14 },
          hovertemplate: '<b>%{x}</b><br>%{y}: %{z:.1f}%<extra></extra>',
          xgap: 1.5, ygap: 1.5,
        } as any]}
        layout={{
          autosize: true,
          paper_bgcolor: THEME.paper, plot_bgcolor: '#000000',
          font: { family: THEME.font, color: THEME.text, size: 9 },
          height: h, margin: { l: 140, r: 80, t: 28, b: 110 },
          xaxis: { tickangle: -40, tickfont: { size: 9 } },
          yaxis: { tickfont: { size: 10 } },
        }}
        config={baseConfig} style={{ width: '100%', height: `${h}px` }} useResizeHandler />
    </ChartCard>
  )
}

// ── Correlation Matrix ─────────────────────────────────────────────────────

export function CorrelationMatrix({ rows, taxa, onExpand }: DiffProps) {
  const { plotRef, download } = usePlotlyDownload('correlation-matrix')
  const short = taxa.map(t => t.replace('aceae', ''))

  const corr = useMemo(() => {
    const rank = (vals: number[]) => {
      const sorted = [...vals].map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v)
      const r = new Array(vals.length)
      let i = 0
      while (i < sorted.length) {
        let j = i
        while (j < sorted.length && sorted[j].v === sorted[i].v) j++
        const avg = (i + j + 1) / 2
        for (let k = i; k < j; k++) r[sorted[k].i] = avg
        i = j
      }
      return r as number[]
    }
    return taxa.map(ta => taxa.map(tb => {
      const xs = rank(rows.map(r => Number(r[ta] ?? 0)))
      const ys = rank(rows.map(r => Number(r[tb] ?? 0)))
      const mx = mean(xs), my = mean(ys)
      const num = xs.reduce((s, x, i) => s + (x - mx) * (ys[i] - my), 0)
      const dx  = Math.sqrt(xs.reduce((s, x) => s + (x - mx) ** 2, 0))
      const dy  = Math.sqrt(ys.reduce((s, y) => s + (y - my) ** 2, 0))
      return dx && dy ? +(num / (dx * dy)).toFixed(3) : 0
    }))
  }, [rows, taxa])

  const n = taxa.length
  const h = Math.max(460, Math.min(720, n * 46 + 160))

  return (
    <ChartCard title="Taxon Correlation Matrix" subtitle="Spearman ρ · +1 = co-occur · −1 = inversely related"
      onExpand={onExpand} onDownload={download}>
      <Plot ref={plotRef}
        data={[{
          type: 'heatmap', z: corr, x: short, y: short,
          colorscale: [
            [0.0, '#2050A0'], [0.35, '#A0C0F8'],
            [0.5, '#FFF8F4'],
            [0.65, '#F8C0A0'], [1.0, '#C03030'],
          ],
          zmid: 0, zmin: -1, zmax: 1,
          colorbar: { title: { text: 'r' }, len: 0.8, thickness: 14 },
          text: corr.map(row => row.map(v => v.toFixed(2))),
          hovertemplate: '<b>%{x} × %{y}</b>: r = %{z:.3f}<extra></extra>',
          xgap: 1.5, ygap: 1.5,
        } as any]}
        layout={{
          autosize: true,
          paper_bgcolor: THEME.paper, plot_bgcolor: '#000000',
          font: { family: THEME.font, color: THEME.text, size: 9 },
          height: h, margin: { l: 120, r: 80, t: 28, b: 120 },
          xaxis: { tickangle: -40, tickfont: { size: 10 } },
          yaxis: { tickfont: { size: 10 } },
        }}
        config={baseConfig} style={{ width: '100%', height: `${h}px` }} useResizeHandler />
    </ChartCard>
  )
}