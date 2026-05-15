/**
 * charts/alpha/RarefactionCurves.tsx
 *
 * Rarefaction curves using the Hurlbert (1971) expected species formula.
 * E[S(n)] = sum_i { 1 - C(N - N_i, n) / C(N, n) }
 * approximated as: sum_i { 1 - (1 - p_i)^n }
 *
 * Each sample has virtual ASVs generated from family-level abundances
 * to produce realistic within-sample richness (60–280 virtual ASVs).
 *
 * Replaces buildRar() from MicroSee.html.
 */

import { useMemo }    from 'react'
import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseLayout, baseConfig, groupColor, mean,
  uniqueGroups, usePlotlyDownload, withAlpha,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  taxa:     string[]
  onExpand?: () => void
}

/** Deterministic hash for per-sample noise (replaces Math.random) */
function hashNoise(s: string): number {
  let h = 5381
  for (const c of s) h = ((h << 5) + h) + c.charCodeAt(0)
  return (Math.abs(h) % 1000) / 1000
}

function buildCurve(row: SampleRow, taxa: string[]): { x: number[]; y: number[] } {
  const props  = taxa.map((t) => Number(row[t] ?? 0))
  const tot    = props.reduce((a, b) => a + b, 0) || 1
  const sh     = Number(row.shannon) || 1.5
  const noise  = hashNoise(String(row.sample_id)) * 0.3 - 0.15

  // Virtual ASV richness: 60–280 based on Shannon diversity
  const richness = Math.max(40, Math.min(280, Math.round(sh * 65 * (1 + noise))))

  // Distribute virtual ASVs across families (log-normal within family)
  const asvProbs: number[] = []
  taxa.forEach((t, fi) => {
    const famP = props[fi] / tot
    const nAsv = Math.max(1, Math.round(famP * richness))
    const base = famP / nAsv
    for (let a = 0; a < nAsv; a++) {
      const lnv = Math.exp((hashNoise(t + a) - 0.5) * 1.2)
      asvProbs.push(base * lnv * (a === 0 ? 2.5 : 1))
    }
  })
  const pTot = asvProbs.reduce((a, b) => a + b, 0)
  const p    = asvProbs.map((x) => x / pTot)

  // Sample scale (40%–160% of max depth) for inter-sample variation
  const sampleScale = 0.4 + hashNoise(String(row.sample_id) + 'scale') * 1.2
  const N    = 50_000
  const steps = 30
  const depths = Array.from({ length: steps + 1 }, (_, i) =>
    Math.round(Math.pow(i / steps, 0.7) * N),
  )

  const y = depths.map((d) => {
    const sd = d * sampleScale
    return p.reduce((s, pi) => (pi > 0 ? s + 1 - Math.pow(1 - pi, sd) : s), 0)
  })

  return { x: depths, y }
}

export function RarefactionCurves({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('rarefaction')
  const groups = uniqueGroups(rows)

  const traces = useMemo(() => {
    const result = []

    // Individual sample curves (thin, semi-transparent)
    for (const row of rows) {
      const { x, y } = buildCurve(row, taxa)
      const c        = groupColor(row.group, groups)
      result.push({
        type:          'scatter'    as const,
        mode:          'lines'     as const,
        x, y,
        line:          { color: withAlpha(c, 0.5), width: 1.2 },
        name:          String(row.sample_id),
        showlegend:    false,
        hovertemplate: `<b>${row.sample_id}</b><br>Depth: %{x:,}<br>Species: %{y:.1f}<extra></extra>`,
      })
    }

    // Group mean lines (bold)
    groups.forEach((g) => {
      const gRows = rows.filter((r) => r.group === g)
      const curves = gRows.map((r) => buildCurve(r, taxa))
      const N      = 50_000
      const steps  = 30
      const depths = Array.from({ length: steps + 1 }, (_, i) =>
        Math.round(Math.pow(i / steps, 0.7) * N),
      )
      const meanY = depths.map((_, di) =>
        mean(curves.map((c) => c.y[di] ?? 0)),
      )
      const c = groupColor(g, groups)
      result.push({
        type:          'scatter'    as const,
        mode:          'lines'     as const,
        x:             depths,
        y:             meanY,
        line:          { color: c, width: 2.5 },
        name:          g,
        showlegend:    true,
        hovertemplate: `<b>${g} mean</b><br>Depth: %{x:,}<br>Expected: %{y:.1f}<extra></extra>`,
      })
    })

    return result
  }, [rows, taxa, groups])

  const layout = baseLayout({
    height:     380,
    xaxis:      { title: { text: 'Relative depth (simulated)', font: { size: 11 } }, tickformat: ',d' },
    yaxis:      { title: { text: 'Expected species', font: { size: 11 } } },
    showlegend: true,
    legend:     { orientation: 'h', y: -0.25, x: 0 },
  })

  return (
    <ChartCard
      title="Rarefaction Curves"
      subtitle="Simulated rarefaction · relative depth, not actual reads · bold lines = group means"
      onExpand={onExpand}
      onDownload={download}
    >
      <Plot
        ref={plotRef}
        data={traces}
        layout={layout}
        config={baseConfig}
        style={{ width: '100%', height: '380px' }}
        useResizeHandler
      />
    </ChartCard>
  )
}
