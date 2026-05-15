/**
 * charts/beta/DeltaHeatmap.tsx
 * Change in relative abundance (T84 − T0) per taxon per patient.
 * Red = increased, Blue = decreased, White = stable.
 */

import { useMemo } from 'react'
import Plot from 'react-plotly.js'
import { ChartCard } from '@/components/ChartCard'
import { baseConfig, usePlotlyDownload, THEME } from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props { rows: SampleRow[]; taxa: string[]; onExpand?: () => void }

export function DeltaHeatmap({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('delta-heatmap')

  const { patients, sortedTaxa, zMatrix, textMatrix } = useMemo(() => {
    const patientIds = [...new Set(rows.map(r => String(r.patient ?? r.sample_id).replace(/_T\d+$/, '')))]
    const deltas = patientIds.map(pid => {
      const r0  = rows.find(r => (String(r.patient ?? r.sample_id).replace(/_T\d+$/, '') === pid) && Number(r.time) === 0)
      const r84 = rows.filter(r => (String(r.patient ?? r.sample_id).replace(/_T\d+$/, '') === pid) && Number(r.time) > 0).sort((a, b) => Number(b.time) - Number(a.time))[0]
      if (!r0 || !r84) return null
      const vals: Record<string, number> = {}
      const tot0  = taxa.reduce((s, t) => s + Number(r0[t]  ?? 0), 0) || 1
      const tot84 = taxa.reduce((s, t) => s + Number(r84[t] ?? 0), 0) || 1
      taxa.forEach(t => {
        vals[t] = +((Number(r84[t] ?? 0) / tot84 * 100) - (Number(r0[t] ?? 0) / tot0 * 100)).toFixed(1)
      })
      return { pid, vals, group: String(r0.base_group ?? r0.group).replace(/_T\d+$/, '') }
    }).filter(Boolean) as { pid: string; vals: Record<string, number>; group: string }[]

    const sorted = [...deltas].sort((a, b) => a.group.localeCompare(b.group))
    const sortedTaxa = [...taxa].sort((a, b) =>
      Math.max(...sorted.map(d => Math.abs(d.vals[b] ?? 0))) -
      Math.max(...sorted.map(d => Math.abs(d.vals[a] ?? 0)))
    )

    const zMatrix = sortedTaxa.map(t => sorted.map(d => d.vals[t] ?? 0))
    const textMatrix = sortedTaxa.map(t => sorted.map(d => `${d.pid}<br>${t}: ${d.vals[t] > 0 ? '+' : ''}${d.vals[t]}%`))

    return { patients: sorted.map(d => d.pid), sortedTaxa, zMatrix, textMatrix }
  }, [rows, taxa])

  return (
    <ChartCard title="Δ Abundance Heatmap — T0 → T84"
      subtitle="Red = family increased · Blue = decreased · White = stable · numbers = Δ%"
      onExpand={onExpand} onDownload={download} fullWidth>
      <Plot ref={plotRef}
        data={[{
          type:          'heatmap',
          z:             zMatrix, x: patients, y: sortedTaxa,
          text:          textMatrix,
          hovertemplate: '%{text}<extra></extra>',
          colorscale: [
            [0.0,  '#2050A0'], [0.35, '#A0C0F8'],
            [0.5,  '#FFF8F4'],
            [0.65, '#F8C0A0'], [1.0,  '#C03030'],
          ],
          zmid: 0, zmin: -12, zmax: 12,
          colorbar: { title: { text: 'Δ%' }, len: 0.8, thickness: 14 },
          showscale: true,
          xgap: 1.5,
          ygap: 1.5,
        } as any]}
        layout={{
          autosize: true,
          paper_bgcolor: THEME.paper, plot_bgcolor: '#000000',
          font: { family: THEME.font, color: THEME.text, size: 10 },
          height: 400,
          margin: { l: 140, r: 80, t: 28, b: 100 },
          xaxis: { tickangle: -35, tickfont: { size: 9 } },
          yaxis: { tickfont: { size: 9 } },
        }}
        config={baseConfig} style={{ width: '100%', height: '400px' }} useResizeHandler />
    </ChartCard>
  )
}
