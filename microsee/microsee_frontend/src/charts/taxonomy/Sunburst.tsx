/**
 * charts/taxonomy/Sunburst.tsx
 *
 * Two-ring sunburst: inner = groups, outer = taxa per group.
 * Replaces buildSunburst() from MicroSee.html.
 */

import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseConfig, taxonColor, groupColor, uniqueGroups,
  mean, taxonValues, usePlotlyDownload, THEME,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  taxa:     string[]
  onExpand?: () => void
}

export function Sunburst({ rows, taxa, onExpand }: Props) {
  const { plotRef, download } = usePlotlyDownload('sunburst')
  const groups = uniqueGroups(rows)

  const ids:     string[] = ['center']
  const labels:  string[] = ['All']
  const parents: string[] = ['']
  const values:  number[] = [0]
  const colors:  string[] = [THEME.bg]

  groups.forEach((g) => {
    const gRows = rows.filter((r) => r.group === g)
    ids.push(g); labels.push(g); parents.push('center')
    values.push(0); colors.push(groupColor(g, groups))

    taxa.forEach((t) => {
      const v = mean(gRows.map((r) => Number(r[t] ?? 0)))
      if (v < 0.5) return
      ids.push(`${g}__${t}`)
      labels.push(t.replace('aceae', ''))
      parents.push(g)
      values.push(+v.toFixed(1))
      colors.push(taxonColor(t))
    })
  })

  const trace = {
    type:          'sunburst' as const,
    ids, labels, parents, values,
    marker:        { colors },
    branchvalues:  'remainder' as const,
    hovertemplate: '<b>%{label}</b>: %{value:.1f}%<extra></extra>',
    insidetextorientation: 'radial' as const,
  }

  return (
    <ChartCard
      title="Sunburst — Group → Taxa"
      subtitle="Inner ring = groups · outer ring = families · hover cells"
      onExpand={onExpand}
      onDownload={download}
    >
      <Plot
        ref={plotRef}
        data={[trace]}
        layout={{
          paper_bgcolor: THEME.paper,
          font:          { family: THEME.font, color: THEME.text, size: 10 },
          height:        380,
          margin:        { l: 0, r: 0, t: 0, b: 0 },
          showlegend:    false,
        }}
        config={baseConfig}
        style={{ width: '100%', height: '380px' }}
        useResizeHandler
      />
    </ChartCard>
  )
}
