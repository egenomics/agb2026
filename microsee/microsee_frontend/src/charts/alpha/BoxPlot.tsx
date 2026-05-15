/**
 * charts/alpha/BoxPlot.tsx
 *
 * Box plot with all points shown.
 * Plotly's built-in box is much better than the custom canvas version.
 * Replaces buildBox() from MicroSee.html.
 */

import { useState } from 'react'
import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseLayout, baseConfig, groupColor, uniqueGroups,
  usePlotlyDownload, withAlpha,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  onExpand?: () => void
}

export function BoxPlot({ rows, onExpand }: Props) {
  const [metric, setMetric] = useState<'shannon' | 'simpson'>('shannon')
  const { plotRef, download } = usePlotlyDownload(`box-plot-${metric}`)
  const groups = uniqueGroups(rows)
  const metricLabel = metric === 'shannon' ? "Shannon H′" : 'Simpson'

  const traces = groups.map((g) => {
    const gd = rows.filter((r) => r.group === g)
    const c  = groupColor(g, groups)
    return {
      type:          'box'       as const,
      name:          g,
      y:             gd.map((r) => Number(r[metric])),
      marker:        { color: c, size: 5 },
      fillcolor:     withAlpha(c, 0.25),
      line:          { color: c, width: 1.5 },
      boxmean:       true,
      boxpoints:     'all'       as const,
      jitter:        0.3,
      pointpos:      0,
      text:          gd.map((r) => r.sample_id as string),
      hovertemplate: `<b>%{text}</b>: %{y:.3f}<extra></extra>`,
    }
  })

  const layout = baseLayout({
    height:     360,
    showlegend: false,
    yaxis: { title: { text: metricLabel, font: { size: 11 } } },
  })

  return (
    <ChartCard
      title="Box Plot"
      subtitle="IQR box · whiskers 1.5× IQR · mean shown"
      onExpand={onExpand}
      onDownload={download}
      controls={
        <select
          value={metric}
          onChange={(e) => setMetric(e.target.value as 'shannon' | 'simpson')}
          style={{ fontSize: '12px', padding: '2px 4px' }}
        >
          <option value="shannon">Shannon H′</option>
          <option value="simpson">Simpson</option>
        </select>
      }
    >
      <Plot
        ref={plotRef}
        data={traces}
        layout={layout}
        config={baseConfig}
        style={{ width: '100%', height: '360px' }}
        useResizeHandler
      />
    </ChartCard>
  )
}
