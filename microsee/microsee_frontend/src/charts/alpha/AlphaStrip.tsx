/**
 * charts/alpha/AlphaStrip.tsx
 *
 * Strip chart: one dot per sample, jittered by group.
 * Mean line per group.
 * Replaces buildAlpha() from MicroSee.html.
 */

import { useState } from 'react'
import Plot            from 'react-plotly.js'
import { ChartCard }  from '@/components/ChartCard'
import {
  baseLayout, baseConfig, groupColor, mean,
  uniqueGroups, usePlotlyDownload, withAlpha,
} from '@/charts/shared/usePlotly'
import type { SampleRow } from '@/types/sample'

interface Props {
  rows:     SampleRow[]
  onExpand?: () => void
}

// Seeded jitter — deterministic so chart doesn't jump on re-render
function jitter(i: number): number {
  return (((i * 1337 + 17) % 100) - 50) / 200
}

export function AlphaStrip({ rows, onExpand }: Props) {
  const [metric, setMetric] = useState<'shannon' | 'simpson'>('shannon')
  const { plotRef, download } = usePlotlyDownload(`alpha-strip-${metric}`)
  const groups = uniqueGroups(rows)
  const metricLabel = metric === 'shannon' ? "Shannon H′" : 'Simpson'

  const traces = groups.map((g) => {
    const gd  = rows.filter((r) => r.group === g)
    const c   = groupColor(g, groups)
    const gi  = groups.indexOf(g)
    const ys  = gd.map((r) => Number(r[metric]))
    const xs  = gd.map((_, i) => gi + jitter(i))
    const avg = mean(ys)

    return [
      // Dots
      {
        type:          'scatter'    as const,
        mode:          'markers'   as const,
        x:             xs,
        y:             ys,
        name:          g,
        marker:        { color: withAlpha(c, 0.85), size: 9, line: { width: 1, color: 'white' } },
        text:          gd.map((r) => r.sample_id as string),
        hovertemplate: `<b>%{text}</b><br>${metricLabel}: %{y:.3f}<extra></extra>`,
        showlegend:    true,
      },
      // Mean line (shape-like scatter)
      {
        type:          'scatter'    as const,
        mode:          'lines'     as const,
        x:             [gi - 0.3, gi + 0.3],
        y:             [avg, avg],
        line:          { color: c, width: 2.5 },
        showlegend:    false,
        hoverinfo:     'skip'      as const,
      },
    ]
  }).flat()

  const layout = baseLayout({
    height:     360,
    showlegend: true,
    xaxis: {
      tickmode:  'array',
      tickvals:  groups.map((_, i) => i),
      ticktext:  groups,
      tickfont:  { size: 10 },
      gridcolor: 'transparent',
      zeroline:  false,
    },
    yaxis: { title: { text: metricLabel, font: { size: 11 } } },
  })

  return (
    <ChartCard
      title="Alpha Diversity — Strip Chart"
      subtitle={`Each dot = one sample · bar = group mean · metric: ${metricLabel}`}
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
