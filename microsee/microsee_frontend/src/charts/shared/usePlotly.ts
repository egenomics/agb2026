/**
 * charts/shared/usePlotly.ts
 *
 * Shared utilities for all Plotly chart components:
 *  - theme config (same palette as MicroSee.html)
 *  - colour helpers
 *  - download hook
 *  - common layout defaults
 */

import { useCallback, useRef } from 'react'
import type { Layout, Config, PlotData } from 'plotly.js'
import type { SampleRow } from '@/types/sample'
import { GROUP_COLORS, TAXA_COLORS } from '@/types/sample'

// ── Theme constants ────────────────────────────────────────────────────────

export const THEME = {
  bg:       '#FEF3EC',
  paper:    '#FFFFFF',
  grid:     'rgba(196,160,140,0.12)',
  text:     '#6B3A2A',
  text2:    '#8B5860',
  text3:    '#C4A0A8',
  font:     'Nunito, system-ui, sans-serif',
} as const

/** Base layout applied to every chart */
export const baseLayout = (overrides: Partial<Layout> = {}): Partial<Layout> => ({
  autosize:      true,
  paper_bgcolor: THEME.paper,
  plot_bgcolor:  THEME.bg,
  font:          { family: THEME.font, color: THEME.text, size: 11 },
  margin:        { l: 52, r: 16, t: 32, b: 40 },
  hoverlabel:    { bgcolor: THEME.paper, font: { size: 11, family: THEME.font } },
  xaxis: {
    gridcolor:     THEME.grid,
    linecolor:     THEME.grid,
    zerolinecolor: THEME.grid,
    tickfont:      { size: 10, color: THEME.text2 },
  },
  yaxis: {
    gridcolor:     THEME.grid,
    linecolor:     THEME.grid,
    zerolinecolor: THEME.grid,
    tickfont:      { size: 10, color: THEME.text2 },
  },
  legend: {
    font:            { size: 10 },
    bgcolor:         'rgba(0,0,0,0)',
    orientation:     'h',
    yanchor:         'bottom',
    y:               -0.28,
    xanchor:         'left',
    x:               0,
  },
  ...overrides,
})

/** Base Plotly config — mode bar without logo */
export const baseConfig: Partial<Config> = {
  displaylogo:               false,
  responsive:                true,
  modeBarButtonsToRemove:    ['lasso2d', 'select2d', 'autoScale2d'],
  toImageButtonOptions:      { format: 'png', scale: 2 },
}

// ── Colour helpers ──────────────────────────────────────────────────────────

export function groupColor(group: string, allGroups: string[]): string {
  const idx = allGroups.sort().indexOf(group)
  return GROUP_COLORS[idx % GROUP_COLORS.length] ?? '#aaa'
}

export function taxonColor(taxon: string): string {
  return TAXA_COLORS[taxon] ?? '#aaa'
}

/** Hex colour → rgba string with given opacity (0–1) */
export function withAlpha(hex: string, alpha: number): string {
  const r = parseInt(hex.slice(1, 3), 16)
  const g = parseInt(hex.slice(3, 5), 16)
  const b = parseInt(hex.slice(5, 7), 16)
  return `rgba(${r},${g},${b},${alpha})`
}

// ── Data helpers ────────────────────────────────────────────────────────────

/** Mean of a numeric array */
export function mean(arr: number[]): number {
  return arr.length ? arr.reduce((a, b) => a + b, 0) / arr.length : 0
}

/** Relative abundance of a taxon across rows (already in %, just return) */
export function taxonValues(rows: SampleRow[], taxon: string): number[] {
  return rows.map((r) => Number(r[taxon] ?? 0))
}

/** All unique groups from rows, sorted */
export function uniqueGroups(rows: SampleRow[]): string[] {
  return [...new Set(rows.map((r) => r.group))].sort()
}

/** Normalise taxon columns so each row sums to 100% */
export function normaliseRows(rows: SampleRow[], taxa: string[]): SampleRow[] {
  return rows.map((r) => {
    const tot = taxa.reduce((s, t) => s + Number(r[t] ?? 0), 0) || 1
    const out = { ...r }
    taxa.forEach((t) => { out[t] = (Number(r[t] ?? 0) / tot) * 100 })
    return out
  })
}

// ── Statistics helpers ──────────────────────────────────────────────────────

/** Lanczos-approximated log-gamma (accurate to ~15 sig figs) */
function lgamma(z: number): number {
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z)
  z -= 1
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
  let x = c[0]
  for (let i = 1; i < 9; i++) x += c[i] / (z + i)
  const t = z + 7.5
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x)
}

/** Lentz continued-fraction kernel for regularised incomplete beta */
function ibetacf(x: number, a: number, b: number): number {
  const fpMin = 1e-30, eps = 3e-7
  let c = 1, d = 1 - (a + b) * x / (a + 1)
  d = Math.abs(d) < fpMin ? fpMin : d; d = 1 / d; let h = d
  for (let m = 1; m <= 200; m++) {
    const m2 = 2 * m
    let aa = m * (b - m) * x / ((a + m2 - 1) * (a + m2))
    d = 1 + aa * d; d = Math.abs(d) < fpMin ? fpMin : d
    c = 1 + aa / c; c = Math.abs(c) < fpMin ? fpMin : c
    d = 1 / d; h *= d * c
    aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1))
    d = 1 + aa * d; d = Math.abs(d) < fpMin ? fpMin : d
    c = 1 + aa / c; c = Math.abs(c) < fpMin ? fpMin : c
    d = 1 / d; const delta = d * c; h *= delta
    if (Math.abs(delta - 1) < eps) break
  }
  return h
}

/** Regularised incomplete beta I(x; a, b) — Numerical Recipes implementation */
function ibeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0
  if (x >= 1) return 1
  const logbt = a * Math.log(x) + b * Math.log(1 - x) - lgamma(a) - lgamma(b) + lgamma(a + b)
  if (x < (a + 1) / (a + b + 2)) return Math.exp(logbt) * ibetacf(x, a, b) / a
  return 1 - Math.exp(logbt) * ibetacf(1 - x, b, a) / b
}

/**
 * Exact two-tailed p-value for a t-statistic with given degrees of freedom.
 * Uses I(df/(df+t²), df/2, 0.5) — the exact t-distribution CDF via incomplete beta.
 */
export function twoTailP(t: number, df: number): number {
  if (df <= 0 || t === 0) return 1
  return Math.max(0.001, Math.min(1, ibeta(df / (df + t * t), df / 2, 0.5)))
}

// ── Download hook ───────────────────────────────────────────────────────────

/**
 * Returns a ref to attach to a <Plot> element and a download function.
 * Usage:
 *   const { plotRef, download } = usePlotlyDownload('my-chart')
 *   <Plot ref={plotRef} ... />
 *   <button onClick={download}>↓ PNG</button>
 */
export function usePlotlyDownload(filename: string) {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const plotRef = useRef<any>(null)

  const download = useCallback(async () => {
    if (!plotRef.current) return
    const Plotly = (await import('plotly.js-dist-min')).default
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    await Plotly.downloadImage(plotRef.current.el as any, {
      format:   'png',
      width:    1200,
      height:   600,
      filename: filename,
      scale:    2,
    })
  }, [filename])

  return { plotRef, download }
}
