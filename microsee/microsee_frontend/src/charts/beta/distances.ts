/**
 * charts/beta/distances.ts
 *
 * Pure TypeScript implementations of Bray-Curtis and Jaccard distances,
 * PCoA (classical MDS), and a stress-minimising NMDS approximation.
 *
 * These run in the browser for demo data.
 * Real QIIME2 data uses pre-computed distance matrices from the backend.
 */

import type { SampleRow } from '@/types/sample'

// ── Distance functions ──────────────────────────────────────────────────────

export function brayCurtis(a: SampleRow, b: SampleRow, taxa: string[]): number {
  let num = 0, den = 0
  for (const t of taxa) {
    const ai = Number(a[t] ?? 0)
    const bi = Number(b[t] ?? 0)
    num += Math.abs(ai - bi)
    den += ai + bi
  }
  return den > 0 ? num / den : 0
}

export function jaccard(a: SampleRow, b: SampleRow, taxa: string[], threshold = 5): number {
  const totA = taxa.reduce((s, t) => s + Number(a[t] ?? 0), 0) || 1
  const totB = taxa.reduce((s, t) => s + Number(b[t] ?? 0), 0) || 1
  let both = 0, either = 0
  for (const t of taxa) {
    const pa = (Number(a[t] ?? 0) / totA * 100) >= threshold
    const pb = (Number(b[t] ?? 0) / totB * 100) >= threshold
    if (pa || pb) { either++; if (pa && pb) both++ }
  }
  return either > 0 ? 1 - both / either : 0
}

// ── Distance matrix ─────────────────────────────────────────────────────────

export function buildDistMatrix(
  rows: SampleRow[],
  taxa: string[],
  distFn: (a: SampleRow, b: SampleRow, taxa: string[]) => number,
): number[][] {
  const n = rows.length
  return Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) =>
      i === j ? 0 : distFn(rows[i], rows[j], taxa),
    ),
  )
}

// ── PCoA (classical MDS) ────────────────────────────────────────────────────

export interface OrdPoint {
  sample_id: string
  x: number
  y: number
  pct1: number
  pct2: number
}

export function pcoa(rows: SampleRow[], matrix: number[]): OrdPoint[] {
  const n = rows.length
  if (n < 3) return rows.map(r => ({ sample_id: String(r.sample_id), x: 0, y: 0, pct1: 0, pct2: 0 }))

  // Double-center
  const D2 = matrix.map((v) => v * v)
  const rm = Array.from({ length: n }, (_, i) => {
    let s = 0; for (let j = 0; j < n; j++) s += D2[i * n + j]; return s / n
  })
  const cm = Array.from({ length: n }, (_, j) => {
    let s = 0; for (let i = 0; i < n; i++) s += D2[i * n + j]; return s / n
  })
  const gm = rm.reduce((a, b) => a + b, 0) / n
  const B  = Array.from({ length: n * n }, (_, k) => {
    const i = Math.floor(k / n), j = k % n
    return -0.5 * (D2[k] - rm[i] - cm[j] + gm)
  })

  // Power iteration for top 2 eigenvectors
  function eigvec(perp: number[] | null): { v: number[]; lam: number } {
    let v = Array.from({ length: n }, (_, i) => Math.sin(i * 1.7 + 0.5))
    const orth = (u: number[]) => {
      if (!perp) return u
      const d = u.reduce((s, x, i) => s + x * perp[i], 0)
      return u.map((x, i) => x - d * perp[i])
    }
    v = orth(v)
    let nm = Math.sqrt(v.reduce((s, x) => s + x * x, 0)) || 1
    v = v.map(x => x / nm)
    for (let it = 0; it < 300; it++) {
      const nv_raw = Array.from({ length: n }, (_, i) =>
        v.reduce((s, vj, j) => s + B[i * n + j] * vj, 0),
      )
      const nv = orth(nv_raw)
      nm = Math.sqrt(nv.reduce((s, x) => s + x * x, 0)) || 1
      const nvu = nv.map(x => x / nm)
      if (Math.sqrt(nvu.reduce((s, x, i) => s + (x - v[i]) ** 2, 0)) < 1e-10) { Object.assign(v, nvu); break }
      Object.assign(v, nvu)
    }
    const lam = v.reduce((s, vi, i) => s + vi * B.slice(i * n, i * n + n).reduce((ss, b, j) => ss + b * v[j], 0), 0)
    return { v, lam }
  }

  const e1 = eigvec(null)
  const e2 = eigvec(e1.v)
  const s1 = Math.sqrt(Math.max(e1.lam, 0))
  const s2 = Math.sqrt(Math.max(e2.lam, 0))
  // Total inertia = trace(B) = sum of ALL eigenvalues (positive + negative).
  // Using trace as denominator gives accurate % even without computing all eigenvectors.
  const traceB = Array.from({ length: n }, (_, i) => B[i * n + i]).reduce((s, v) => s + v, 0)
  const tot = Math.max(traceB, Math.max(e1.lam, 0) + Math.max(e2.lam, 0) + 1e-10)
  const pct1 = +((Math.max(e1.lam, 0) / tot) * 100).toFixed(1)
  const pct2 = +((Math.max(e2.lam, 0) / tot) * 100).toFixed(1)

  return rows.map((r, i) => ({
    sample_id: String(r.sample_id),
    x: +(e1.v[i] * s1).toFixed(4),
    y: +(e2.v[i] * s2).toFixed(4),
    pct1, pct2,
  }))
}

export function computePCoA(rows: SampleRow[], taxa: string[], distFn: typeof brayCurtis): OrdPoint[] {
  const mat = buildDistMatrix(rows, taxa, distFn)
  const flat = mat.flat()
  return pcoa(rows, flat)
}
