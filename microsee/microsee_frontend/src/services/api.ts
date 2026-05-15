/**
 * services/api.ts
 *
 * Typed wrappers around the FastAPI backend.
 * Every function returns a typed result or throws an ApiError.
 *
 * Usage:
 *   const result = await api.integrate(featFile, taxFile, metaFile)
 */

import type { IntegrateResult, ParseError } from '@/types/sample'

const BASE = '/api'  // proxied to http://localhost:8000 via vite.config.ts

// ── Error type ────────────────────────────────────────────────────────────────

export class ApiError extends Error {
  constructor(
    public status:     number,
    public detail:     ParseError | string,
    public file_type?: string,
  ) {
    const msg =
      typeof detail === 'string'
        ? detail
        : `[${detail.file_type}] ${detail.message}`
    super(msg)
    this.name = 'ApiError'
  }
}

async function handleResponse<T>(res: Response): Promise<T> {
  if (res.ok) return res.json() as Promise<T>

  let detail: ParseError | string
  try {
    const body = await res.json()
    detail = body.detail ?? body
  } catch {
    detail = await res.text()
  }
  throw new ApiError(res.status, detail)
}

// ── Endpoints ─────────────────────────────────────────────────────────────────

export const api = {

  /**
   * Health check — use to verify the backend is running.
   */
  async health(): Promise<{ status: string; version: string }> {
    const res = await fetch(`${BASE.replace('/api', '')}/health`)
    return handleResponse(res)
  },

  /**
   * Full integration endpoint.
   * Sends all uploaded files and returns chart-ready SampleRows.
   */
  async integrate(
    featureTable:    File,
    taxonomy:        File,
    metadata:        File,
    alphaDiversity?: File | null,
    distanceMatrix?: File | null,
  ): Promise<IntegrateResult> {
    const form = new FormData()
    form.append('feature_table',   featureTable)
    form.append('taxonomy',        taxonomy)
    form.append('metadata',        metadata)
    if (alphaDiversity) form.append('alpha_diversity', alphaDiversity)
    if (distanceMatrix) form.append('distance_matrix', distanceMatrix)

    const res = await fetch(`${BASE}/parse/integrate`, {
      method: 'POST',
      body:   form,
    })
    return handleResponse<IntegrateResult>(res)
  },

  /**
   * Parse feature table only (for validation / preview).
   */
  async parseFeatureTable(file: File) {
    const form = new FormData()
    form.append('file', file)
    const res = await fetch(`${BASE}/parse/feature-table`, {
      method: 'POST', body: form,
    })
    return handleResponse(res)
  },

  /**
   * Parse taxonomy only (for validation / preview).
   */
  async parseTaxonomy(file: File) {
    const form = new FormData()
    form.append('file', file)
    const res = await fetch(`${BASE}/parse/taxonomy`, {
      method: 'POST', body: form,
    })
    return handleResponse(res)
  },

  /**
   * Parse metadata only (for validation / preview).
   */
  async parseMetadata(file: File) {
    const form = new FormData()
    form.append('file', file)
    const res = await fetch(`${BASE}/parse/metadata`, {
      method: 'POST', body: form,
    })
    return handleResponse(res)
  },
}
