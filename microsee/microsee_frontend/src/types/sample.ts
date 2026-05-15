/**
 * types/sample.ts
 *
 * TypeScript interfaces that exactly mirror the Pydantic models
 * in backend/app/models/sample.py.
 *
 * When a backend model changes, update it here too.
 */

export const TAXA = [
  'Bacteroidaceae',
  'Lachnospiraceae',
  'Ruminococcaceae',
  'Prevotellaceae',
  'Rikenellaceae',
  'Enterobacteriaceae',
  'Oscillospiraceae',
  'Tannerellaceae',
  'Akkermansiaceae',
] as const

export type TaxonName = (typeof TAXA)[number]

/** Taxa colour map — same order as TC in MicroSee.html */
export const TAXA_COLORS: Record<string, string> = {
  Bacteroidaceae:     '#4A7ED4',
  Lachnospiraceae:    '#D84E6A',
  Ruminococcaceae:    '#2FA896',
  Prevotellaceae:     '#D97A3A',
  Rikenellaceae:      '#9058C4',
  Enterobacteriaceae: '#C4960A',
  Oscillospiraceae:   '#8CAD70',
  Tannerellaceae:     '#E87AA0',
  Akkermansiaceae:    '#5B8C5A',
}

/** Group colour palette — same order as GPAL in MicroSee.html */
export const GROUP_COLORS = [
  '#D97A3A',
  '#C4960A',
  '#4A7ED4',
  '#89BAF0',
  '#2FA896',
  '#9058C4',
]

export function groupColor(group: string, allGroups: string[]): string {
  const idx = allGroups.indexOf(group)
  return GROUP_COLORS[idx % GROUP_COLORS.length] ?? '#aaa'
}

// ── API response types ────────────────────────────────────────────────────────

export interface SampleMetadata {
  sample_id:  string
  group:      string
  base_group: string
  timepoint:  string
  time:       number | null
  patient:    string
  sixmwt:     number
  il18:       number
}

/**
 * A fully integrated, chart-ready sample row.
 * Taxon abundances are dynamic extra fields (string → number).
 */
export interface SampleRow extends SampleMetadata {
  shannon:  number
  simpson:  number
  // taxon relative abundances as dynamic keys
  [taxon: string]: number | string | null
}

export interface IntegrateResult {
  rows:         SampleRow[]
  taxa:         string[]
  n_samples:    number
  n_taxa:       number
  groups:       string[]
  has_clinical: boolean
  warnings:     string[]
}

export interface ParseError {
  file_type: string
  message:   string
  line:      number | null
  hint:      string | null
}

// ── Upload state ──────────────────────────────────────────────────────────────

export type UploadSlotId = 'feat' | 'tax' | 'meta' | 'alpha' | 'dist'

export type SlotStatus = 'idle' | 'loading' | 'ok' | 'error'

export interface UploadSlotState {
  status:   SlotStatus
  filename: string | null
  file:     File | null
}
