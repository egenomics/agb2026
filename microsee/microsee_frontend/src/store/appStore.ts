/**
 * store/appStore.ts
 *
 * Global app state using Zustand.
 * Replaces the 6 global variables in MicroSee.html:
 *   data, rawTaxa, activeGroups, activeTaxa, CHS, CHART_CONFIGS
 */

import { create } from 'zustand'
import type { IntegrateResult, SampleRow, UploadSlotId, UploadSlotState, SlotStatus } from '@/types/sample'
import { DEMO_RESULT } from '@/data/demo'

type AppScreen = 'splash' | 'landing' | 'results'

interface AppStore {
  // ── Screen state ──────────────────────────────────────────────────────────
  screen:       AppScreen
  setScreen:    (s: AppScreen) => void

  // ── Data ──────────────────────────────────────────────────────────────────
  result:       IntegrateResult | null
  isDemo:       boolean
  loadDemo:     () => void
  loadResult:   (r: IntegrateResult) => void
  clearData:    () => void

  // ── Filters ───────────────────────────────────────────────────────────────
  activeGroups: Set<string>
  activeTaxa:   Set<string>
  toggleGroup:  (g: string) => void
  toggleTaxon:  (t: string) => void
  setAllGroups: (groups: string[]) => void
  setAllTaxa:   (taxa: string[]) => void

  // ── Metric ────────────────────────────────────────────────────────────────
  metric:       'shannon' | 'simpson'
  setMetric:    (m: 'shannon' | 'simpson') => void

  // ── Upload slots ──────────────────────────────────────────────────────────
  slots:        Record<UploadSlotId, UploadSlotState>
  setSlot:      (id: UploadSlotId, patch: Partial<UploadSlotState>) => void
  resetSlots:   () => void

  // ── Error message ─────────────────────────────────────────────────────────
  error:        string | null
  setError:     (e: string | null) => void

  // ── Derived helpers ───────────────────────────────────────────────────────
  filteredRows: () => SampleRow[]
  filteredTaxa: () => string[]
}

const emptySlot = (): UploadSlotState => ({
  status:   'idle',
  filename: null,
  file:     null,
})

const SLOT_IDS: UploadSlotId[] = ['feat', 'tax', 'meta', 'alpha', 'dist']

export const useAppStore = create<AppStore>((set, get) => ({
  // ── Screen ────────────────────────────────────────────────────────────────
  screen:    'splash',
  setScreen: (screen) => set({ screen }),

  // ── Data ──────────────────────────────────────────────────────────────────
  result: null,
  isDemo: false,

  loadDemo: () => {
    set({
      result:       DEMO_RESULT,
      isDemo:       true,
      screen:       'results',
      activeGroups: new Set(DEMO_RESULT.groups),
      activeTaxa:   new Set(DEMO_RESULT.taxa),
      error:        null,
    })
  },

  loadResult: (result) => {
    set({
      result,
      isDemo:       false,
      screen:       'results',
      activeGroups: new Set(result.groups),
      activeTaxa:   new Set(result.taxa),
      error:        null,
    })
  },

  clearData: () => {
    set({
      result:  null,
      isDemo:  false,
      screen:  'landing',
      error:   null,
    })
    get().resetSlots()
  },

  // ── Filters ───────────────────────────────────────────────────────────────
  activeGroups: new Set<string>(),
  activeTaxa:   new Set<string>(),

  toggleGroup: (g) =>
    set((state) => {
      const next = new Set(state.activeGroups)
      next.has(g) ? next.delete(g) : next.add(g)
      return { activeGroups: next }
    }),

  toggleTaxon: (t) =>
    set((state) => {
      const next = new Set(state.activeTaxa)
      next.has(t) ? next.delete(t) : next.add(t)
      return { activeTaxa: next }
    }),

  setAllGroups: (groups) => set({ activeGroups: new Set(groups) }),
  setAllTaxa:   (taxa)   => set({ activeTaxa: new Set(taxa) }),

  // ── Metric ────────────────────────────────────────────────────────────────
  metric:    'shannon',
  setMetric: (metric) => set({ metric }),

  // ── Upload slots ──────────────────────────────────────────────────────────
  slots: Object.fromEntries(SLOT_IDS.map((id) => [id, emptySlot()])) as Record<
    UploadSlotId,
    UploadSlotState
  >,

  setSlot: (id, patch) =>
    set((state) => ({
      slots: { ...state.slots, [id]: { ...state.slots[id], ...patch } },
    })),

  resetSlots: () =>
    set({
      slots: Object.fromEntries(
        SLOT_IDS.map((id) => [id, emptySlot()])
      ) as Record<UploadSlotId, UploadSlotState>,
    }),

  // ── Error ─────────────────────────────────────────────────────────────────
  error:    null,
  setError: (error) => set({ error }),

  // ── Derived ───────────────────────────────────────────────────────────────
  filteredRows: () => {
    const { result, activeGroups } = get()
    if (!result) return []
    return result.rows.filter((r) => activeGroups.has(r.group))
  },

  filteredTaxa: () => {
    const { result, activeTaxa } = get()
    if (!result) return []
    return result.taxa.filter((t) => activeTaxa.has(t))
  },
}))
