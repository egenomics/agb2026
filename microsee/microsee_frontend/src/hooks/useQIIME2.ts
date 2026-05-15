/**
 * hooks/useQIIME2.ts
 *
 * Manages the QIIME2 file upload flow.
 * Tracks per-slot state and triggers integration when all required files present.
 */

import { useCallback } from 'react'
import { useAppStore } from '@/store/appStore'
import { api, ApiError }  from '@/services/api'
import type { UploadSlotId } from '@/types/sample'

const REQUIRED: UploadSlotId[] = ['feat', 'tax', 'meta']

export function useQIIME2() {
  const { slots, setSlot, setError, loadResult } = useAppStore()

  const handleFile = useCallback(
    async (id: UploadSlotId, file: File) => {
      setSlot(id, { status: 'loading', filename: file.name, file })
      setError(null)

      // Optional slots just store the file, no individual validation call
      if (!REQUIRED.includes(id)) {
        setSlot(id, { status: 'ok' })
        return
      }

      // For required slots, quickly validate by calling the individual endpoint
      try {
        if (id === 'feat')  await api.parseFeatureTable(file)
        if (id === 'tax')   await api.parseTaxonomy(file)
        if (id === 'meta')  await api.parseMetadata(file)
        setSlot(id, { status: 'ok' })
      } catch (err) {
        const msg = err instanceof ApiError ? err.message : 'Unknown parse error'
        setSlot(id, { status: 'error' })
        setError(msg)
        return
      }

      // Check if all required files are now OK — if so, integrate
      const updatedSlots = { ...slots, [id]: { ...slots[id], status: 'ok', file } }
      const allReady     = REQUIRED.every((r) => updatedSlots[r].status === 'ok' && updatedSlots[r].file)

      if (allReady) {
        await triggerIntegration(updatedSlots, setError, loadResult, setSlot)
      }
    },
    [slots, setSlot, setError, loadResult],
  )

  const requiredReady = REQUIRED.every((r) => slots[r].status === 'ok')
  const anyLoading    = Object.values(slots).some((s) => s.status === 'loading')

  return { handleFile, slots, requiredReady, anyLoading }
}

async function triggerIntegration(
  slots:      ReturnType<typeof useAppStore.getState>['slots'],
  setError:   (e: string | null) => void,
  loadResult: ReturnType<typeof useAppStore.getState>['loadResult'],
  setSlot:    ReturnType<typeof useAppStore.getState>['setSlot'],
) {
  const feat  = slots.feat.file!
  const tax   = slots.tax.file!
  const meta  = slots.meta.file!
  const alpha = slots.alpha.file ?? null
  const dist  = slots.dist.file  ?? null

  try {
    const result = await api.integrate(feat, tax, meta, alpha, dist)

    if (result.warnings.length > 0) {
      setError(`Loaded with warnings: ${result.warnings.join(' | ')}`)
    }

    loadResult(result)
  } catch (err) {
    const msg = err instanceof ApiError ? err.message : 'Integration failed'
    setError(msg)
  }
}
