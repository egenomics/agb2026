/**
 * components/FlashCard.tsx
 *
 * Blurred-background modal that shows an expanded version of any chart.
 * Replaced the PNG-screenshot approach from MicroSee.html with:
 *  - A live re-rendered Plotly figure at full size
 *  - A plain-language "what it shows" description
 *  - A highlighted insight paragraph
 *  - Dynamic stat pills computed from current data
 */

import { useEffect, type ReactNode } from 'react'
import styles from './FlashCard.module.css'

interface FlashCardInfo {
  title:   string
  what:    string
  insight: string
  pills?:  string[]
}

interface Props {
  isOpen:  boolean
  onClose: () => void
  info:    FlashCardInfo
  children: ReactNode  // the live Plotly chart
}

export function FlashCard({ isOpen, onClose, info, children }: Props) {
  // Close on Escape
  useEffect(() => {
    if (!isOpen) return
    const handler = (e: KeyboardEvent) => { if (e.key === 'Escape') onClose() }
    window.addEventListener('keydown', handler)
    document.body.style.overflow = 'hidden'
    return () => {
      window.removeEventListener('keydown', handler)
      document.body.style.overflow = ''
    }
  }, [isOpen, onClose])

  if (!isOpen) return null

  return (
    <div
      className={styles.overlay}
      onClick={(e) => { if (e.target === e.currentTarget) onClose() }}
    >
      <div className={styles.bg} onClick={onClose} />

      <div className={styles.card}>
        {/* Header */}
        <div className={styles.head}>
          <span className={styles.title}>{info.title}</span>
          <button className={styles.close} onClick={onClose}>✕</button>
        </div>

        {/* Live chart */}
        <div className={styles.chartWrap}>
          {children}
        </div>

        {/* Info panel */}
        <div className={styles.info}>
          <p className={styles.what}>{info.what}</p>
          <div className={styles.insight}>💡 {info.insight}</div>
          {info.pills && info.pills.length > 0 && (
            <div className={styles.pills}>
              {info.pills.map((p) => (
                <span key={p} className={styles.pill}>{p}</span>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  )
}

// ── Hook to manage open/close state ──────────────────────────────────────────

import { useState, useCallback } from 'react'

export function useFlashCard() {
  const [open, setOpen] = useState(false)
  const show  = useCallback(() => setOpen(true),  [])
  const close = useCallback(() => setOpen(false), [])
  return { open, show, close }
}
