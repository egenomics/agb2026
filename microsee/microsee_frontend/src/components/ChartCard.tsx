/**
 * components/ChartCard.tsx
 * Reusable card wrapper for every chart.
 * Handles: title, subtitle, ⤢ expand button, ↓ download button.
 * Phase 3 will add the interactive flashcard modal.
 */

import type { ReactNode } from 'react'
import styles from './ChartCard.module.css'

interface Props {
  title:     string
  subtitle:  string
  children:  ReactNode
  onExpand?: () => void
  onDownload?: () => void
  controls?: ReactNode  // dropdowns, selects, etc.
  fullWidth?: boolean
}

export function ChartCard({
  title, subtitle, children,
  onExpand, onDownload, controls,
  fullWidth = false,
}: Props) {
  return (
    <div className={`${styles.card} ${fullWidth ? styles.full : ''}`}>
      <div className={styles.header}>
        <div className={styles.headerLeft}>
          <h2 className={styles.title}>{title}</h2>
          <p  className={styles.subtitle}>{subtitle}</p>
        </div>
        <div className={styles.headerRight}>
          {controls}
          {onExpand && (
            <button className={styles.expandBtn} onClick={onExpand} title="Expand">
              ⤢ Expand
            </button>
          )}
          {onDownload && (
            <button className={styles.dlBtn} onClick={onDownload} title="Download PNG">
              ↓ PNG
            </button>
          )}
        </div>
      </div>
      <div className={styles.body}>{children}</div>
    </div>
  )
}
