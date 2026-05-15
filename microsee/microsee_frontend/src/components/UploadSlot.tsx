/**
 * components/UploadSlot.tsx
 * Single QIIME2 file upload slot with drag-and-drop.
 */

import { useRef, type DragEvent, type ChangeEvent } from 'react'
import type { UploadSlotId, UploadSlotState } from '@/types/sample'
import styles from './UploadSlot.module.css'

interface Props {
  id:       UploadSlotId
  icon:     string
  label:    string
  desc:     string
  required: boolean
  state:    UploadSlotState
  onFile:   (id: UploadSlotId, file: File) => void
  accept?:  string
}

export function UploadSlot({ id, icon, label, desc, required, state, onFile, accept = '.tsv,.txt' }: Props) {
  const inputRef = useRef<HTMLInputElement>(null)

  const handleClick = () => inputRef.current?.click()

  const handleChange = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (file) onFile(id, file)
    // reset so same file can be re-uploaded
    e.target.value = ''
  }

  const handleDragOver = (e: DragEvent) => {
    e.preventDefault()
    e.currentTarget.classList.add(styles.dragOver)
  }

  const handleDragLeave = (e: DragEvent) => {
    e.currentTarget.classList.remove(styles.dragOver)
  }

  const handleDrop = (e: DragEvent) => {
    e.preventDefault()
    e.currentTarget.classList.remove(styles.dragOver)
    const file = e.dataTransfer.files?.[0]
    if (file) onFile(id, file)
  }

  const statusIcon: Record<string, string> = {
    idle:    '—',
    loading: '⋯',
    ok:      '✓',
    error:   '✗',
  }

  const statusColor: Record<string, string> = {
    idle:    'var(--text3)',
    loading: 'var(--text2)',
    ok:      '#2FA896',
    error:   '#C03030',
  }

  return (
    <div
      className={`${styles.slot} ${state.status === 'ok' ? styles.loaded : ''}`}
      onClick={handleClick}
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
      title={state.filename ?? `Drop ${label} here`}
      role="button"
      tabIndex={0}
      onKeyDown={(e) => e.key === 'Enter' && handleClick()}
    >
      <span className={styles.icon}>{icon}</span>
      <div className={styles.info}>
        <span className={styles.name}>
          {label}
          {required && <span className={styles.req}>*</span>}
        </span>
        <span className={styles.descText}>
          {state.filename ?? desc}
        </span>
      </div>
      <span
        className={styles.status}
        style={{ color: statusColor[state.status] }}
      >
        {statusIcon[state.status]}
      </span>
      <input
        ref={inputRef}
        type="file"
        accept={accept}
        style={{ display: 'none' }}
        onChange={handleChange}
      />
    </div>
  )
}
