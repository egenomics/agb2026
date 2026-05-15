/**
 * components/StatsBar.tsx
 * Top bar showing sample count, groups, mean Shannon.
 */

import { useAppStore } from '@/store/appStore'
import styles from './StatsBar.module.css'

export function StatsBar() {
  const { result, isDemo, filteredRows } = useAppStore()
  if (!result) return null

  const rows    = filteredRows()
  const shannon = rows.length
    ? (rows.reduce((s, r) => s + (r.shannon as number), 0) / rows.length).toFixed(3)
    : '—'

  return (
    <div className={styles.bar}>
      <div className={styles.left}>
        <span className={styles.title}>MicroSee</span>
        <span className={styles.sub}>
          Advanced Genome Bioinformatics 2026 · UPF · Group D
        </span>
      </div>
      <div className={styles.stats}>
        <Stat value={rows.length}                    label="Samples"      />
        <Stat value={result.groups.length}           label="Groups"       />
        <Stat value={result.n_taxa}                  label="Taxa"         />
        <Stat value={`H′ ${shannon}`}                label="Mean Shannon" />
        {isDemo && <span className={styles.demoTag}>demo data</span>}
      </div>
    </div>
  )
}

function Stat({ value, label }: { value: string | number; label: string }) {
  return (
    <div className={styles.stat}>
      <span className={styles.statVal}>{value}</span>
      <span className={styles.statLbl}>{label}</span>
    </div>
  )
}
