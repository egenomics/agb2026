/**
 * components/Sidebar.tsx
 * Left sidebar: logo, upload slots, filters.
 */

import { useAppStore } from '@/store/appStore'
import { useQIIME2 }   from '@/hooks/useQIIME2'
import { UploadSlot }  from './UploadSlot'
import styles from './Sidebar.module.css'

const SLOTS = [
  { id: 'feat'  as const, icon: '📊', label: 'feature-table.tsv', desc: 'ASV counts per sample',          required: true  },
  { id: 'tax'   as const, icon: '🔬', label: 'taxonomy.tsv',       desc: 'ASV → family classification',   required: true  },
  { id: 'meta'  as const, icon: '📋', label: 'metadata.tsv',        desc: 'Sample groups & timepoints',    required: true  },
  { id: 'alpha' as const, icon: '📈', label: 'alpha-diversity.tsv', desc: 'Shannon, Faith PD, Pielou',     required: false },
  { id: 'dist'  as const, icon: '🌐', label: 'distance-matrix.tsv', desc: 'Beta diversity matrix',         required: false },
]

export function Sidebar() {
  const {
    result, isDemo, loadDemo, clearData,
    activeGroups, activeTaxa,
    toggleGroup, toggleTaxon,
    error,
  } = useAppStore()

  const { handleFile, slots } = useQIIME2()

  const groups = result?.groups ?? []
  const taxa   = result?.taxa   ?? []

  return (
    <aside className={styles.sidebar}>
      {/* Logo */}
      <div className={styles.logo}>
        <span>🍦</span>
        <h2>MicroSee</h2>
      </div>

      {/* Data pills */}
      {result && (
        <div className={styles.pills}>
          <span className={styles.pill}>{result.n_samples} samples</span>
          <span className={styles.pill}>{groups.length} groups</span>
          <span className={styles.pill}>{result.n_taxa} taxa</span>
          {isDemo && <span className={`${styles.pill} ${styles.demoTag}`}>demo</span>}
        </div>
      )}

      {/* Upload section */}
      <div className={styles.section}>
        <h3 className={styles.sectionTitle}>QIIME2 Files</h3>
        {SLOTS.map((s) => (
          <UploadSlot
            key={s.id}
            {...s}
            state={slots[s.id]}
            onFile={handleFile}
          />
        ))}
        <div className={styles.hint}>* required · click or drag & drop</div>
      </div>

      {/* Error */}
      {error && <div className={styles.error}>{error}</div>}

      {/* Actions */}
      <div className={styles.actions}>
        <button className={styles.btnDemo}  onClick={loadDemo}>↺ Demo</button>
        {result && (
          <button className={styles.btnClear} onClick={clearData}>✕ Clear</button>
        )}
      </div>

      {/* Group filter */}
      {result && (
        <div className={styles.section}>
          <h3 className={styles.sectionTitle}>Groups</h3>
          <div className={styles.checkList}>
            {groups.map((g) => (
              <label key={g} className={styles.checkItem}>
                <input
                  type="checkbox"
                  checked={activeGroups.has(g)}
                  onChange={() => toggleGroup(g)}
                />
                <span>{g}</span>
              </label>
            ))}
          </div>
        </div>
      )}

      {/* Taxa filter */}
      {result && (
        <div className={styles.section}>
          <h3 className={styles.sectionTitle}>Taxa</h3>
          <div className={styles.checkList}>
            {taxa.map((t) => (
              <label key={t} className={styles.checkItem}>
                <input
                  type="checkbox"
                  checked={activeTaxa.has(t)}
                  onChange={() => toggleTaxon(t)}
                />
                <span className={styles.taxonLabel}>{t}</span>
              </label>
            ))}
          </div>
        </div>
      )}
    </aside>
  )
}
