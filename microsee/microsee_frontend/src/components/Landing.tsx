/**
 * components/Landing.tsx
 * Landing page shown after the splash, before data is loaded.
 */

import { useAppStore } from '@/store/appStore'
import styles from './Landing.module.css'

export function Landing() {
  const { loadDemo } = useAppStore()

  return (
    <div className={styles.landing}>
      <div className={styles.icon}>🍦</div>
      <h2 className={styles.title}>Ready when your data is</h2>
      <p className={styles.sub}>
        Drop your QIIME2 output files into the sidebar slots.
        MicroSee will instantly generate all visualisations from your results.
      </p>

      <div className={styles.steps}>
        <div className={styles.step}>
          <span>📊</span>
          <span>feature-table.tsv</span>
        </div>
        <div className={styles.arrow}>+</div>
        <div className={styles.step}>
          <span>🔬</span>
          <span>taxonomy.tsv</span>
        </div>
        <div className={styles.arrow}>+</div>
        <div className={styles.step}>
          <span>📋</span>
          <span>metadata.tsv</span>
        </div>
        <div className={styles.arrow}>→</div>
        <div className={`${styles.step} ${styles.stepResult}`}>
          <span>📈</span>
          <span>Visualisations</span>
        </div>
      </div>

      <div className={styles.or}>— or —</div>

      <button className={styles.demoBtn} onClick={loadDemo}>
        ✦ Explore Demo Data
      </button>

      <p className={styles.note}>
        Demo: EAA vs Whey protein · 12 weeks · older adults · 16S rRNA (PMID 38426663)
      </p>
    </div>
  )
}
