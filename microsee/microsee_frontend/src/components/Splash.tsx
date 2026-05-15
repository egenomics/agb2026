/**
 * components/Splash.tsx
 *
 * Ice cream splash screen — exact same animation as MicroSee.html,
 * rebuilt as a React component with CSS modules.
 *
 * Dismisses itself after 6.5 s OR as soon as data is ready.
 */

import { useEffect, useRef } from 'react'
import { useAppStore } from '@/store/appStore'
import styles from './Splash.module.css'

export function Splash() {
  const { setScreen, screen } = useAppStore()
  const timerRef = useRef<ReturnType<typeof setTimeout> | null>(null)

  useEffect(() => {
    if (screen !== 'splash') return
    timerRef.current = setTimeout(() => setScreen('landing'), 6500)
    return () => { if (timerRef.current) clearTimeout(timerRef.current) }
  }, [screen, setScreen])

  if (screen !== 'splash') return null

  return (
    <div className={styles.splash}>
      <div className={styles.iceWrap}>
        <svg
          className={styles.iceFloat}
          viewBox="0 0 100 175"
          xmlns="http://www.w3.org/2000/svg"
          aria-hidden="true"
        >
          {/* Cone */}
          <g className={styles.coneG}>
            <polygon points="50,165 18,95 82,95" fill="#C4803A" />
            <line x1="50" y1="165" x2="34" y2="95" stroke="#A86C28" strokeWidth="1.2" opacity=".6" />
            <line x1="50" y1="165" x2="50" y2="95" stroke="#A86C28" strokeWidth="1.2" opacity=".6" />
            <line x1="50" y1="165" x2="66" y2="95" stroke="#A86C28" strokeWidth="1.2" opacity=".6" />
            <line x1="23" y1="115" x2="77" y2="115" stroke="#A86C28" strokeWidth="1" opacity=".5" />
            <line x1="27" y1="130" x2="73" y2="130" stroke="#A86C28" strokeWidth="1" opacity=".5" />
          </g>
          {/* Ice cream layers — each pops in with a staggered animation */}
          <ellipse className={styles.l6} cx="50" cy="87" rx="34" ry="14" fill="#5A2614" />
          <ellipse className={styles.l5} cx="50" cy="73" rx="28" ry="12" fill="#4A1E10" />
          <ellipse className={styles.l4} cx="50" cy="60" rx="22" ry="10" fill="#D84E6A" />
          <ellipse className={styles.l3} cx="50" cy="49" rx="16" ry="8"  fill="#4A7ED4" />
          <ellipse className={styles.l2} cx="50" cy="40" rx="10" ry="6"  fill="#2FA896" />
          <ellipse className={styles.l1} cx="50" cy="33" rx="6"  ry="4"  fill="#D97A3A" />
          <ellipse className={styles.l0} cx="50" cy="29" rx="3.5" ry="2.5" fill="#9058C4" />
          {/* Drip */}
          <path
            className={styles.dripP}
            d="M 77,82 Q 82,96 79,110 Q 77,120 80,128"
            stroke="#3A1408"
            strokeWidth="3"
            fill="none"
            strokeLinecap="round"
          />
        </svg>
      </div>

      <h1 className={styles.logo}>MicroSee</h1>
      <p className={styles.sub}>
        Your microbiome visualisation suite &middot; AGB 2026 &middot; UPF
      </p>
      <button
        className={styles.btn}
        onClick={() => setScreen('landing')}
      >
        Start exploring →
      </button>
    </div>
  )
}
