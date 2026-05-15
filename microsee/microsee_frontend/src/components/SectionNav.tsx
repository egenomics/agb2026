/**
 * components/SectionNav.tsx
 * Section navigation tabs — same sections as MicroSee.html.
 */

import { useState } from 'react'
import styles from './SectionNav.module.css'

const SECTIONS = [
  { key: 'all',        label: '✦ All'         },
  { key: 'taxonomy',   label: '🔬 Taxonomy'    },
  { key: 'alpha',      label: '📊 Alpha'       },
  { key: 'beta',       label: '🌐 Beta'        },
  { key: 'individual', label: '👤 Individual'  },
  { key: 'compare',    label: '⚖️ Comparative' },
  { key: 'clinical',   label: '💊 Clinical'    },
  { key: 'longitudinal',label: '📈 Longitudinal'},
  { key: 'stats',      label: '📋 Stats'       },
]

interface Props {
  activeSection: string
  onChange: (key: string) => void
}

export function SectionNav({ activeSection, onChange }: Props) {
  return (
    <nav className={styles.nav}>
      {SECTIONS.map((s) => (
        <button
          key={s.key}
          className={`${styles.btn} ${activeSection === s.key ? styles.active : ''}`}
          onClick={() => onChange(s.key)}
        >
          {s.label}
        </button>
      ))}
    </nav>
  )
}

export function useSectionNav() {
  const [activeSection, setActiveSection] = useState('all')
  const isVisible = (key: string) => activeSection === 'all' || activeSection === key
  return { activeSection, setActiveSection, isVisible }
}
