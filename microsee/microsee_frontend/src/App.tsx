/**
 * App.tsx
 * Root component. Manages splash → landing → results transitions.
 */

import { useAppStore }   from '@/store/appStore'
import { Splash }        from '@/components/Splash'
import { Sidebar }       from '@/components/Sidebar'
import { Landing }       from '@/components/Landing'
import { ResultsPage }   from '@/pages/ResultsPage'
import styles from './App.module.css'

export function App() {
  const { screen } = useAppStore()

  return (
    <>
      {/* Splash — renders on top of everything */}
      <Splash />

      {/* Main app layout — fades in after splash */}
      <div className={`${styles.app} ${screen !== 'splash' ? styles.visible : ''}`}>
        <Sidebar />
        <main className={styles.main}>
          {screen === 'landing' && <Landing />}
          {screen === 'results' && <ResultsPage />}
        </main>
      </div>
    </>
  )
}
