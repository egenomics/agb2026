/**
 * data/demo.ts
 *
 * Typed demo dataset — EAA vs Whey protein supplementation in older adults.
 * PMID: 38426663 | Mol. Nutr. Food Res. 2024
 *
 * Same 24 samples as MicroSee.html, now typed as SampleRow[].
 * This path never hits the backend — pure frontend, works offline.
 */

import type { SampleRow, IntegrateResult } from '@/types/sample'
import { TAXA } from '@/types/sample'

export const DEMO_ROWS: SampleRow[] = [
  // ── EAA group ──────────────────────────────────────────────────────────────
  { sample_id:'EAA01_T0',  group:'EAA_T0',  base_group:'EAA', patient:'EAA01', timepoint:'T0',  time:0,  Bacteroidaceae:36, Lachnospiraceae:22, Ruminococcaceae:14, Prevotellaceae:9,  Rikenellaceae:6,  Enterobacteriaceae:3, Oscillospiraceae:4, Tannerellaceae:1, Akkermansiaceae:5, shannon:1.791, simpson:0.786, sixmwt:310, il18:248 },
  { sample_id:'EAA01_T84', group:'EAA_T84', base_group:'EAA', patient:'EAA01', timepoint:'T84', time:84, Bacteroidaceae:34, Lachnospiraceae:23, Ruminococcaceae:14, Prevotellaceae:10, Rikenellaceae:7,  Enterobacteriaceae:3, Oscillospiraceae:4, Tannerellaceae:2, Akkermansiaceae:3, shannon:1.814, simpson:0.793, sixmwt:372, il18:192 },
  { sample_id:'EAA02_T0',  group:'EAA_T0',  base_group:'EAA', patient:'EAA02', timepoint:'T0',  time:0,  Bacteroidaceae:28, Lachnospiraceae:24, Ruminococcaceae:16, Prevotellaceae:15, Rikenellaceae:5,  Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:3, Akkermansiaceae:1, shannon:1.835, simpson:0.809, sixmwt:295, il18:262 },
  { sample_id:'EAA02_T84', group:'EAA_T84', base_group:'EAA', patient:'EAA02', timepoint:'T84', time:84, Bacteroidaceae:27, Lachnospiraceae:25, Ruminococcaceae:17, Prevotellaceae:14, Rikenellaceae:5,  Enterobacteriaceae:5, Oscillospiraceae:4, Tannerellaceae:2, Akkermansiaceae:1, shannon:1.829, simpson:0.809, sixmwt:358, il18:208 },
  { sample_id:'EAA03_T0',  group:'EAA_T0',  base_group:'EAA', patient:'EAA03', timepoint:'T0',  time:0,  Bacteroidaceae:25, Lachnospiraceae:30, Ruminococcaceae:16, Prevotellaceae:10, Rikenellaceae:6,  Enterobacteriaceae:3, Oscillospiraceae:5, Tannerellaceae:3, Akkermansiaceae:2, shannon:1.838, simpson:0.804, sixmwt:340, il18:235 },
  { sample_id:'EAA03_T84', group:'EAA_T84', base_group:'EAA', patient:'EAA03', timepoint:'T84', time:84, Bacteroidaceae:26, Lachnospiraceae:29, Ruminococcaceae:16, Prevotellaceae:11, Rikenellaceae:6,  Enterobacteriaceae:3, Oscillospiraceae:5, Tannerellaceae:2, Akkermansiaceae:2, shannon:1.826, simpson:0.803, sixmwt:395, il18:182 },
  { sample_id:'EAA04_T0',  group:'EAA_T0',  base_group:'EAA', patient:'EAA04', timepoint:'T0',  time:0,  Bacteroidaceae:22, Lachnospiraceae:24, Ruminococcaceae:22, Prevotellaceae:11, Rikenellaceae:7,  Enterobacteriaceae:4, Oscillospiraceae:5, Tannerellaceae:3, Akkermansiaceae:2, shannon:1.900, simpson:0.823, sixmwt:280, il18:270 },
  { sample_id:'EAA04_T84', group:'EAA_T84', base_group:'EAA', patient:'EAA04', timepoint:'T84', time:84, Bacteroidaceae:23, Lachnospiraceae:23, Ruminococcaceae:21, Prevotellaceae:12, Rikenellaceae:7,  Enterobacteriaceae:4, Oscillospiraceae:5, Tannerellaceae:3, Akkermansiaceae:2, shannon:1.906, simpson:0.825, sixmwt:338, il18:215 },
  { sample_id:'EAA05_T0',  group:'EAA_T0',  base_group:'EAA', patient:'EAA05', timepoint:'T0',  time:0,  Bacteroidaceae:30, Lachnospiraceae:22, Ruminococcaceae:15, Prevotellaceae:8,  Rikenellaceae:10, Enterobacteriaceae:3, Oscillospiraceae:5, Tannerellaceae:5, Akkermansiaceae:2, shannon:1.894, simpson:0.816, sixmwt:325, il18:255 },
  { sample_id:'EAA05_T84', group:'EAA_T84', base_group:'EAA', patient:'EAA05', timepoint:'T84', time:84, Bacteroidaceae:29, Lachnospiraceae:23, Ruminococcaceae:15, Prevotellaceae:8,  Rikenellaceae:10, Enterobacteriaceae:3, Oscillospiraceae:5, Tannerellaceae:5, Akkermansiaceae:2, shannon:1.897, simpson:0.818, sixmwt:385, il18:200 },
  { sample_id:'EAA06_T0',  group:'EAA_T0',  base_group:'EAA', patient:'EAA06', timepoint:'T0',  time:0,  Bacteroidaceae:28, Lachnospiraceae:26, Ruminococcaceae:16, Prevotellaceae:10, Rikenellaceae:7,  Enterobacteriaceae:4, Oscillospiraceae:5, Tannerellaceae:2, Akkermansiaceae:2, shannon:1.851, simpson:0.809, sixmwt:305, il18:244 },
  { sample_id:'EAA06_T84', group:'EAA_T84', base_group:'EAA', patient:'EAA06', timepoint:'T84', time:84, Bacteroidaceae:27, Lachnospiraceae:27, Ruminococcaceae:16, Prevotellaceae:10, Rikenellaceae:7,  Enterobacteriaceae:4, Oscillospiraceae:5, Tannerellaceae:2, Akkermansiaceae:2, shannon:1.852, simpson:0.809, sixmwt:362, il18:195 },
  // ── Whey group ─────────────────────────────────────────────────────────────
  { sample_id:'WHY01_T0',  group:'Whey_T0',  base_group:'Whey', patient:'WHY01', timepoint:'T0',  time:0,  Bacteroidaceae:42, Lachnospiraceae:20, Ruminococcaceae:13, Prevotellaceae:8,  Rikenellaceae:5, Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:3, Akkermansiaceae:1, shannon:1.712, simpson:0.754, sixmwt:315, il18:242 },
  { sample_id:'WHY01_T84', group:'Whey_T84', base_group:'Whey', patient:'WHY01', timepoint:'T84', time:84, Bacteroidaceae:43, Lachnospiraceae:19, Ruminococcaceae:13, Prevotellaceae:8,  Rikenellaceae:5, Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:3, Akkermansiaceae:1, shannon:1.704, simpson:0.749, sixmwt:330, il18:240 },
  { sample_id:'WHY02_T0',  group:'Whey_T0',  base_group:'Whey', patient:'WHY02', timepoint:'T0',  time:0,  Bacteroidaceae:23, Lachnospiraceae:32, Ruminococcaceae:17, Prevotellaceae:10, Rikenellaceae:6, Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:3, Akkermansiaceae:1, shannon:1.812, simpson:0.798, sixmwt:290, il18:258 },
  { sample_id:'WHY02_T84', group:'Whey_T84', base_group:'Whey', patient:'WHY02', timepoint:'T84', time:84, Bacteroidaceae:24, Lachnospiraceae:31, Ruminococcaceae:17, Prevotellaceae:10, Rikenellaceae:6, Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:3, Akkermansiaceae:1, shannon:1.815, simpson:0.800, sixmwt:305, il18:255 },
  { sample_id:'WHY03_T0',  group:'Whey_T0',  base_group:'Whey', patient:'WHY03', timepoint:'T0',  time:0,  Bacteroidaceae:27, Lachnospiraceae:22, Ruminococcaceae:22, Prevotellaceae:13, Rikenellaceae:5, Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:2, Akkermansiaceae:1, shannon:1.817, simpson:0.807, sixmwt:330, il18:238 },
  { sample_id:'WHY03_T84', group:'Whey_T84', base_group:'Whey', patient:'WHY03', timepoint:'T84', time:84, Bacteroidaceae:28, Lachnospiraceae:21, Ruminococcaceae:22, Prevotellaceae:13, Rikenellaceae:5, Enterobacteriaceae:4, Oscillospiraceae:4, Tannerellaceae:2, Akkermansiaceae:1, shannon:1.814, simpson:0.806, sixmwt:345, il18:242 },
  { sample_id:'WHY04_T0',  group:'Whey_T0',  base_group:'Whey', patient:'WHY04', timepoint:'T0',  time:0,  Bacteroidaceae:30, Lachnospiraceae:23, Ruminococcaceae:15, Prevotellaceae:10, Rikenellaceae:5, Enterobacteriaceae:7, Oscillospiraceae:5, Tannerellaceae:3, Akkermansiaceae:2, shannon:1.883, simpson:0.813, sixmwt:285, il18:265 },
  { sample_id:'WHY04_T84', group:'Whey_T84', base_group:'Whey', patient:'WHY04', timepoint:'T84', time:84, Bacteroidaceae:31, Lachnospiraceae:22, Ruminococcaceae:15, Prevotellaceae:10, Rikenellaceae:5, Enterobacteriaceae:7, Oscillospiraceae:5, Tannerellaceae:3, Akkermansiaceae:2, shannon:1.880, simpson:0.812, sixmwt:298, il18:262 },
  { sample_id:'WHY05_T0',  group:'Whey_T0',  base_group:'Whey', patient:'WHY05', timepoint:'T0',  time:0,  Bacteroidaceae:32, Lachnospiraceae:21, Ruminococcaceae:14, Prevotellaceae:9,  Rikenellaceae:6, Enterobacteriaceae:3, Oscillospiraceae:8, Tannerellaceae:5, Akkermansiaceae:2, shannon:1.888, simpson:0.812, sixmwt:320, il18:248 },
  { sample_id:'WHY05_T84', group:'Whey_T84', base_group:'Whey', patient:'WHY05', timepoint:'T84', time:84, Bacteroidaceae:33, Lachnospiraceae:20, Ruminococcaceae:14, Prevotellaceae:9,  Rikenellaceae:6, Enterobacteriaceae:3, Oscillospiraceae:8, Tannerellaceae:5, Akkermansiaceae:2, shannon:1.884, simpson:0.810, sixmwt:335, il18:250 },
  { sample_id:'WHY06_T0',  group:'Whey_T0',  base_group:'Whey', patient:'WHY06', timepoint:'T0',  time:0,  Bacteroidaceae:27, Lachnospiraceae:25, Ruminococcaceae:17, Prevotellaceae:11, Rikenellaceae:7, Enterobacteriaceae:4, Oscillospiraceae:5, Tannerellaceae:2, Akkermansiaceae:2, shannon:1.865, simpson:0.814, sixmwt:300, il18:255 },
  { sample_id:'WHY06_T84', group:'Whey_T84', base_group:'Whey', patient:'WHY06', timepoint:'T84', time:84, Bacteroidaceae:26, Lachnospiraceae:26, Ruminococcaceae:17, Prevotellaceae:11, Rikenellaceae:7, Enterobacteriaceae:4, Oscillospiraceae:5, Tannerellaceae:2, Akkermansiaceae:2, shannon:1.866, simpson:0.814, sixmwt:315, il18:253 },
]

export const DEMO_RESULT: IntegrateResult = {
  rows:         DEMO_ROWS,
  taxa:         [...TAXA],
  n_samples:    24,
  n_taxa:       9,
  groups:       ['EAA_T0', 'EAA_T84', 'Whey_T0', 'Whey_T84'],
  has_clinical: true,
  warnings:     [],
}
