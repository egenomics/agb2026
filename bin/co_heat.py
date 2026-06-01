#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load data from Nextflow arguments
metadata_file = sys.argv[1]
output_dir = sys.argv[2]

df = pd.read_csv(metadata_file, sep="\t")

# 3. Select columns that represent clinical conditions
disease_cols = [
    'acid_reflux', 'autoimmune', 'cardiovascular_disease', 'diabetes', 
    'ibd', 'ibs', 'kidney_disease', 'liver_disease', 'lung_disease', 
    'mental_illness', 'migraine', 'skin_condition', 'thyroid', 'sibo'
]

disease_cols = [col for col in disease_cols if col in df.columns]

# 4. Clean and Binarize the Data
df_clean = pd.DataFrame()
for col in disease_cols:
    df_clean[col] = df[col].astype(str).str.contains('Diagnosed|Yes', case=False, na=False).astype(int)

df_clean = df_clean.loc[:, df_clean.sum() > 0]

# 5. Calculate Co-occurrence
co_occurrence = df_clean.T.dot(df_clean)
np.fill_diagonal(co_occurrence.values, 0)

# 6. Build the Heatmap Plot
plt.figure(figsize=(10, 8))
sns.heatmap(co_occurrence, annot=True, fmt="d", cmap="YlOrRd", linewidths=.5, cbar_kws={'label': 'Number of Co-occurring Patients'})

plt.title('Clinical Co-occurrence Heatmap\n(How conditions overlap in the cohort)', fontsize=14, pad=15)
plt.xlabel('Condition A')
plt.ylabel('Condition B')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save using the Nextflow argument (Force PNG format)
os.makedirs(output_dir, exist_ok=True)

plt.savefig(os.path.join(output_dir, "clinical_association_map.png"), dpi=150)