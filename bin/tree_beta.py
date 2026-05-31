import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# 1. Load data from Nextflow arguments
dm_file = sys.argv[1]
metadata_file = sys.argv[2]
output_img = sys.argv[3]

dm_df = pd.read_csv(dm_file, sep="\t", index_col=0)
meta_df = pd.read_csv(metadata_file, sep="\t")
status_map = dict(zip(meta_df['sra_id'], meta_df['healthy']))

# 3. CLEAN AND CONDENSE MATRIX
matrix = dm_df.values
matrix = (matrix + matrix.T) / 2
np.fill_diagonal(matrix, 0)
condensed_dm = squareform(matrix)

# 4. PERFORM HIERARCHICAL CLUSTERING (UPGMA)
Z = linkage(condensed_dm, method='average')

# 5. PLOT DENDROGRAM TREE
fig, ax = plt.subplots(figsize=(13, 7))
dend_results = dendrogram(Z, labels=dm_df.index, leaf_rotation=90, leaf_font_size=10, color_threshold=0.45, ax=ax)

for lbl in ax.get_xticklabels():
    sample_id = lbl.get_text()
    if status_map.get(sample_id) == "yes":
        lbl.set_color("green")
    else:
        lbl.set_color("red")

# 6. ADD LEGEND FOR SAMPLES
green_patch = mpatches.Patch(color='green', label='Controls')
red_patch = mpatches.Patch(color='red', label='Cases')
ax.legend(handles=[green_patch, red_patch], loc='upper left', frameon=False, fontsize=11)

# 7. APPLY ACADEMIC STYLING
ax.set_title("Beta Diversity Hierarchical Clustering (UPGMA Tree)", fontsize=16, weight='bold', pad=15)
ax.set_ylabel("Distance", fontsize=12)
ax.grid(axis='y', linestyle='--', alpha=0.5)

for spine in ['top', 'right', 'left', 'bottom']:
    ax.spines[spine].set_visible(False)

# 8. SAVE OUTPUT IMAGE (Force PNG format)
plt.savefig(output_img, dpi=300, bbox_inches='tight', format='png')