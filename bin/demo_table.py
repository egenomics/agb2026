import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# 1. Load data from Nextflow arguments
metadata_file = sys.argv[1]
output_dir = sys.argv[2]

df = pd.read_csv(metadata_file, sep="\t")

# 3. Define the demographic columns we want to summarize
demo_cols = ['age_cat', 'sex', 'bmi_cat', 'diet_type', 'exercise_frequency', 'antibiotic_history']

# Clean up the 'healthy' column for clear labeling
df['healthy_status'] = df['healthy'].map({'yes': 'Healthy', 'no': 'Has Condition'})

# 4. Function to build the Table data
table_rows = []
total_healthy = len(df[df['healthy_status'] == 'Healthy'])
total_condition = len(df[df['healthy_status'] == 'Has Condition'])
total_all = len(df)

# Loop through each demographic variable
for col in demo_cols:
    table_rows.append({
        "Characteristic": f"{col.replace('_', ' ').upper()}", 
        "Overall (N=81)": "", f"Healthy (N={total_healthy})": "", f"Has Condition (N={total_condition})": ""
    })
    
    categories = df[col].value_counts(dropna=False).index
    
    for cat in categories:
        cat_name = "Not Provided" if pd.isna(cat) or cat == 'not provided' else str(cat)
        
        overall_count = len(df[df[col] == cat]) if pd.notna(cat) else df[col].isna().sum()
        overall_perc = (overall_count / total_all) * 100
        overall_str = f"{overall_count} ({overall_perc:.1f}%)"
        
        h_df = df[df['healthy_status'] == 'Healthy']
        h_count = len(h_df[h_df[col] == cat]) if pd.notna(cat) else h_df[col].isna().sum()
        h_perc = (h_count / total_healthy * 100) if total_healthy > 0 else 0
        h_str = f"{h_count} ({h_perc:.1f}%)"
        
        c_df = df[df['healthy_status'] == 'Has Condition']
        c_count = len(c_df[c_df[col] == cat]) if pd.notna(cat) else c_df[col].isna().sum()
        c_perc = (c_count / total_condition * 100) if total_condition > 0 else 0
        c_str = f"{c_count} ({c_perc:.1f}%)"
        
        table_rows.append({
            "Characteristic": f"    {cat_name}",
            "Overall (N=81)": overall_str,
            f"Healthy (N={total_healthy})": h_str,
            f"Has Condition (N={total_condition})": c_str
        })

df_table1 = pd.DataFrame(table_rows)

# 5. Plot as a clean paper-style table using Matplotlib
fig, ax = plt.subplots(figsize=(14, 10))
ax.axis('off') 

col_widths = [0.45, 0.18, 0.18, 0.19]
table = ax.table(cellText=df_table1.values, colLabels=df_table1.columns, 
                 loc='center', cellLoc='left', colWidths=col_widths, bbox=[0, 0, 1, 0.93]) 

table.auto_set_font_size(False)
table.set_fontsize(10)

# 6. Apply academic styling
for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor('#d3d3d3') 
    if row == 0:
        cell.set_text_props(weight='bold', color='white', size=11)
        cell.set_facecolor('#2b579a') 
    else:
        text_val = str(df_table1.iloc[row-1, 0])
        if not text_val.startswith("    "):
            cell.set_text_props(weight='bold')
            cell.set_facecolor('#e8e8e8') 
        else:
            cell.set_facecolor('#ffffff' if row % 2 == 0 else '#f9f9f9')

plt.title("Table: Cohort Baseline Characteristics", fontsize=16, weight='bold', y=0.96)

# Save using the Nextflow argument (Force PNG format)
os.makedirs(output_dir, exist_ok=True)

plt.savefig(os.path.join(output_dir, "demographic_table.png"), dpi=300, bbox_inches='tight', format='png')