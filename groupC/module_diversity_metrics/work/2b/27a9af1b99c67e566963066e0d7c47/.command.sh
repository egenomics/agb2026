#!/bin/bash -ue
qiime tools export \
        --input-path table.qza \
        --output-path exported_table

    biom summarize-table \
        -i exported_table/feature-table.biom > summary.txt

    python3 << 'EOF'
import re

counts = []
in_detail = False

with open('summary.txt', 'r') as f:
    for line in f:
        line = line.strip()

        if 'Counts/sample detail:' in line:
            in_detail = True
            continue

        if not in_detail:
            continue

        match = re.match(r'^\S+:\s+([\d.,]+)', line)
        if match:
            raw = match.group(1)
            raw = raw.replace('.', '')
            raw = raw.replace(',', '.')
            val = int(float(raw))
            counts.append(val)

if not counts:
    raise ValueError("No sample counts found")

counts.sort()
q1 = counts[len(counts) // 4]

filtered = [c for c in counts if c >= q1]

depth = int(min(filtered) * 0.9)

print(f"All counts: {counts}")
print(f"Q1 threshold: {q1}")
print(f"Filtered counts: {filtered}")
print(f"Sampling depth (90% of min filtered): {depth}")

with open("sampling_depth.txt", "w") as f:
    f.write(str(depth))
EOF
