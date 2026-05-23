process get_sampling_depth {

    input:
    path table

    output:
    path "sampling_depth.txt", emit: depth

    script:
    """
    qiime tools export \\
        --input-path ${table} \\
        --output-path exported_table

    biom summarize-table \\
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

        match = re.match(r'^\\S+:\\s+([\\d.,]+)', line)
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
    """
}

process diversity_analysis {
    tag "alpha_beta_diversity"

    input:
    path phylogeny
    path table
    path metadata
    val  sampling_depth

    output:
    path "diversity_table/shannon/*",          emit: shannon
    path "diversity_table/observed/*",         emit: observed_features
    path "diversity_table/faith/*",            emit: faith_pd
    path "diversity_table/simpson/*",          emit: simpson
    path "diversity_table/weighted_unifrac/*", emit: weighted_unifrac
    path "diversity_table/bray_curtis/*",      emit: bray_curtis
    path "core_metrics_output/*.qza",          emit: qza_files
    path "core_metrics_output/*.qzv",          emit: qzv_files

    script:
    """
    rm -rf core_metrics_output
    mkdir -p core_metrics_output \\
             diversity_table/shannon \\
             diversity_table/observed \\
             diversity_table/faith \\
             diversity_table/simpson \\
             diversity_table/weighted_unifrac \\
             diversity_table/bray_curtis

    qiime diversity core-metrics-phylogenetic \\
        --i-phylogeny ${phylogeny} \\
        --i-table ${table} \\
        --p-sampling-depth ${sampling_depth} \\
        --m-metadata-file ${metadata} \\
        --o-rarefied-table                     core_metrics_output/rarefied_table.qza \\
        --o-faith-pd-vector                    core_metrics_output/faith_pd_vector.qza \\
        --o-observed-features-vector           core_metrics_output/observed_features_vector.qza \\
        --o-shannon-vector                     core_metrics_output/shannon_vector.qza \\
        --o-evenness-vector                    core_metrics_output/evenness_vector.qza \\
        --o-unweighted-unifrac-distance-matrix core_metrics_output/unweighted_unifrac_distance_matrix.qza \\
        --o-weighted-unifrac-distance-matrix   core_metrics_output/weighted_unifrac_distance_matrix.qza \\
        --o-jaccard-distance-matrix            core_metrics_output/jaccard_distance_matrix.qza \\
        --o-bray-curtis-distance-matrix        core_metrics_output/bray_curtis_distance_matrix.qza \\
        --o-unweighted-unifrac-pcoa-results    core_metrics_output/unweighted_unifrac_pcoa_results.qza \\
        --o-weighted-unifrac-pcoa-results      core_metrics_output/weighted_unifrac_pcoa_results.qza \\
        --o-jaccard-pcoa-results               core_metrics_output/jaccard_pcoa_results.qza \\
        --o-bray-curtis-pcoa-results           core_metrics_output/bray_curtis_pcoa_results.qza \\
        --o-unweighted-unifrac-emperor         core_metrics_output/unweighted_unifrac_emperor.qzv \\
        --o-weighted-unifrac-emperor           core_metrics_output/weighted_unifrac_emperor.qzv \\
        --o-jaccard-emperor                    core_metrics_output/jaccard_emperor.qzv \\
        --o-bray-curtis-emperor                core_metrics_output/bray_curtis_emperor.qzv

    qiime diversity alpha \\
        --i-table ${table} \\
        --p-metric simpson \\
        --o-alpha-diversity core_metrics_output/simpson_vector.qza

    qiime tools export \\
        --input-path core_metrics_output/shannon_vector.qza \\
        --output-path diversity_table/shannon

    qiime tools export \\
        --input-path core_metrics_output/observed_features_vector.qza \\
        --output-path diversity_table/observed

    qiime tools export \\
        --input-path core_metrics_output/faith_pd_vector.qza \\
        --output-path diversity_table/faith

    qiime tools export \\
        --input-path core_metrics_output/simpson_vector.qza \\
        --output-path diversity_table/simpson

    sed -i 's/functools.partial(<function _simpsons_dominance.*>)/simpson/' \\
        diversity_table/simpson/alpha-diversity.tsv

    qiime tools export \\
        --input-path core_metrics_output/weighted_unifrac_distance_matrix.qza \\
        --output-path diversity_table/weighted_unifrac

    qiime tools export \\
        --input-path core_metrics_output/bray_curtis_distance_matrix.qza \\
        --output-path diversity_table/bray_curtis
    """
}

workflow {
    phylogeny = file("${params.data_dir}/rooted-tree.qza")
    table     = file("${params.data_dir}/table.qza")
    metadata  = file("${params.data_dir}/sample-metadata.tsv")

    get_sampling_depth(table)

    depth = get_sampling_depth.out.depth
                .map { it.text.trim().toInteger() }

    diversity_analysis(phylogeny, table, metadata, depth)
}
