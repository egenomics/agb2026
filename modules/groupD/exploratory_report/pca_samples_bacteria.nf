process PCA_MICROBIOME {

    tag "${metadata.baseName}"

    publishDir "${params.outdir}/pca",
        mode: 'copy'

    // Load the needed containers
    container "containers/r_groupD.sif"

    input:
    path(asv_table)
    path(taxonomy)
    path(metadata)

    output:
    path("pca_results")
    path("pca_results/pca_samples_bacteria.rds"),
        emit: pca_object
    path("pca_results/genus_counts.rds"),
        emit: pca_counts

    script:
    """
    mkdir -p pca_results

    Rscript ${moduleDir}/PCA_samples_bacteria.R \
        --asv ${asv_table} \
        --taxonomy ${taxonomy} \
        --metadata ${metadata} \
        --outdir pca_results
    """
}