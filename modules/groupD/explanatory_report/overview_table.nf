process overview_table{

    publishDir "${params.outdir}/explanatory", mode: 'copy'
    
    container "containers/groupD.sif"

    input:
    tuple path(metadata), path(pca), path(z_scores), path(genus_counts)

    output:
    path "explanatory/overview_sample", emit: overview_sample_table_dir  

    script:
    """
    overview_table.sh $metadata $pca $z_scores $genus_counts "explanatory/overview_sample"
    """
}