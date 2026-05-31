process overvoew_table{

    publishDir "${params.outdir}/explanatory_report", mode: 'copy'
    
    input:
    tuple path(metadata), path(pca), path(z_scores), path(genus_counts)

    output:
    path "explanatory_report/overview_sample", emit: overview_sample_table_dir  

    script:
    """
    PCoA.sh $metadata $pca $z_scores $genus_counts "explanatory_report/overview_sample"
    """
}