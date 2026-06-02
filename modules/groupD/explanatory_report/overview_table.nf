process overview_table{

    publishDir "${params.outdir}/explanatory", mode: 'copy'
    
    label "groupD"

    input:
    tuple path(metadata), path(pca), path(z_scores), path(genus_counts)

    output:
    path "overview_sample", emit: overview_sample_table_dir  

    script:
    """
    overview_table.sh $metadata $pca $z_scores $genus_counts "overview_sample"
    """
}