process alpha_diversity_status{

    publishDir "${params.outdir}/explanatory_report", mode: 'copy'

    input:
    tuple path(metadata), path(faith), path(observed), path(shannon), path(simpson)

    output:
    path "explanatory_report/alpha_patient", emit: alpha_div_dist_dir  
    path "alpha_zscore_status.rds", emit: zscore_status

    script:
    """
    alpha_diversity_status.sh $metadata $faith $observed $shannon $simpson "explanatory_report/alpha_patient"
    """
}