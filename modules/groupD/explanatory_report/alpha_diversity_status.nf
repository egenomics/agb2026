process alpha_diversity_status{

    publishDir "${params.outdir}/explanatory", mode: 'copy'

    label "groupD"

    input:
    tuple path(metadata),
          path(faith, stageAs: "faith.tsv"),
          path(observed, stageAs: "observed.tsv"),
          path(shannon, stageAs: "shannon.tsv"),
          path(simpson, stageAs: "simpson.tsv")

    output:
    path "alpha_patient", emit: alpha_div_dist_dir  
    path "alpha_zscore_status.rds", emit: zscore_status

    script:
    """
    alpha_diversity_status.sh $metadata $faith $observed $shannon $simpson "alpha_patient"
    """
}