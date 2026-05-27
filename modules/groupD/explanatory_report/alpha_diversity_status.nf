process alpha_diversity_status{

    publishDir "${params.outdir}", mode: 'copy'
    conda 'conda-forge::r-tidyverse=2.0.0 conda-forge::r-optparse=1.8.2'
    
    input:
    tuple path(metadata), path(bray), path(unifrac), path(faith), path(observed), path(shannon), path(simpson)

    output:
    path "alpha_patient", emit: alpha_div_dist_dir  

    script:
    """
    alpha_diversity_status.sh $metadata $bray $unifrac $faith $observed $shannon $simpson "alpha_patient"
    """
}