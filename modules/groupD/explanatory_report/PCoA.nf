process PCoA_plots{

    publishDir "${params.outdir}/explanatory", mode: 'copy'
    
    container "containers/groupD.sif"

    input:
    tuple path(pca), path(metadata)

    output:
    path "explanatory/PCoA_patient", emit: PCoA_patient_plots_dir  

    script:
    """
    PCoA.sh $pca $metadata "explanatory/PCoA_patient"
    """
}