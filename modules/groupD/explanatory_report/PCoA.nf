process PCoA_plots{

    publishDir "${params.outdir}/explanatory", mode: 'copy'
    
    label "groupD"

    input:
    tuple path(pca), path(metadata)

    output:
    path "PCoA_patient", emit: PCoA_patient_plots_dir  

    script:
    """
    PCoA.sh $pca $metadata "PCoA_patient"
    """
}