process PCoA_plots{

    publishDir "${params.outdir}/explanatory_report", mode: 'copy'
    
    input:
    tuple path(pca), path(metadata)

    output:
    path "explanatory_report/PCoA_patient", emit: PCoA_patient_plots_dir  

    script:
    """
    PCoA.sh $pca $metadata "explanatory_report/PCoA_patient"
    """
}