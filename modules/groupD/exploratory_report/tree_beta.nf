process BETA_DIVERSITY_TREE {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    container "containers/groupD.sif"

    input:
    
    tuple path (metadata), path (bray)
    

    output:
    path "beta_diversity_tree", emit: tree_png

    script:
    """
    python tree_beta.py ${bray} ${metadata} "beta_diversity_tree"
    """
}