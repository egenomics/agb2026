process DADA2_EXPORT {
    label 'process_low'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(seqtab)

    output:
    path "rep_seqs.fasta", emit: fasta
    path "asv_table.tsv",  emit: table

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    seqtab <- readRDS("$seqtab")
    seqs <- colnames(seqtab)
    ids  <- paste0("ASV", seq_along(seqs))

    # Write FASTA of representative sequences
    fasta_lines <- c(rbind(paste0(">", ids), seqs))
    writeLines(fasta_lines, "rep_seqs.fasta")

    # Write TSV abundance table (PICRUSt2 format: ASV rows, sample columns)
    abund <- t(seqtab)
    rownames(abund) <- ids
    df <- data.frame(cbind(ids, abund), check.names = FALSE)
    colnames(df)[1] <- "#OTU ID"
    write.table(df, "asv_table.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)
    """
}
