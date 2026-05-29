process DADA2_MERGE_CHIMERA {
    tag "merge_chimera"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    path(seqtabs, stageAs: 'seqtabs/*')

    output:
    path "seqtab_final.rds"            , emit: seqtab
    path "merge_chimera.log"           , emit: log
    path "removeBimeraDenovo.args.txt" , emit: args
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'method = "consensus"'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    sink(file = "merge_chimera.log")

    seqtab_files <- sort(list.files("seqtabs", pattern = "\\\\.seqtab\\\\.rds\$", full.names = TRUE))
    cat("Per-sample seqtab files found:", length(seqtab_files), "\\n")
    stopifnot(length(seqtab_files) > 0)

    # DADA2_DENOISING runs one sample at a time and calls makeSequenceTable() on a
    # single dada object, which produces a 1-row matrix with NULL rownames.
    # mergeSequenceTables rejects tables without sample names ("invalid table"),
    # so assign the sample ID — derived from the <sample_id>.seqtab.rds filename
    # set by ext.prefix in modules.config — before merging.
    seqtabs <- lapply(seqtab_files, function(f) {
        st <- readRDS(f)
        sample_id <- sub("\\\\.seqtab\\\\.rds\$", "", basename(f))
        rownames(st) <- sample_id
        st
    })

    # Fan-in: combine all per-sample 1-row matrices into one study-wide matrix.
    # 'repeats="sum"' guards against the rare case of duplicate row names.
    if (length(seqtabs) == 1) {
        merged <- seqtabs[[1]]
    } else {
        merged <- do.call(mergeSequenceTables, c(seqtabs, list(repeats = "sum")))
    }
    cat("Merged seqtab dim (samples x ASVs):", dim(merged), "\\n")

    # Global chimera removal (needs full abundance distribution to work)
    seqtab_nochim <- removeBimeraDenovo(merged, $args, multithread = $task.cpus, verbose = TRUE)
    cat("Post-chimera dim (samples x ASVs):", dim(seqtab_nochim), "\\n")
    cat("ASVs kept fraction:", round(sum(seqtab_nochim) / sum(merged), 4), "\\n")

    saveRDS(seqtab_nochim, "seqtab_final.rds")

    sink(file = NULL)

    write.table('removeBimeraDenovo\\t$args', file = "removeBimeraDenovo.args.txt",
                row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":",
                 paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                 paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}
