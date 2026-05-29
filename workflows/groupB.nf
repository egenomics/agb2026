// workflows/groupB.nf
//
// Canonical Group B pipeline. Composes three subworkflows
// (QC_CHECKS, TAXONOMIC_PROFILING, FUNCTIONAL_ANNOTATION) — see
// subworkflows/local/groupB/.

include { QC_CHECKS             } from '../subworkflows/local/groupB/qc_checks/main'
include { TAXONOMIC_PROFILING   } from '../subworkflows/local/groupB/taxonomic_profiling/main'
include { FUNCTIONAL_ANNOTATION } from '../subworkflows/local/groupB/functional_annotation/main'

workflow GROUPB {

    take:
    ch_samplesheet            // channel: [ val(meta), [ reads ] ]
    multiqc_config            // path
    multiqc_logo              // path
    ch_collated_versions      // channel: [ path(versions.yml) ]
    ch_methods_description    // channel: [ path(methods_description_mqc.yaml) ]
    ch_workflow_summary       // channel: [ path(workflow_summary_mqc.yaml) ]
    dada2_train_set           // path
    dada2_species_set         // path

    main:

    // STEP 1: QUALITY CONTROL — disabled here; Group A already runs FastQC/MultiQC
    // on the same reads upstream, so re-running in Group B is redundant. Re-enable
    // (and re-sync modules/nf-core/multiqc/ to the iteration-tree version) if Group B
    // ever needs its own QC pass independent of Group A.
    // QC_CHECKS(
    //     ch_samplesheet,
    //     multiqc_config,
    //     multiqc_logo,
    //     ch_collated_versions,
    //     ch_methods_description,
    //     ch_workflow_summary
    // )

    // STEP 2: TAXONOMIC PROFILING (DADA2 map-reduce: filtntrim/err/denoising
    // per-sample, then merge+chimera+taxonomy once on the merged seqtab)
    TAXONOMIC_PROFILING(
        ch_samplesheet,
        dada2_train_set,
        dada2_species_set
    )

    // STEP 3: FUNCTIONAL ANNOTATION (DADA2_EXPORT + PICRUSt2)
    FUNCTIONAL_ANNOTATION(
        TAXONOMIC_PROFILING.out.seqtab
    )

    emit:
    // DELIVERABLES FOR GROUP C — all study-wide, plain `path` channels (no meta)
    table_counts   = FUNCTIONAL_ANNOTATION.out.table_counts   // asv_table.tsv
    rep_seqs       = FUNCTIONAL_ANNOTATION.out.rep_seqs       // rep_seqs.fasta
    taxonomy       = TAXONOMIC_PROFILING.out.taxonomy         // ASV_taxonomy.tsv

    // GENERAL OUTPUTS
    // multiqc_report = QC_CHECKS.out.multiqc_report  // re-enable with QC_CHECKS call above
    pathways       = FUNCTIONAL_ANNOTATION.out.pathways
    versions       = channel.empty()
}
