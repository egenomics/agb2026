// subworkflows/local/groupB/functional_annotation/main.nf
//
// Convert the merged DADA2 seqtab to rep_seqs.fasta + asv_table.tsv, then
// run PICRUSt2 on those study-wide inputs. Both steps run exactly once.

include { DADA2_EXPORT } from '../../../../modules/groupB/dada2/dada2_export'
include { PICRUST      } from '../../../../modules/groupB/picrust'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_seqtab                 // path: seqtab_final.rds (study-wide, no meta)

    main:
    ch_versions = channel.empty()

    DADA2_EXPORT(ch_seqtab)

    PICRUST(
        DADA2_EXPORT.out.fasta,
        DADA2_EXPORT.out.table,
        "metagenome",
        true
    )

    emit:
    table_counts = DADA2_EXPORT.out.table      // path: asv_table.tsv
    rep_seqs     = DADA2_EXPORT.out.fasta      // path: rep_seqs.fasta
    pathways     = PICRUST.out.pathways
    versions     = ch_versions
}
