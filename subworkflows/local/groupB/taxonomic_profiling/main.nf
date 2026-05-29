// subworkflows/local/groupB/taxonomic_profiling/main.nf
//
// DADA2 taxonomic profiling. Map (per-sample filtntrim → err → denoising)
// then reduce (collect → merge_chimera → taxonomy).

include { DADA2_FILTNTRIM     } from '../../../../modules/groupB/dada2/filtandtrim/main'
include { DADA2_ERR           } from '../../../../modules/groupB/dada2/learnerrors/main'
include { DADA2_DENOISING     } from '../../../../modules/groupB/dada2/dada2_denoising'
include { DADA2_MERGE_CHIMERA } from '../../../../modules/groupB/dada2/dada2_merge_chimera'
include { DADA2_TAXONOMY      } from '../../../../modules/groupB/dada2/dada2_taxonomy'

workflow TAXONOMIC_PROFILING {

    take:
    ch_samplesheet            // channel: [ val(meta), [ reads ] ]
    dada2_train_set           // path:    SILVA toGenus train set (.fa)
    dada2_species_set         // path:    SILVA assignSpecies reference (.fa)

    main:
    ch_versions = channel.empty()

    // Package truncation thresholds as lists so DADA2_FILTNTRIM's internal [1]
    // indexing works for both single-end and paired-end inputs.
    ch_filt_input = ch_samplesheet.map { meta, reads ->
        tuple(meta, reads, [0, 140], [0, 0])
    }

    DADA2_FILTNTRIM(ch_filt_input)

    ch_filtered_reads = DADA2_FILTNTRIM.out.reads_logs_args
        .map { meta, reads, _stats, _args -> tuple(meta, reads) }

    DADA2_ERR(ch_filtered_reads)

    ch_denoising_input = ch_filtered_reads.join(DADA2_ERR.out.errormodel)

    DADA2_DENOISING(ch_denoising_input)

    // Map → reduce boundary: drop meta and collect all per-sample seqtabs into
    // a single emission. Without .collect(), MERGE_CHIMERA would fire N times.
    ch_seqtabs = DADA2_DENOISING.out.asv_table
        .map { _meta, rds -> rds }
        .collect()

    DADA2_MERGE_CHIMERA(ch_seqtabs)

    // Taxonomy runs ONCE on the merged seqtab (SILVA loaded once, ~90 s vs
    // ~85 s per sample).
    DADA2_TAXONOMY(
        DADA2_MERGE_CHIMERA.out.seqtab,
        dada2_train_set,
        dada2_species_set,
        "ASV_taxonomy.tsv",
        "Kingdom,Phylum,Class,Order,Family,Genus,Species"
    )

    emit:
    seqtab   = DADA2_MERGE_CHIMERA.out.seqtab   // path: seqtab_final.rds (study-wide, no meta)
    taxonomy = DADA2_TAXONOMY.out.tsv           // path: ASV_taxonomy.tsv
    versions = ch_versions
}
