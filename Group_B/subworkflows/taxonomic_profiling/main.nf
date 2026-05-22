// subworkflows/taxonomic_profiling/main.nf

include { DADA2_DENOISING } from '../../modules/local/dada2/dada2_denoising'
include { DADA2_TAXONOMY  } from '../../modules/local/dada2/dada2_taxonomy'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'

workflow TAXONOMIC_PROFILING {
    take:
    ch_samplesheet 

    main:
    ch_versions = channel.empty()

    // Step 1: Denoise reads into ASVs
    DADA2_DENOISING(ch_samplesheet)
    ch_versions = ch_versions.mix(DADA2_DENOISING.out.versions)

    // Step 2: Assign Taxonomy using the Databases you defined
    DADA2_TAXONOMY(
        DADA2_DENOISING.out.asv_table, 
        params.dada2_train_set, 
        params.dada2_species_set
    )
    ch_versions = ch_versions.mix(DADA2_TAXONOMY.out.versions)

    // Step 3: Run Kraken2 for a second opinion
    KRAKEN2(ch_samplesheet)
    ch_versions = ch_versions.mix(KRAKEN2.out.versions)

    emit:
    asv_table = DADA2_DENOISING.out.asv_table // CRITICAL: This goes to PICRUSt2 
    versions  = ch_versions                   // This goes to the Master versions collector [cite: 8]
}