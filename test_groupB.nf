// test_groupB.nf
nextflow.enable.dsl=2

// --- MANUALLY ADD DATABASE PATHS HERE ---
params.dada2_train_set   = "${projectDir}/biodb/dada2/silva_nr99_v138.2_toGenus_trainset.fa"
params.dada2_species_set = "${projectDir}/biodb/dada2/silva_v138.2_assignSpecies.fa"
params.kraken2_db        = "${projectDir}/biodb" // Points to your biodb directory as a placeholder

// Import unified workflow
include { GROUPB } from './workflows/groupB'

workflow {
    // 1. Define standard channels from params (matching nextflow.config)
    //
    // STANDALONE MODE: build the channel from a CSV of pre-trimmed
    // samplesheet.csv generated individually using generate_samplesheet.sh.
    // Used for testing Group B in isolation
    // data (e.g. paquete_cluster/seqs/trimmed/*.fastq.gz)
    ch_samplesheet = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [ id: row.sample, single_end: row.fastq_2 ? false : true, strandedness: row.strandedness ]
            def reads = row.fastq_2 ? [ file(row.fastq_1), file(row.fastq_2) ] : [ file(row.fastq_1) ]
            return [ meta, reads ]
        }

    // ────────────────────────────────────────────────────────────────────────
    // INTEGRATED MODE (future, in main.nf): when Group A's workflow exists,
    // replace the CSV-parsing block above with a direct channel handoff.
    // Group A's `trimmed_reads` output has the SAME shape — [meta, reads] —
    // so it drops straight in as Group B's input.
    //
    //
    //
    //     include { GROUPA } from './workflows/groupA'
    //
    //     GROUPA(ch_raw_samplesheet, ch_metadata)
    //     ch_samplesheet = GROUPA.out.trimmed_reads
    //         // emits: [ [id:'ERR1074192', single_end:true, strandedness:'forward'],
    //         //          [file('.../ERR1074192_trim.fastq.gz')] ]
    //         // ...one tuple per sample, identical topology to the CSV path above.
    // ────────────────────────────────────────────────────────────────────────

    // 2. Mock the standard empty or placeholder assets that MultiQC expects
    ch_multiqc_config      = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: false)
    ch_multiqc_logo        = file("${projectDir}/assets/multiqc_logo.png", checkIfExists: false)
    ch_versions            = Channel.empty()
    ch_methods_description = Channel.empty()
    ch_workflow_summary    = Channel.empty()

    // 3. Fire up Group B's workflow completely standalone!
    GROUPB (
        ch_samplesheet,            // ← in integrated mode this is GROUPA.out.trimmed_reads
        ch_multiqc_config,
        ch_multiqc_logo,
        ch_versions,
        ch_methods_description,
        ch_workflow_summary,
        file(params.dada2_train_set),
        file(params.dada2_species_set)
    )
}
