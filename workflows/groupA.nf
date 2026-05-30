nextflow.enable.dsl=2

// Module imports
include { FASTQC as FASTQC_RAW }  from '../modules/groupA/a_fastqc.nf'
include { FASTQC as FASTQC_TRIM } from '../modules/groupA/a_fastqc.nf'
include { CUTADAPT }              from '../modules/groupA/a_cutadapt.nf'
include { MULTIQC }               from '../modules/groupA/a_multiqc.nf'
include { CLEAN_MULTIQC }         from '../modules/groupA/a_clean_multiqc.nf'
// include { NEW_SAMPLE_SHEET }      from '../modules/groupA/a_new_sample_sheet.nf'   // disabled: real handoff is the trimmed_reads channel

// Group A: raw-read QC + primer trimming. Single-end only (CUTADAPT trims
// reads[0] with the forward primer)
// meta.single_end is carried through so the rest of the pipeline stays paired-end capable.
workflow GROUPA {

    take:
    ch_samplesheet // channel: [ val(meta), [ reads ] ] --> raw reads
    ch_metadata // file: patient metadata .tsv (passthrough to Group C)

    main:
    ch_versions = channel.empty()

    // Read CSV and check the files
    // The CSV is parsed in PIPELINE_INITIALISATION and arrives as ch_samplesheet [meta, [reads]].
    ch_valid = ch_samplesheet
        .filter { meta, reads ->
            def fastq_file = reads[0]
            // First, does the file exist?
            if (fastq_file == null || !fastq_file.exists()) {
                log.warn "Skipping sample [${meta.id}]: File does not exist"
                return false
            }
            // Second, is the file empty?
            if (fastq_file.size() == 0) {
                log.warn "Skipping sample [${meta.id}]: File is empty (0 bytes)"
                return false
            }
            return true
        }

    // Audit log of skipped samples
    ch_samplesheet
        .map { meta, reads ->
            def fastq_file = reads[0]
            if (fastq_file == null || !fastq_file.exists()) return "${meta.id}\tMissing_File\t${fastq_file}"
            if (fastq_file.size() == 0)                     return "${meta.id}\tEmpty_File\t${fastq_file}"
            return null
        }
        .filter { line -> line != null }
        .collectFile(
            name:     'skipped_samples_log.tsv',
            storeDir: "${params.outdir}/groupA/errors",
            seed:     "Sample_ID\tError_Type\tExpected_Path",
            newLine:  true,
            sort:     true
        )

    // Apply the FASTQC on the raw reads
    // CHANGED: input is now [meta, [reads]]; mapped to the module's (id, stage, read) shape.
    FASTQC_RAW( ch_valid.map { meta, reads -> tuple(meta.id, 'raw', reads[0]) } )

    // Apply cutadapt to trim the reads
    CUTADAPT( ch_valid.map { meta, reads -> tuple(meta.id, reads[0]) } )

    // Apply the FASTQC on the trimmed reads
    FASTQC_TRIM( CUTADAPT.out.trimmed_reads.map { id, read -> tuple(id, 'trimmed', read) } )

    // Aggregate all the FASTQC results to pipe them into the MULTIQC
    all_qc_files_ch = FASTQC_RAW.out.qc_files
        .mix(FASTQC_TRIM.out.qc_files)
        .map { _id, _stage, files -> files }
        .collect()

    // Apply the MULTIQC to all the FASTQC
    MULTIQC( all_qc_files_ch )

    // Clean the MULTIQC file: extract specific data and apply some basic interpretation
    CLEAN_MULTIQC( MULTIQC.out.fastqc_txt )

    // NEW_SAMPLE_SHEET disabled: the real A→B handoff
    // is the trimmed_reads channel below.
    //
    // successful_ids_file = CUTADAPT.out.trimmed_reads
    //     .map { sample_id, _fastq -> sample_id }
    //     .collectFile(name: 'successful_samples.txt', newLine: true)
    // NEW_SAMPLE_SHEET( successful_ids_file )

    // Generate the channel for the following group (Group B).
    // Mofication: the old code was `CUTADAPT.out.trimmed_reads.map{ id, fastq -> fastq }.collect()`,
    // which DROPPED per-sample meta and collapsed everything into ONE emission. Group B fans
    // out per sample and needs meta (especially single_end), so we rejoin meta on id and keep
    // one [meta, [reads]] tuple per sample.
    ch_trimmed_reads = ch_valid
        .map { meta, _reads -> tuple(meta.id, meta) }
        .join( CUTADAPT.out.trimmed_reads )
        .map { _id, meta, fastq -> tuple(meta, [ fastq ]) }

    emit:
    // CHANGED emit names to match main.nf. Old -> new:
    //   ch_trimmed        -> trimmed_reads   (now [meta,[reads]] per sample, not a single collected folder)
    //   ch_metadata       -> metadata        (now received via `take` from main.nf, not re-read from params)
    //   ch_quality_report -> quality_report
    //   ch_sample_sheet_B -> (dropped: NEW_SAMPLE_SHEET still runs, but the CSV is no longer emitted)
    trimmed_reads  = ch_trimmed_reads                                                  // [meta, [reads]] — for Group B
    metadata       = ch_metadata                                                       // file — for Group C
    quality_report = MULTIQC.out.report.mix(CLEAN_MULTIQC.out.final_report).collect()  // MultiQC + custom report
    versions       = ch_versions                                                       // software versions
}
