nextflow.enable.dsl=2

// Module imports
include { FASTQC as FASTQC_RAW }  from './modules/a_fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/a_fastqc.nf'
include { CUTADAPT }              from './modules/a_cutadapt.nf'
include { MULTIQC }               from './modules/a_multiqc.nf'
include { CLEAN_MULTIQC }         from './modules/a_clean_multiqc.nf'
include { NEW_SAMPLE_SHEET }      from './modules/a_new_sample_sheet.nf'

// Initialize error log file
file("${params.out_dir}/errors").mkdirs()
def error_log = file("${params.out_dir}/errors/skipped_samples_log.tsv")
error_log.text = "Sample_ID\tError_Type\tExpected_Path\n"

// Main workflow engine
workflow {

    // Read CSV and check the files
    raw_reads_ch = Channel.fromPath(params.sample_sheet)
        // Splits the file by colons
        .splitCsv(sep: ',', header: true) 
        // Map the file to extract the sample_id and fastq file location
        .map { row ->
            def sra_id = row.sample_id
            def fastq_file = file("${row.directory}/${sra_id}.fastq")
            return tuple(sra_id, fastq_file)
        }
        // Apply error detection
            // First, does the file exist?
        .filter { sra_id, fastq_file ->
            if (!sra_id) return false
            if (!fastq_file.exists()) {
                log.warn "Skipping sample [${sra_id}]: File does not exist"
                error_log.append("${sra_id}\tMissing_File\t${fastq_file}\n")
                return false
            }
            // Second, is the file empty?
            if (fastq_file.size() == 0) {
                log.warn "Skipping sample [${sra_id}]: File is empty (0 bytes)"
                error_log.append("${sra_id}\tEmpty_File\t${fastq_file}\n")
                return false
            }
            return true
        }

    // Apply the FASTQC on the raw reads
    fastqc_raw_in = raw_reads_ch.map { id, read -> tuple(id, 'raw', read) }
    FASTQC_RAW(fastqc_raw_in)

    // Apply cutadapt to trim the reads
    CUTADAPT(raw_reads_ch)

    // Apply the FASTQC on the trimmed reads
    fastqc_trim_in = CUTADAPT.out.trimmed_reads.map { id, read -> tuple(id, 'trimmed', read) }
    FASTQC_TRIM(fastqc_trim_in)

    // Aggregate all the FASTQC results to pipe them into the MULTIQC
    all_qc_files_ch = FASTQC_RAW.out.qc_files
        .mix(FASTQC_TRIM.out.qc_files)
        .map { id, stage, files -> files }
        .collect()

    // Apply the MULTIQC to all the FASTQC
    MULTIQC(all_qc_files_ch)

    // Apply the cleaning of the MULTIQC file, extracting specific data and applying some basic interpretation
    CLEAN_MULTIQC(MULTIQC.out.fastqc_txt)

    // Gather all the data needed to create the new sample sheet
    successful_ids_file = CUTADAPT.out.trimmed_reads
        .map { sample_id, fastq -> sample_id }
        .collectFile(name: 'successful_samples.txt', newLine: true)

    // Generatee a new sample sheet for the group B
    NEW_SAMPLE_SHEET(successful_ids_file)

    // Generate the channels for the following group
    
    // Group all the trimmed files into a single folder channel
    pre_ch_trimmed_seqs = CUTADAPT.out.trimmed_reads.map { id, fastq -> fastq }.collect()

    emit:

    // ch_trimmed: folder of all the trimmed sequences (not published)
    ch_trimmed = pre_ch_trimmed_seqs
    
    // ch_sample_sheet_B: the generated sample sheet (not published)
    ch_sample_sheet_B = NEW_SAMPLE_SHEET.out.final_sheet
    
    // ch_metadata: reads the file from params without publishing it
    ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)
    
    // ch_quality_report: passes the multiqc report and your custom report
    ch_quality_report = MULTIQC.out.report.mix(CLEAN_MULTIQC.out.final_report).collect()
}