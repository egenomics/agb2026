nextflow.enable.dsl=2

// --- MODULE IMPORTS ---
include { FASTQC as FASTQC_RAW }  from './modules/a_fastqc.nf'
include { FASTQC as FASTQC_TRIM } from './modules/a_fastqc.nf'
include { CUTADAPT }              from './modules/a_cutadapt.nf'
include { MULTIQC }               from './modules/a_multiqc.nf'
include { CLEAN_MULTIQC }         from './modules/a_clean_multiqc.nf'

// --- INITIALIZE ERROR LOG ---
// Create a log file in the output directory to track skipped samples
def error_log = file("${params.out_dir}/skipped_samples_log.tsv")
// Initialize with headers (overwrites any previous log to keep it fresh)
error_log.text = "Sample_ID\tError_Type\tExpected_Path\n"

// --- MAIN WORKFLOW ENGINE ---
workflow {

    // 1. Read the CSV file, locate fastq, filter, and log errors
    raw_reads_ch = Channel.fromPath(params.sample_sheet)
        .splitCsv(sep: ',', header: true) // Using comma separator for your new format
        .map { row ->
            def sra_id = row.sample_id
            def fastq_file = file("${row.directory}/${sra_id}.fastq")
            return tuple(sra_id, fastq_file)
        }
        .filter { sra_id, fastq_file ->
            if (!sra_id) return false // Ignore completely blank rows
            
            if (!fastq_file.exists()) {
                log.warn "Skipping sample [${sra_id}]: File does not exist"
                // Append the error to our log file
                error_log.append("${sra_id}\tMissing_File\t${fastq_file}\n")
                return false
            }
            if (fastq_file.size() == 0) {
                log.warn "Skipping sample [${sra_id}]: File is empty (0 bytes)"
                // Append the error to our log file
                error_log.append("${sra_id}\tEmpty_File\t${fastq_file}\n")
                return false
            }
            return true
        }

    // 2. FASTQC on Raw Reads
    fastqc_raw_in = raw_reads_ch.map { id, read -> tuple(id, 'raw', read) }
    FASTQC_RAW(fastqc_raw_in)

    // 3. CUTADAPT to trim primers (Single-end version)
    CUTADAPT(raw_reads_ch)

    // 4. FASTQC on Trimmed Reads
    fastqc_trim_in = CUTADAPT.out.trimmed_reads.map { id, read -> tuple(id, 'trimmed', read) }
    FASTQC_TRIM(fastqc_trim_in)

    // 5. Aggregate ALL FastQC reports for MultiQC
    all_qc_files_ch = FASTQC_RAW.out.qc_files
        .mix(FASTQC_TRIM.out.qc_files)
        .map { id, stage, files -> files }
        .collect()

    // 6. MULTIQC generates the merged report
    MULTIQC(all_qc_files_ch)

    // 7. CLEAN_MULTIQC runs your Python script on the summary report
    CLEAN_MULTIQC(MULTIQC.out.fastqc_txt)
}