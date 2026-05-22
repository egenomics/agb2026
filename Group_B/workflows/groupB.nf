// workflows/groupB.nf

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT ALL GROUP B MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// QC Modules
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

// Taxonomic Profiling Modules
include { DADA2_FILTNTRIM   } from '../modules/nf-core/dada2/filtandtrim/main'
include { DADA2_ERR         } from '../modules/nf-core/dada2/learnerrors/main'
include { DADA2_DENOISING             } from '../modules/local/dada2/dada2_denoising'
include { DADA2_TAXONOMY              } from '../modules/local/dada2/dada2_taxonomy'
include { KRAKEN2_KRAKEN2 as KRAKEN2  } from '../modules/nf-core/kraken2/kraken2/main'

// Functional Annotation Modules
include { DADA2_EXPORT                } from '../modules/local/dada2/dada2_export'
include { PICRUST                     } from '../modules/local/picrust'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN GROUP B WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow GROUPB {

    take:
    ch_samplesheet             // channel: [ val(meta), [ reads ] ]
    multiqc_config             // path: multiqc_config
    multiqc_logo               // path: multiqc_logo
    ch_collated_versions       // channel: [ path(versions.yml) ]
    ch_methods_description     // channel: [ path(methods_description_mqc.yaml) ]
    ch_workflow_summary        // channel: [ path(workflow_summary_mqc.yaml) ]
    dada2_train_set            // path: SILVA train set
    dada2_species_set          // path: SILVA species set

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // ==========================================
    // STEP 1: QUALITY CONTROL (QC_CHECKS)
    // ==========================================
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })
    // ch_versions = ch_versions.mix(FASTQC.out.versions)

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions, ch_workflow_summary, ch_methods_description)
    
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [[id: 'alphaflow_qc'], files, multiqc_config, multiqc_logo, [], []]
        }
    )

    // ==========================================
    // STEP 2: TAXONOMIC PROFILING
    // ==========================================
    // 1. Package the truncation thresholds as lists to satisfy the module's internal [1] indexing
    ch_filt_input = ch_samplesheet.map { meta, reads ->
        tuple(meta, reads, [0, 140], [0, 0])
    }

    // 2. Filter and Trim using the updated, structurally perfect data stream
    DADA2_FILTNTRIM(ch_filt_input)
    
    ch_filtered_reads = DADA2_FILTNTRIM.out.reads_logs_args
        .map { meta, reads, stats, args -> tuple(meta, reads) }

    // 3. Learn the error rates from the filtered reads
    DADA2_ERR(ch_filtered_reads)

    // Combine the reads and the error models by their 'meta' key into ONE channel
    ch_denoising_input = ch_filtered_reads.join(DADA2_ERR.out.errormodel)    

    // 3. Denoise with BOTH the filtered reads and the learned error model!
    DADA2_DENOISING(ch_denoising_input)

    // Pass the denoised ASV table and the cluster database paths to taxonomy
    DADA2_TAXONOMY(
        DADA2_DENOISING.out.asv_table.map { meta, rds -> rds },
        file(dada2_train_set),
        "ASV_taxonomy.tsv",
        "Kingdom,Phylum,Class,Order,Family,Genus,Species"
    )

    // Example of standard 4-argument alignment for nf-core Kraken2
    //KRAKEN2(
	//ch_samplesheet,
	//file(params.kraken2_db), // Or path to your Kraken2 database directory
	//false,                   // save_output_fastqs
	//false                    // open_reading_frames
	//)
    // ch_versions = ch_versions.mix(KRAKEN2.out.versions)

    // ==========================================
    // STEP 3: FUNCTIONAL ANNOTATION
    // ==========================================
    // Convert DADA2 seqtab RDS to FASTA + TSV for PICRUSt2
    DADA2_EXPORT(DADA2_DENOISING.out.asv_table)

    PICRUST(
        DADA2_EXPORT.out.fasta,
        DADA2_EXPORT.out.table,
        "metagenome",
        true
    )
    // ch_versions = ch_versions.mix(PICRUST.out.versions)

    emit:
    // DELIVERABLES FOR GROUP C
    table_counts = DADA2_DENOISING.out.asv_table     // Becomes table_counts.tsv
    rep_seqs     = DADA2_DENOISING.out.denoised      // Becomes rep_seqs.fasta
    
    // GENERAL OUTPUTS
    multiqc_report = MULTIQC.out.report
    pathways       = PICRUST.out.pathways
    versions       = ch_versions
}
