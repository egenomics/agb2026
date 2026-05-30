// main.nf: phase2_diversity_analysis

include {DIVERSITY_METRICS} from '../modules/groupC/phase2_diversity_analysis/modules/module_diversity_metrics/diversity_analysis.nf'
include {RAREFACTION_THRESHOLD;
        SKIP_RAREFACTION} from '../modules/groupC/phase2_diversity_analysis/modules/module_rarefaction/select_rarefaction_batches.nf'

// Subworkflow: PHASE2_DIVERSITY
// Runs module 1 (diversity metrics) then decides whether to run module 2 (rarefaction) based on the minimum sample count across alpha diversity TSVs.

workflow PHASE2_DIVERSITY {

    take:
    asv_table       // path — asv_table.tsv
    rep_seqs_fasta  // path — rep-seqs.fasta
    metadata        // path — sample-metadata.tsv

    main:

    // Module 1: diversity metrics 
    DIVERSITY_METRICS(asv_table, rep_seqs_fasta, metadata)

    // Sample count check (done in the driver workflow, not here)

    // The rarefied_table.qza produced by core-metrics is what module 2 operates on
    rarefied_table = DIVERSITY_METRICS.out.qza_files
        .flatten()
        .filter { it.name == 'rarefied_table.qza' }

    // Pass the parent directory of the shannon TSV as the diversity_table dir.
    diversity_table_dir = DIVERSITY_METRICS.out.shannon
        .flatten()
        .first()
        .map { it.parent.parent }

    input_ch = rarefied_table
        .combine(DIVERSITY_METRICS.out.rooted_tree)
        .combine(metadata)
        .combine(diversity_table_dir)
        .combine(DIVERSITY_METRICS.out.sampling_depth)
        .map { table, phylo, meta_file, div_dir, depth_file ->
            tuple([id: "rarefaction"], table, phylo, meta_file, div_dir, depth_file)
        }

    // Module 2: rarefaction (conditional on sample count)
    // Count the minimum samples across all alpha diversity TSVs.
    // Uses the on-disk TSVs emitted by module 1 so the count reflects exactly what made it through the rarefaction + filtering steps.
    tsv_subpaths = ['shannon/alpha-diversity.tsv',
                    'observed/alpha-diversity.tsv',
                    'faith/alpha-diversity.tsv',
                    'simpson/alpha-diversity.tsv']

    n_samples = diversity_table_dir.map { dir ->
        def counts = tsv_subpaths.collect { sub ->
            def f = file("${dir}/${sub}")
            f.exists()
                ? f.readLines().findAll { !it.startsWith('#') && !it.trim().isEmpty() }.size() - 1
                : null
        }.findAll { it != null }
        counts.min()
    }

    min_size = params.min_batch_size

    // Branch channel: enough samples: RAREFACTION_THRESHOLD, too few: SKIP
    input_ch
        .combine(n_samples)
        .branch {
            run:  it[-1] >= min_size
            skip: it[-1] <  min_size
        }
        .set { branched }

    RAREFACTION_THRESHOLD(
        branched.run.map { it[0..-2] }
    )

    SKIP_RAREFACTION(
        branched.skip.map { row ->
            def meta = row[0]
            def n    = row[-1]
            tuple(meta, n, min_size)
        }
    )

    emit:
    // Module 1 outputs
    rooted_tree      = DIVERSITY_METRICS.out.rooted_tree
    sampling_depth   = DIVERSITY_METRICS.out.sampling_depth
    // module 1 alpha diversity (pre-rarefaction threshold selection)
    m1_shannon       = DIVERSITY_METRICS.out.shannon
    m1_observed      = DIVERSITY_METRICS.out.observed
    m1_faith         = DIVERSITY_METRICS.out.faith
    m1_simpson       = DIVERSITY_METRICS.out.simpson
    m1_weighted_unifrac = DIVERSITY_METRICS.out.weighted_unifrac
    m1_bray_curtis   = DIVERSITY_METRICS.out.bray_curtis
    m1_qza_files     = DIVERSITY_METRICS.out.qza_files
    m1_qzv_files     = DIVERSITY_METRICS.out.qzv_files

    // Module 2 outputs (only populated when enough samples)
    rarefaction_threshold = RAREFACTION_THRESHOLD.out.threshold
    rarefaction_plots     = RAREFACTION_THRESHOLD.out.plots
    rarefaction_curves    = RAREFACTION_THRESHOLD.out.curves
    rarefaction_qc        = RAREFACTION_THRESHOLD.out.qc
    rarefaction_report    = RAREFACTION_THRESHOLD.out.report
    depth_used            = RAREFACTION_THRESHOLD.out.depth_used
    // final alpha diversity (post-rarefaction threshold selection)
    shannon               = RAREFACTION_THRESHOLD.out.shannon
    observed              = RAREFACTION_THRESHOLD.out.observed_features
    faith                 = RAREFACTION_THRESHOLD.out.faith_pd
    simpson               = RAREFACTION_THRESHOLD.out.simpson
    weighted_unifrac      = RAREFACTION_THRESHOLD.out.weighted_unifrac
    bray_curtis           = RAREFACTION_THRESHOLD.out.bray_curtis
    qza_files             = RAREFACTION_THRESHOLD.out.qza_files
    qzv_files             = RAREFACTION_THRESHOLD.out.qzv_files
    // skip report (only populated when not enough samples)
    skip_report           = SKIP_RAREFACTION.out.report
}


workflow {
    asv_table      = Channel.value(file("${params.data_dir}/asv_table.tsv"))
    rep_seqs_fasta = Channel.value(file("${params.data_dir}/rep-seqs.fasta"))
    metadata       = Channel.value(file("${params.data_dir}/sample-metadata.tsv"))

    PHASE2_DIVERSITY(asv_table, rep_seqs_fasta, metadata)

    PHASE2_DIVERSITY.out.sampling_depth
        .view { f -> "Sampling depth         : ${f.text.trim()}" }

    PHASE2_DIVERSITY.out.depth_used
        .mix( Channel.value( [null, null] ) )
        .view { meta, f -> f ? "Rarefaction depth used : ${f.text.trim()}" : "Rarefaction depth used : [skipped]" }

    PHASE2_DIVERSITY.out.rarefaction_plots
        .mix( Channel.value( [null, null] ) )
        .view { meta, f -> f ? "Rarefaction plots      : ${f}" : "Rarefaction plots      : [skipped]" }

    PHASE2_DIVERSITY.out.skip_report
        .mix( Channel.value( [null, null] ) )
        .view { meta, f -> f ? "Skip report            : ${f}" : "Skip report            : [not triggered]" }
}
