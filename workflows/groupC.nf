// groupC.nf: phase2_diversity_analysis

include { CONTAMINATION     } from '../modules/groupC/phase1_system_validation/contamination.nf'
include { DIVERSITY_METRICS } from '../modules/groupC/phase2_diversity_analysis/modules/module_diversity_metrics/diversity_analysis.nf'
include { RAREFACTION_THRESHOLD;
          SKIP_RAREFACTION  } from '../modules/groupC/phase2_diversity_analysis/modules/module_rarefaction/select_rarefaction_batches.nf'

// Subworkflow: PHASE2_DIVERSITY
workflow PHASE2_DIVERSITY {

    take:
    asv_table       // path — asv_table.tsv        (from GROUPB.out.table_counts)
    rep_seqs_fasta  // path — rep-seqs.fasta        (from GROUPB.out.rep_seqs)
    metadata        // path — sample-metadata.tsv   (from GROUPA.out.metadata)

    main:

    // Module 1: diversity metrics
    DIVERSITY_METRICS(asv_table, rep_seqs_fasta, metadata)

    // The rarefied_table.qza produced by core-metrics is what module 2 operates on
    rarefied_table = DIVERSITY_METRICS.out.qza_files
        .flatten()
        .filter { it.name == 'rarefied_table.qza' }

    // Pass the parent directory of the shannon TSV as the diversity_table dir
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

    // ── Group D output channels ───────────────────────────────────
    ch_alpha = RAREFACTION_THRESHOLD.out.shannon
        .join(RAREFACTION_THRESHOLD.out.observed_features)
        .join(RAREFACTION_THRESHOLD.out.faith_pd)
        .join(RAREFACTION_THRESHOLD.out.simpson)
        .map { meta, shannon, observed, faith, simpson ->
            tuple(shannon, observed, faith, simpson)
        }

    ch_beta = RAREFACTION_THRESHOLD.out.weighted_unifrac
        .join(RAREFACTION_THRESHOLD.out.bray_curtis)
        .map { meta, unifrac, bray ->
            tuple(unifrac, bray)
        }

    ch_rarefaction = RAREFACTION_THRESHOLD.out.plots
        .map { meta, plots -> plots }

    emit:
    // ── CHANNEL OUTPUTS FOR GROUP D ──────────────────────────────
    OutputMetricResultsAlpha       = ch_alpha
    OutputMetricResultsBeta        = ch_beta
    OutputMetricResultsRarefaction = ch_rarefaction

    // ── MODULE 1 OUTPUTS (pre-rarefaction) ───────────────────────
    rooted_tree          = DIVERSITY_METRICS.out.rooted_tree
    sampling_depth       = DIVERSITY_METRICS.out.sampling_depth
    m1_shannon           = DIVERSITY_METRICS.out.shannon
    m1_observed          = DIVERSITY_METRICS.out.observed
    m1_faith             = DIVERSITY_METRICS.out.faith
    m1_simpson           = DIVERSITY_METRICS.out.simpson
    m1_weighted_unifrac  = DIVERSITY_METRICS.out.weighted_unifrac
    m1_bray_curtis       = DIVERSITY_METRICS.out.bray_curtis
    m1_qza_files         = DIVERSITY_METRICS.out.qza_files
    m1_qzv_files         = DIVERSITY_METRICS.out.qzv_files

    // ── MODULE 2 OUTPUTS (post-rarefaction) ──────────────────────
    rarefaction_threshold = RAREFACTION_THRESHOLD.out.threshold
    rarefaction_plots     = RAREFACTION_THRESHOLD.out.plots
    rarefaction_curves    = RAREFACTION_THRESHOLD.out.curves
    rarefaction_qc        = RAREFACTION_THRESHOLD.out.qc
    rarefaction_report    = RAREFACTION_THRESHOLD.out.report
    depth_used            = RAREFACTION_THRESHOLD.out.depth_used
    shannon               = RAREFACTION_THRESHOLD.out.shannon
    observed              = RAREFACTION_THRESHOLD.out.observed_features
    faith                 = RAREFACTION_THRESHOLD.out.faith_pd
    simpson               = RAREFACTION_THRESHOLD.out.simpson
    weighted_unifrac      = RAREFACTION_THRESHOLD.out.weighted_unifrac
    bray_curtis           = RAREFACTION_THRESHOLD.out.bray_curtis
    qza_files             = RAREFACTION_THRESHOLD.out.qza_files
    qzv_files             = RAREFACTION_THRESHOLD.out.qzv_files
    skip_report           = SKIP_RAREFACTION.out.report
}


workflow GROUPC {

    take:
    ch_table_counts  // path: asv_table.tsv       from GROUPB.out.table_counts
    ch_rep_seqs      // path: rep_seqs.fasta       from GROUPB.out.rep_seqs
    ch_taxonomy      // path: ASV_taxonomy.tsv     from GROUPB.out.taxonomy
    ch_metadata      // path: sample-metadata.tsv  from GROUPA.out.metadata

    main:

    // Phase 1: contamination filtering
    CONTAMINATION(
        ch_table_counts,
        ch_taxonomy,
        ch_rep_seqs,
        ch_metadata
    )

    // Phase 2: diversity analysis
    PHASE2_DIVERSITY(
        ch_table_counts,  
        ch_rep_seqs,
        ch_metadata
    )

    emit:
    // ── OUTPUTS FOR GROUP D ───────────────────────────────────────
    contamination_summary          = CONTAMINATION.out.summary_png
    OutputMetricResultsAlpha       = PHASE2_DIVERSITY.out.OutputMetricResultsAlpha
    OutputMetricResultsBeta        = PHASE2_DIVERSITY.out.OutputMetricResultsBeta
    OutputMetricResultsRarefaction = PHASE2_DIVERSITY.out.OutputMetricResultsRarefaction

    // ── ALL PHASE 2 OUTPUTS ───────────────────────────────────────
    rooted_tree           = PHASE2_DIVERSITY.out.rooted_tree
    sampling_depth        = PHASE2_DIVERSITY.out.sampling_depth
    m1_shannon            = PHASE2_DIVERSITY.out.m1_shannon
    m1_observed           = PHASE2_DIVERSITY.out.m1_observed
    m1_faith              = PHASE2_DIVERSITY.out.m1_faith
    m1_simpson            = PHASE2_DIVERSITY.out.m1_simpson
    m1_weighted_unifrac   = PHASE2_DIVERSITY.out.m1_weighted_unifrac
    m1_bray_curtis        = PHASE2_DIVERSITY.out.m1_bray_curtis
    m1_qza_files          = PHASE2_DIVERSITY.out.m1_qza_files
    m1_qzv_files          = PHASE2_DIVERSITY.out.m1_qzv_files
    rarefaction_threshold = PHASE2_DIVERSITY.out.rarefaction_threshold
    rarefaction_plots     = PHASE2_DIVERSITY.out.rarefaction_plots
    rarefaction_curves    = PHASE2_DIVERSITY.out.rarefaction_curves
    rarefaction_qc        = PHASE2_DIVERSITY.out.rarefaction_qc
    rarefaction_report    = PHASE2_DIVERSITY.out.rarefaction_report
    depth_used            = PHASE2_DIVERSITY.out.depth_used
    shannon               = PHASE2_DIVERSITY.out.shannon
    observed              = PHASE2_DIVERSITY.out.observed
    faith                 = PHASE2_DIVERSITY.out.faith
    simpson               = PHASE2_DIVERSITY.out.simpson
    weighted_unifrac      = PHASE2_DIVERSITY.out.weighted_unifrac
    bray_curtis           = PHASE2_DIVERSITY.out.bray_curtis
    qza_files             = PHASE2_DIVERSITY.out.qza_files
    qzv_files             = PHASE2_DIVERSITY.out.qzv_files
    skip_report           = PHASE2_DIVERSITY.out.skip_report
}


workflow {

    asv_table      = Channel.value(file("${params.data_dir}/asv_table.tsv"))
    rep_seqs_fasta = Channel.value(file("${params.data_dir}/rep-seqs.fasta"))
    metadata       = Channel.value(file("${params.data_dir}/sample-metadata.tsv"))

    // Testing PHASE2_DIVERSITY standalone (no contamination step)
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
