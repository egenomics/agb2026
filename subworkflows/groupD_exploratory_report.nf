//Process Import
include { DEMOGRAPHIC_SUMMARY }         from '../modules/groupD/exploratory_report/demo_table' 
include { METADATA_REPORT }             from '../modules/groupD/exploratory_report/co_heat'
include { PCA_MICROBIOME }              from '../modules/groupD/exploratory_report/pca_samples_bacteria'
include { PARALLEL_PLOT_MICROBIOME }    from '../modules/groupD/exploratory_report/parallel_plot_relabu_healthygrouped'
include { HEATMAP_MICROBIOME }          from '../modules/groupD/exploratory_report/heatmap_samples_bacteria'
include { VOLCANO_PLOT }                from '../modules/groupD/exploratory_report/volcano'
include { BETA_DIVERSITY_TREE }         from '../modules/groupD/exploratory_report/tree_beta'

workflow EXPLORATORY_REPORT{
    take:
    ch_metadata     // channel: Metadata from Group A
    ch_counts       // channel: Path to "annotated_table_counts.tsv"
    ch_taxonomy     // channel: Path to "annotated_taxonomy.tsv"
    ch_bray         // channel: Path to beta diversity"bray.tsv"

    main:
    // ---------------------------------------------------------
    // 1. Prepare and Demographics summary
    // ---------------------------------------------------------
    // Only metadata required
    DEMOGRAPHIC_SUMMARY(ch_metadata)

    // ---------------------------------------------------------
    // 2. Prepare and Run Cohort Report
    // ---------------------------------------------------------
    // Only metadata required
    METADATA_REPORT(ch_metadata)

    // ---------------------------------------------------------
    // 3. Prepare and Run PCA of the Microbiome
    // ---------------------------------------------------------
    // It takes 3 arguments (counts, taxonomy, metadata) in that order
    PCA_MICROBIOME(ch_counts, ch_taxonomy, ch_metadata)

    // ---------------------------------------------------------
    // 4. Prepare and Run Parallel Microbiome Plot
    // ---------------------------------------------------------
    // It takes 3 arguments (counts, taxonomy, metadata) in that order
    PARALLEL_PLOT_MICROBIOME(ch_counts, ch_taxonomy, ch_metadata)

    // ---------------------------------------------------------
    // 5. Prepare and Run Microbiome Heatmap Plot
    // ---------------------------------------------------------
    // It takes 3 arguments (counts, taxonomy, metadata) in that order
    HEATMAP_MICROBIOME(ch_counts, ch_taxonomy, ch_metadata)

    // ---------------------------------------------------------
    // 6. Prepare and Run Volcano sample vs bacteria
    // ---------------------------------------------------------
    // It takes 3 arguments (metadata, counts, taxonomy) in that order
    VOLCANO_PLOT(ch_metadata, ch_counts, ch_taxonomy)

    // ---------------------------------------------------------
    // 7. Prepare and Run Microbiome Heatmap Plot
    // ---------------------------------------------------------
    // It takes a tuple with metadata and Bray
    ch_tree = ch_metadata.combine(ch_bray)
    BETA_DIVERSITY_TREE(ch_tree)

    emit:
    // Demographic Summary
    demo_table_png   = DEMOGRAPHIC_SUMMARY.out.demo_table_png

    // Cohort Report
    metadata_report  = METADATA_REPORT.out.co_heat 

    // PCA Outputs
    pca_results      = PCA_MICROBIOME.out.pca_results
    pca_object       = PCA_MICROBIOME.out.pca_object
    pca_genus_counts = PCA_MICROBIOME.out.pca_counts

    // Parallel Output
    parallel_plot    = PARALLEL_PLOT_MICROBIOME.out.parallel_results

    // Heatmap Output
    heatmap_plot     = HEATMAP_MICROBIOME.out.heatmap_results

    // Volcano Plot
    volcano_plot     = VOLCANO_PLOT.out.volcano_png

    // Beta Tree
    beta_tree        = BETA_DIVERSITY_TREE.out.tree_png

}