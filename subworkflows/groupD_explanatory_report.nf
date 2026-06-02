include { alpha_diversity_status } from '../modules/groupD/explanatory_report/alpha_diversity_status'
include { PCoA_plots }             from '../modules/groupD/explanatory_report/PCoA' 
include { overview_table }         from '../modules/groupD/explanatory_report/overview_table' 
include { VIRULENCE_PLOT }         from '../modules/groupD/explanatory_report/virulence_report'
include { VIOLIN_PLOTS }           from '../modules/groupD/explanatory_report/violin_plot'

workflow EXPLANATORY_REPORT{
    take:
    ch_metadata
    ch_counts
    ch_taxonomy
    ch_alphas
    ch_pca
    ch_genus_counts

    main:
    // ---------------------------------------------------------
    // 1. Prepare and Run Alpha Diversity
    // ---------------------------------------------------------
    // Combined with metadata
    ch_input_for_alpha = ch_metadata.combine(ch_alphas)
    
    alpha_diversity_status(ch_input_for_alpha)

    // ---------------------------------------------------------
    // 2. Prepare and Run PCoA
    // ---------------------------------------------------------
    // PCoA process expects: tuple path(pca), path(metadata)
    //Change for the PCA
    ch_input_for_pcoa = ch_pca.combine(ch_metadata)
    
    PCoA_plots(ch_input_for_pcoa)

    // ---------------------------------------------------------
    // 3. Prepare and run the overview table
    // ---------------------------------------------------------
    // overview process expects: tuple path(metadata), path(pca), path(z_scores), path(genus_counts)
    //Change for the PCA and annotated genus
    ch_overview = ch_metadata
        .combine(ch_pca)
        .combine(alpha_diversity_status.out.zscore_status)
        .combine(ch_genus_counts)

    overview_table(ch_overview)

    // ---------------------------------------------------------
    // 4. Prepare and run the Virulence plot
    // ---------------------------------------------------------
    // 3 arguments (meta, counts, taxonomy)
    VIRULENCE_PLOT(ch_metadata, ch_counts, ch_taxonomy)

    // ---------------------------------------------------------
    // 5. Prepare and run the Violin plot
    // ---------------------------------------------------------
    // 3 arguments (meta, counts, taxonomy)
    VIOLIN_PLOTS(ch_metadata, ch_counts, ch_taxonomy)

    emit:
    overview_table        = overview_table.out.overview_sample_table_dir
    alpha_plots           = alpha_diversity_status.out.alpha_div_dist_dir
    pcoa_plots            = PCoA_plots.out.PCoA_patient_plots_dir
    patogeny_plots        = VIRULENCE_PLOT.out.virulence_report
    violin_plots          = VIOLIN_PLOTS.out.violin_report
}