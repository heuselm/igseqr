#' @title differentialAbundanceTestMS
#'
#' @description test differential abundance of proteins, clones and clone groups in DIA-MS data on precursor level
#'
#' @param DIA_resultset DIA-MS resultset imported by importDIAresults function.
#' @param study_design_external Alternative study design overwriting the study design in DIA_resultset containing columns filename
#' @param comparison_matrix A matrix defining the comparisons to be run.
#' Column one contains the condition in the counter of each comparison, column two condition in the denominator of the comparison.
#' Positive log2FCs will mean higher signal in counter condition listed in column 1.
#' For possible values, check unique(DIA_resultset$study_design$condition).
#' The data will be subset to the conditions specified in the comparison matrix.
#' @param write_preprocessing_results Whether the intermediate results of each pairwise comparison shall be written. Separate subfolders will be generated. Default: FALSE
#' @param write_csv_protein Whether to write differential abundance testing result table summarized to protein level (csv). Default: TRUE
#' @param write_csv_precursor Whether to write differential abundance testing result table containing precursor level information (csv). Default: FALSE
#' @param imputation_percentile Percentile of non-0 quantitative values that the normal distribution for missing value imputation shall be centered on.
#'  Default: 0.001. Check imputation_histogram.pdf output for parameter selection.
#' @param imputation_rnorm_sd Standard deviation (~width) of the normal distribution for missing value imputation. Check imputation_histogram.pdf output for parameter selection.
#' @param overview_label_prefix Prefix label for the differential abundance overview .pdf and interactive .html plots
#' @param plot_pdf Whether to plot summary of the differential tests.
#' @param plot_html If plotting summary, whether to also generate interactive summary of the differential tests.
#' @param prot_highlight_tag tag that determines which protein.groups get highlighted
#' @return DIA_resultset with appended element $differentialAbundanceTestResults  with tables
#' $diffTestRes_prec
#' $diffTestRes_prot
#' $comparison_matrix
#' $data_source
#' $study_design
#' $data_long
#' $data_matrix_log2
#' $data_matrix_log2_imp
#'
#' @import data.table ggplot2 DiffTestR
#'
#' @export

differentialAbundanceTestMS = function(DIA_resultset,
                                       study_design_external = NULL,
                                       comparison_matrix = matrix(c("AP_m11_tF", "AP_m11_t0",
                                                                      "AP_m12_tF", "AP_m12_t0",
                                                                      "AP_m1_tF", "AP_m1_t0",
                                                                      "AP_m2_tF", "AP_m2_t0"),
                                                                    ncol = 2, byrow = TRUE),
                                      write_preprocessing_results = TRUE,
                                      write_diffTable_csv_protein = TRUE,
                                      write_diffTable_csv_precursor = TRUE,
                                      imputation_percentile = 0.001,
                                      imputation_rnorm_sd = 0.2,
                                      overview_label_prefix = "IGSeq_M1imm_APMS_batch1_t3t0",
                                      plot_pdf = TRUE,
                                      plot_html = TRUE,
                                      prot_highlight_tag = "IGSeq"){

  # Initialize table to collect Results across the differential tests
  cRes = data.table()

  # If needed, overwrite experimental design
  if(!is.null(study_design_external)){
    study_des = study_design_external
  } else {
    study_des = DIA_resultset$study_design
  }

  # Subset study design to those of interest in the current comparisons:
  study_des_oi = study_des[condition %in% comparison_matrix]
  pairs_oi = comparison_matrix

  for (i in 1:nrow(pairs_oi)){

    # say what's happening
    message("Running comparison ",i," of ", nrow(pairs_oi),":\n",
            pairs_oi[i,1], " vs ",  pairs_oi[i,2])

    diffTestRes = testDifferentialAbundance(input_dt = DIA_resultset$data_wide,
                                              study_design = study_des_oi,
                                              condition_1 = pairs_oi[i,1],
                                              condition_2 = pairs_oi[i,2],

                                              # Define normalization behavior
                                              normalize_data = TRUE,
                                              normalization_function = limma::normalizeQuantiles,

                                              # filtering (only global filtering implemented)
                                              min_n_obs = 4,
                                              # imputation of missing values
                                              imp_percentile = imputation_percentile,
                                              imp_sd = imputation_rnorm_sd,

                                              # plots?
                                              plot_pdf = if(i==1){TRUE}else{FALSE},
                                              # tsv result table?
                                              write_tsv_tables = FALSE,
                                              # highlight protein in volcano?
                                              target_protein = prot_highlight_tag)

    if(i==1 & write_preprocessing_results){
      saveRDS(diffTestRes, file = "differentialAbundanceTestMS_Example.rds")
      fwrite(diffTestRes$diffExpr_result_dt, file = "differentialAbundanceTestMS_ScaledIntensities.csv")
    }

    cRes = rbind(cRes, diffTestRes$diffExpr_result_dt)

  }
  if (write_diffTable_csv_precursor){
    fwrite(cRes, "differentialAbundanceTestMSResults_PrecursorLevel.csv")
  }

  cRes_prot = unique(cRes[, .(Protein.Group, Protein.Names, comparison, target_prot,
                              log2_fold_change_protein,n_precursors,p_value_protein,p_value_BHadj_protein)])
  if (write_diffTable_csv_protein){
    fwrite(cRes_prot, "differentialAbundanceTestMSResults_ProteinLevel.csv")
  }

  if (plot_pdf){
    diffOverview = DiffTestR::plotDifferentialAbundanceOverview(cRes,
                                                                label_prefix = overview_label_prefix,
                                                                browsable_html = plot_html,
                                                                protein_highlight_tag = prot_highlight_tag)

  }

  # Assemble results and append to input object
  res = list("diffTestRes_prec" = cRes,
             "diffTestRes_prot" = cRes_prot,
             "comparison_matrix" = comparison_matrix,
             "data_source" = diffTestRes$data_source,
             "study_design" = diffTestRes$study_design,
             "data_long" = diffTestRes$data_long,
             "data_matrix_log2" = diffTestRes$data_matrix_log2,
             "data_matrix_log2_imp" = diffTestRes$data_matrix_log2_imp)

  if ("differentialAbundanceTestMSResults$" %in% names(DIA_resultset)){
    stop("DIA MS resultset already contains differential abundance test results.
         Rename or remove these before running a new differential test.")
  } else{
    DIA_resultset$differentialAbundanceTestMSResults = res
    return(DIA_resultset)
  }

}
