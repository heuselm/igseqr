#' @title differentialAbundanceTestMS
#'
#' @description test differential abundance of proteins, clones and clone groups in DIA-MS data on precursor level
#'
#' @param DIA_resultset DIA-MS resultset imported by importDIAresults function.
#' @param study_design_external Alternative study design overwriting the study design in DIA_resultset containing columns filename
#' @param comparison_matrix A matrix defining the comparisons to be run. column one contains the condition in the counter of each comparison, column two the denominator condition.
#' in other words, positive log2FCs will mean higher signal in the conditions listed in column 1. For possible values, check unique(DIA_resultset$study_design$condition). Default: TRUE.
#' If you would like to run all possible pairwise comparisons, enter "all_pairs". The data will be subset to the conditions specified in the comparison matrix.
#' @param write_intermediate_results Whether the intermediate results of each pairwise comparison shall be written. Separate subfolders will be generated. Default: FALSE
#' @param write_csv_protein Whether to write differential abundance testing result table summarized to protein level (csv). Default: TRUE
#' @param write_csv_precursor Whether to write differential abundance testing result table containing precursor level information (csv). Default: FALSE
#' @param plot_pdf Whether to plot summary of the differential tests.
#' @param plot_html Whether to plot interactive summary of the differential tests.
#'
#' @return DIA_resultset with added element $differentialAbundanceTestResults with sub-tables $experimental_design $comparison_matrix, $diffRes_pep and $diffRes_prot
#'
#' @import data.table ggplot2 ggrepel seqinr
#'
#' @export

differentialAbundanceTestMS = function(DIA_resultset,
                                       study_design_external = NULL,
                                       comparison_matrix = matrix(c("AP_m11_tF", "AP_m11_t0",
                                                                    "AP_m12_tF", "AP_m12_t0",
                                                                    "AP_m1_tF", "AP_m1_t0",
                                                                    "AP_m2_tF", "AP_m2_t0"),
                                                                  ncol = 2, byrow = TRUE),
                                       write_intermediate_results = TRUE,
                                       write_diffTable_csv_protein = TRUE,
                                       write_diffTable_csv_precursor = TRUE,
                                       plot_pdf = TRUE,
                                       plot_html = TRUE){


}
