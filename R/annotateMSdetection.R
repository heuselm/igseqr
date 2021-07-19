#' @title annotateMSdetection
#'
#' @description annotate MS detection of IGSeq clones based on DIA_resultset
#'
#' @param IGSeq_resultset IgSeq resultset to be annotated with MS detection
#' @param DIA_resultset DIA MS resultset based on which to annotate MS detection of the clones in IGSeq_resultset.
#' @param write_table Whether to write .csv table of MS detection. Default: TRUE
#' @param n_prot_uniqueness_threshold Uniqueness threshold to define clones detected with sufficient specificity.
#' This may be useful to quantify clonal families with high homology with a certain degree of sharedness of the detectable peptide set.
#' Defaults to 1, i.e. only peptides mapping to only one of the clones in the context of the current .fasta will be maintained for
#'  unique analysis
#'
#' @return IGSeq_resultset with added element step7_clonesMSdetection with tables $clones_detected
#' and quantitative matrices across the MS runs in the DIA_resultset calculated from all mapped or only unique peptides
#' according to the n_prot_uniqueness threshold set as parameter:
#' $MS_abundance_allprecursors_mean, $MS_abundance_uniqueprecursors_mean
#'
#' @import data.table
#'
#' @export

annotateMSdetection = function(IGSeq_resultset = seq_r06_fasta,
                               DIA_resultset = ms_r06,

                               write_table = TRUE,
                               n_prot_uniqueness_threshold = 1){

  message("IGSeq-res: ", deparse(substitute(IGSeq_resultset)))
  message("DIA-res: ", deparse(substitute(DIA_resultset)))
  file_out_name = paste0("MSdetection_", deparse(substitute(IGSeq_resultset)), "_", deparse(substitute(DIA_resultset)), ".csv")

  # Create target object name depending on input
  res_list_name = paste0("step7_MS_detection_", deparse(substitute(DIA_resultset)))
  IGSeq_resultset[[res_list_name]] = list()

  # First, subset clone table to those detected in MS at all:
  ms_detected_all = unique(DIA_resultset$data_wide$Protein.Group)
  ms_detected_all_groupsresolved = unique(unlist(lapply(ms_detected_all, function(x){strsplit(x, split = ";")})))
  IGSeq_resultset[[res_list_name]]$clones_detected = IGSeq_resultset$step6_clonestofasta[cloneIdGlobal %in% ms_detected_all_groupsresolved]
  IGSeq_resultset[[res_list_name]]$clones_detected

  # Then, for each of these clones, produce the above metrics and store in appropriate containers
  msdetected_clones = IGSeq_resultset[[res_list_name]]$clones_detected$cloneIdGlobal

  # Target containers
  IGSeq_resultset[[res_list_name]]$MS_abundance_allprecursors_mean = matrix(nrow = length(msdetected_clones), ncol = ncol(DIA_resultset$data_wide)-2)
  IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean = matrix(nrow = length(msdetected_clones), ncol = ncol(DIA_resultset$data_wide)-2)
  rownames(IGSeq_resultset[[res_list_name]]$MS_abundance_allprecursors_mean) = msdetected_clones
  rownames(IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean) = msdetected_clones
  colnames(IGSeq_resultset[[res_list_name]]$MS_abundance_allprecursors_mean) = colnames(DIA_resultset$data_wide)[3:ncol(DIA_resultset$data_wide)]
  colnames(IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean) = colnames(DIA_resultset$data_wide)[3:ncol(DIA_resultset$data_wide)]

  # Run @TODO: Parallelize/refactor to speed up
  pb <- txtProgressBar(max = length(msdetected_clones), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)

  message("Retrieving MS detection and -quantification information for ", length(msdetected_clones),
          " clones across ", ncol(DIA_resultset$data_wide)-2, " MS runs")

  for (i in seq_along(msdetected_clones)){
    # print(i)
    setTxtProgressBar(pb, i)
    grep_pattern = paste0(msdetected_clones[i], "$|", msdetected_clones[i], ";")
    quant_wide = DIA_resultset$data_wide[grep(grep_pattern , Protein.Group)]
    quant_prec = as.matrix(quant_wide[,3:ncol(quant_wide)])
    # Protein group members and size
    n_proteins = sapply(unique(quant_wide$Protein.Group), FUN = function(x){length(unlist(strsplit(x, split = ";")))})
    proteins = unique(quant_wide$Protein.Group)

    # All precursors
    n_precursors_all = length(unique(quant_wide$Precursor.Id))
    n_peptides_all = length(unique(sapply(quant_wide$Precursor.Id, function(x){substr(x, 1, nchar(x)-1)})))
    precursors_all = paste(unique(quant_wide$Precursor.Id), collapse = ";")
    peptides_all = paste(unique(sapply(quant_wide$Precursor.Id, function(x){substr(x, 1, nchar(x)-1)})), collapse = ";")
    quant_prot_all = colSums(quant_prec, na.rm = TRUE)

    # Unique precursors only
    quant_wide_unique = quant_wide[n_proteins == n_prot_uniqueness_threshold]
    quant_prec_unique = as.matrix(quant_wide_unique[,3:ncol(quant_wide_unique)])
    n_precursors_unique = length(unique(quant_wide_unique$Precursor.Id))
    n_peptides_unique = length(unique(sapply(quant_wide_unique$Precursor.Id, function(x){substr(x, 1, nchar(x)-1)})))
    precursors_unique = paste(unique(quant_wide_unique$Precursor.Id), collapse = ";")
    peptides_unique = paste(unique(sapply(quant_wide_unique$Precursor.Id, function(x){substr(x, 1, nchar(x)-1)})), collapse = ";")
    quant_prot_unique = colSums(quant_prec_unique, na.rm = TRUE)

    # Write global stats to clones_detected table
    IGSeq_resultset[[res_list_name]]$clones_detected$MSdetection_protein_group_sizes[i] = paste(n_proteins, collapse = ";")
    IGSeq_resultset[[res_list_name]]$clones_detected$MSdetection_protein_groups[i] = paste(proteins, collapse = ":")
    IGSeq_resultset[[res_list_name]]$clones_detected$n_precursors_all[i] = n_precursors_all
    IGSeq_resultset[[res_list_name]]$clones_detected$n_peptides_all[i] = n_peptides_all
    IGSeq_resultset[[res_list_name]]$clones_detected$precursors_all[i] = precursors_all
    IGSeq_resultset[[res_list_name]]$clones_detected$peptides_all[i] = peptides_all
    IGSeq_resultset[[res_list_name]]$clones_detected$n_precursors_unique[i] = n_precursors_unique
    IGSeq_resultset[[res_list_name]]$clones_detected$n_peptides_unique[i] = n_peptides_unique
    IGSeq_resultset[[res_list_name]]$clones_detected$precursors_unique[i] = precursors_unique
    IGSeq_resultset[[res_list_name]]$clones_detected$peptides_unique[i] = peptides_unique

    # Write quantification to quant matrices
    IGSeq_resultset[[res_list_name]]$MS_abundance_allprecursors_mean[i,] = quant_prot_all
    IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean[i,] = quant_prot_unique

  }
  close(pb)

  ## aggregate by study design
  MS_abundance_uniqueprecursors_mean_long = reshape2::melt(IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean,
                                                           value.name = "MS_abundance_uniqueprecursors_mean")
  names(MS_abundance_uniqueprecursors_mean_long)[1:2] = c("cloneIdGlobal", "filename")
  MS_abundance_uniqueprecursors_mean_long = merge(MS_abundance_uniqueprecursors_mean_long, DIA_resultset$study_design,
                                                  by = "filename")
  MS_abundance_uniqueprecursors_mean_long = as.data.table(MS_abundance_uniqueprecursors_mean_long)
  MS_abundance_uniqueprecursors_mean_long[, nrep := 0]
  MS_abundance_uniqueprecursors_mean_long[MS_abundance_uniqueprecursors_mean > 0, nrep:=length(unique(replicate)), .(cloneIdGlobal,condition)]

  # Additional matrices + table:
  # MS_abundance per condition and replicate
  MS_abundance_uniqueprecursors_mean_cond_rep = dcast(MS_abundance_uniqueprecursors_mean_long, cloneIdGlobal~condition+replicate,
                                                      value.var = "MS_abundance_uniqueprecursors_mean")
  MS_abundance_uniqueprecursors_mean_cond_rep.m = as.matrix(MS_abundance_uniqueprecursors_mean_cond_rep[, 2:ncol(MS_abundance_uniqueprecursors_mean_cond_rep)])
  row.names(MS_abundance_uniqueprecursors_mean_cond_rep.m) = MS_abundance_uniqueprecursors_mean_cond_rep$cloneIdGlobal

  # MS abundance per condition
  MS_abundance_uniqueprecursors_mean_cond = dcast(MS_abundance_uniqueprecursors_mean_long, cloneIdGlobal~condition,
                                                  value.var = "MS_abundance_uniqueprecursors_mean", fun.aggregate = mean, fill = 0)
  MS_abundance_uniqueprecursors_mean_cond.m = as.matrix(MS_abundance_uniqueprecursors_mean_cond[, 2:ncol(MS_abundance_uniqueprecursors_mean_cond)])
  row.names(MS_abundance_uniqueprecursors_mean_cond.m) = MS_abundance_uniqueprecursors_mean_cond$cloneIdGlobal

  # MS detection (# of replicates) per condition
  MS_detection_uniqueprecursors_condition = dcast(MS_abundance_uniqueprecursors_mean_long, cloneIdGlobal~condition, value.var = "nrep",
                                                  fun.aggregate = unique, fill = 0)
  MS_detection_uniqueprecursors_condition.m = as.matrix(MS_detection_uniqueprecursors_condition[, 2:ncol(MS_detection_uniqueprecursors_condition)])
  row.names(MS_detection_uniqueprecursors_condition.m) = MS_detection_uniqueprecursors_condition$cloneIdGlobal

  # MS_abundance with these metrics in long format

  # write to result object:
  IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean_cond_rep = MS_abundance_uniqueprecursors_mean_cond_rep.m
  IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean_cond = MS_abundance_uniqueprecursors_mean_cond.m
  IGSeq_resultset[[res_list_name]]$MS_detection_uniqueprecursors_condition = MS_detection_uniqueprecursors_condition.m
  IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean_long = MS_abundance_uniqueprecursors_mean_long
  IGSeq_resultset[[res_list_name]]$study_design = DIA_resultset$study_design

  if (write_table){
    # Assemble & write combined table from MSannotation results
    clone_detection_and_quant = cbind(IGSeq_resultset[[res_list_name]]$clones_detected,
                                      IGSeq_resultset[[res_list_name]]$MS_abundance_allprecursors_mean,
                                      IGSeq_resultset[[res_list_name]]$MS_abundance_uniqueprecursors_mean)
    fwrite(clone_detection_and_quant,
           file = file_out_name)

    fwrite(MS_detection_uniqueprecursors_condition,
           file = gsub(".csv", "_cond_nrep.csv", file_out_name))

  }
  return(IGSeq_resultset)
}
