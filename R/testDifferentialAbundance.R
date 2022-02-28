#' Test Differential Protein Expression in MS proteomics data
#'
#' @description Test Differential Protein Expression in MS proteomics data starting small: From the precursor level.
#'
#' @param input_dt Input data table either in tsv/txt format or already in R as data.table or data.frame with the following columns:
#' #' The table should have the following columns:
#' \itemize{
#' \item Protein.Group: Semicolon-separated Uniprot IDs (or similar, as long as it matches)
#' \item Precursor.Id: Unique Precursor Id for which the quantitative values are contained
#' \item "filename": The file names of the MS raw data, must be identical to the entries in study_design$filename
#' }
#' Note: The data will be log2-transformed internally.
#'
#' @param protein_group_annotation Protein annotation table with columns Protein.Group and Protein.Names (and others if desired)
#' that will be used to annotate the results. By default it is assumed to be a subset of and and an attempt will be made
#' to extract it from the input_dt.
#' @param study_design Study design in tab-separated .txt with mandatory columns:
#' \itemize{
#' \item filename: Must match quantitative data-containing column headers in the input_dt
#' \item condition: String, biological condition (e.g. "treated" and "untreated")
#' \item replicate: Replicate number (integer). Minimally 3 replicates are needed per condition for this type of analysis.
#' }
#' @param normalize_data Whether or not data is scaled/normalized before differential testing. In some cases
#' it might be preferable not to scale the datasets, e.g. when comparing pulldowns vs. input samples! Defaults to TRUE.
#' @param normalization_function Normalization function to use that transforms a matrix of quantities where columns are
#' samples and rows are analytes. Defaults to limma:normalizeQuantiles, but can be replaced with any such function. You may want
#' to try limma::normalizeVSN or limma::normalizeMedianValues.
#' @param condition_1 Manual override to the condition 1 for the differential comparison. By default it is guessed from unique(study_design$condition)
#' @param condition_2 Manual override to the condition 2 for the differential comparison. By default it is guessed from unique(study_design$condition)
#' @param min_n_obs Minimum number of observations per precursor (number of runs it was identified in)
#' in order to keep in in the analysis
#' @param imp_percentile Percentile of the total distribution of values on which the random
#' distribution for sampling will be centered
#' @param imp_sd standard deviation of the normal distribution from which values are sampled to impute missing values
#' @param plot_pdf Document processing steps in a string of pdf graphs
#' @param write_tsv_tables Write out final quant table with differential expression testing results
#' @param target_protein Optional string with protein identifier to highlight in volcano plots
#'
#' @return A diffExpr object (list) containing (access by x$ or by "x[[name]]")
#' \itemize{
#' \item data_source: input_dt path or input R object name
#' \item data_long: Data in long format
#' \item data_matrix_log2: Data, filtered and log2 transformed, in wide format matrix
#' \item data_matrix_log2_imp: Data, filtered, log2 transformed and with missing values imputed, in wide format matrix
#' \item study_design: study design table
#' \item annotation_col: column annotation
#' \item diffExpr_result_dt: Result table with intensities and differential expression testing results
#' \item candidates_condition1: Proteins that appear higher abundant in condition 1
#' \item candidates_condition2: Proteins that appear higher abundant in condition 2
#' \item
#' }
#'
#' @author Moritz Heusel
#'
#' @import data.table ggplot2 corrplot pheatmap ggrepel limma
#' @importFrom reshape2 melt
#'
#' @export

testDifferentialAbundance <- function(input_dt = "path/to/DIANN_matrix.tsv",
                                      protein_group_annotation = NULL,
                                      study_design = "path/to/Study_design_filled.tsv",

                                      # toggle normalization & -function
                                      normalize_data = TRUE,
                                      normalization_function = limma::normalizeQuantiles,
                                      # also try

                                      # select conditions to be compared
                                      condition_1 = unique(fread(study_design)$condition)[2],
                                      condition_2 = unique(fread(study_design)$condition)[1],

                                      # filtering options
                                      min_n_obs = 4,

                                      # imputation of missing values options
                                      imp_percentile = 0.001,
                                      imp_sd = 0.2,

                                      # output options
                                      plot_pdf = TRUE,
                                      write_tsv_tables = TRUE,

                                      # target protein highlight
                                      target_protein = "O08760")
{
  #####################################################################################################

  ## Set seed to ensure reproducibility
  set.seed(123)

  ## Processing
  # load data
  if (is.character(input_dt)){
    data = fread(input_dt)
  }else {
    data = as.data.table(input_dt)
  }

  if (is.character(study_design)){
    annotation = fread(study_design)
  }else {
    annotation = as.data.table(study_design)
  }

  if (is.null(protein_group_annotation)){
    # Try to retrieve protein annotation from input table
    if("Protein.Names" %in% names(data)){
      protein_group_annotation = unique(data[, .(Protein.Group, Protein.Names)])
    } else {
      data[, Protein.Names:=Protein.Group]
      protein_group_annotation = unique(data[, .(Protein.Group, Protein.Names)])
      warning("Protein.Names column is missing in input data\n
      Protein.Group identifiers will be used, unless Protein.Group-to-Protein.Name mapping is provided via
      protein_group_annotation table")
    }
  }

  # Stop if annotation in study design is insufficient
  if (!(all(c("filename","condition","replicate") %in% names(annotation)))){
    stop("study design must contain columns filename,condition and replicate")
  }

  # subset to runs in annotation, allowing some flexibility for colname structure
  if (sum(grep("\\\\", names(data)))>0)
    {
    names_split = sapply(names(data), function(x){unlist(strsplit(x, split = "\\\\"))[length(unlist(strsplit(x, split = "\\\\")))]})
    if (sum(grep(".raw", names(data)))>0){
      names_split = gsub(".raw", "", names_split)
    }
    names(data) = names_split
  }

  # Now get only the relevant columns
  data.s.wide = cbind(data[, .(Protein.Group, Precursor.Id)],
                      data[, which(names(data) %in% annotation$filename), with = F])

  if(ncol(data.s.wide)<4){
    stop("filename in study design must match column name in input_dt!")
  }
  annotation_col = data.table("cn" = names(data.s.wide)[3:ncol(data.s.wide)])
  annotation_col = merge(annotation_col, annotation, by.x = "cn", by.y = "filename")
  annotation_col[, new_cn:=paste(cn,condition,replicate)]
  names(data.s.wide)[3:ncol(data.s.wide)] = annotation_col$new_cn
  annotation_col = as.data.frame(annotation_col)
  row.names(annotation_col) = annotation_col$new_cn
  annotation_col$cn = NULL
  annotation_col$new_cn = NULL

  # Convert to matrix format
  data.s.wide.quant = as.matrix(data.s.wide[,3:ncol(data.s.wide)])
  row.names(data.s.wide.quant) = data.s.wide$Precursor.Id

  # log transform and remove inf
  data.s.wide.quant.log2 = log2(data.s.wide.quant)
  data.s.wide.quant.log2[is.infinite(data.s.wide.quant.log2)] <- 0

  # Plot abundance distributions before and after normalization (if applied)
  boxplot(data.s.wide.quant.log2, las = 2, main = "Input intensities, log2-transformed")

  if (normalize_data == FALSE) {
    normalization_function = function(x){return(x)}
  }

  data.s.wide.quant.log2.qnorm = normalization_function(data.s.wide.quant.log2)
  boxplot(data.s.wide.quant.log2.qnorm, las = 2, main = "log2-transformed and quantile-normalized intensities")
  data.s.wide.log2.qnorm = cbind(data.s.wide[,1:2],data.s.wide.quant.log2.qnorm)

  # Plot normalization effect in density plot
  data.s.long.prenorm = reshape2::melt(data.s.wide.quant.log2)
  data.s.long.prenorm$normalization = "before (input)"
  data.s.long.qNorm = reshape2::melt(data.s.wide.quant.log2.qnorm)
  data.s.long.qNorm$normalization = "after normalization"
  data.s.long = rbind(data.s.long.prenorm, data.s.long.qNorm)

  ggplot(data.s.long, aes(value, group = Var2, color = normalization)) + geom_density() + facet_wrap(~normalization, ncol = 1) +
    ggtitle(paste("Impact of normalization, normalize_data = ",normalize_data)) + theme_bw() + xlab("log2 Precursor.Quantity")
  if(plot_pdf){
    ggsave("Normalization_impact_abundance_density.pdf", height = 5, width = 5)
  }

  ## Impute missing values
  # Visualize input data before imputation
  data.s.wide.quant.log2.qnorm.noNa = copy(data.s.wide.quant.log2.qnorm)
  data.s.wide.quant.log2.qnorm.noNa[is.na(data.s.wide.quant.log2.qnorm.noNa)] = 0
  pheatmap::pheatmap(data.s.wide.quant.log2.qnorm.noNa[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa)),], show_rownames = F,
           main = "Quantile-normalized log2 precursor-level intensities, non-imputed", cluster_rows = F)
  if(plot_pdf){
    dev.copy(pdf, "Heatmap_unfiltered_pre_imputation.pdf", width = ncol(data.s.wide.quant.log2.qnorm.noNa)/2, height = 8)
    dev.off()
  }
  # too big for hclust, thus cluster_rows = F

  # Check correlation to remove outliers / select data subset
  dev.off()
  corrplot::corrplot(cor(data.s.wide.quant.log2.qnorm.noNa),
           method = "color",
           is.corr = F,
           order = "hclust",
           mar = c(4,4,4,4),
           tl.cex = 1-(ncol(data.s.wide.quant.log2.qnorm.noNa)/100),
           title = "Sample-sample correlation before min_obs filtering")
  if(plot_pdf){
    dev.copy(pdf, "Correlation_sample_sample_unfiltered.pdf")
    dev.off()
  }

  # Filter dataset based on minimum n_obs to avoid imputation around noise..
  n_obs = apply(data.s.wide.quant.log2.qnorm, MARGIN = 1, function(x) length(x)-sum(is.na(x)))
  hist(n_obs, n = 100)
  abline(v = min_n_obs, col = "red", lty = 2)
  data.s.wide.quant.log2.qnorm.noNa.minObs = data.s.wide.quant.log2.qnorm.noNa[n_obs >= min_n_obs,]
  nrow(data.s.wide.quant.log2.qnorm.noNa.minObs)

  plot(n_obs[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa))], main = "n_observations over abundance rank")
  abline(h = min_n_obs, col = "red", lty = 2)

  if(plot_pdf){
    par(mfrow = c(2,1))
    pdf("Filtering_min_obs.pdf")
    hist(n_obs, n = 100)
    abline(v = min_n_obs, col = "red", lty = 2)
    plot(n_obs[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa))],
         main = "n_observations over abundance rank")
    abline(h = min_n_obs, col = "red", lty = 2)
    dev.off()
  }

  # Plot heatmap after min obs filtering, before imputation
  pheatmap::pheatmap(data.s.wide.quant.log2.qnorm.noNa.minObs[order(-rowSums(data.s.wide.quant.log2.qnorm.noNa.minObs)),], show_rownames = F,
           ann_col = annotation_col,
           main = "Quantile-normalized log2 precursor-level intensities, min obs. filtered, before imputation",
           cluster_rows = F,
           cluster_cols = F)
  if(plot_pdf){
    dev.copy(pdf, "Heatmap_min_obs_filtered_pre_imputation.pdf", width = ncol(data.s.wide.quant.log2.qnorm.noNa)/2, height = 8)
    dev.off()
  }

  # Sample-sample correlation after filtering to min obs
  dev.off()
  corrplot::corrplot(cor(data.s.wide.quant.log2.qnorm.noNa.minObs),
           method = "color",
           is.corr = F,
           order = "hclust",
           mar = c(4,4,4,4),
           tl.cex = 1-(ncol(data.s.wide.quant.log2.qnorm.noNa)/100),
           title = "Sample-sample correlation after min_obs filtering")
  if(plot_pdf){
    dev.copy(pdf, "Correlation_sample_sample_min_obs_filtered_pre_imputation.pdf")
    dev.off()
  }

    ## Impute missing values
  # distribution of all values before imputation
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs)

  # impute missing values
  imp_percentile = quantile(unique(data.s.wide.quant.log2.qnorm.noNa.minObs[data.s.wide.quant.log2.qnorm.noNa.minObs>0]),
                        na.rm = T, imp_percentile)
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp = copy(data.s.wide.quant.log2.qnorm.noNa.minObs)
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs.imp == 0] =
    rnorm(length(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs.imp == 0]),
          mean = imp_percentile,  sd = imp_sd)

  # visualize distribution of imputed values
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs)
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp)
  hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")

  if(plot_pdf){
    par(mfrow = c(2,1))
    pdf("Imputation_histograms.pdf")
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs)
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp)
    hist(data.s.wide.quant.log2.qnorm.noNa.minObs.imp[data.s.wide.quant.log2.qnorm.noNa.minObs == 0], add = T, col = "red")
    dev.off()
  }

  # heatmap after filtering to min obs and imputation
  visualize_imputation_heatmap = cbind(data.s.wide.quant.log2.qnorm.noNa.minObs, data.s.wide.quant.log2.qnorm.noNa.minObs.imp)
  colnames(visualize_imputation_heatmap)[(ncol(data.s.wide.quant.log2.qnorm.noNa.minObs)+1):ncol(visualize_imputation_heatmap)] =
    paste(colnames(data.s.wide.quant.log2.qnorm.noNa.minObs.imp), "imputed")

  # To avoid ram overflow the plottable number of precursors is limited to 20000
  if(nrow(visualize_imputation_heatmap)>20000){
    visualize_imputation_heatmap =
      visualize_imputation_heatmap[sample(1:nrow(visualize_imputation_heatmap), 20000, replace = F),]
  }

  pheatmap::pheatmap(visualize_imputation_heatmap[order(-rowSums(visualize_imputation_heatmap)),],
           main = "Quantile-normalized log2 precursor-level intensities, min obs.\nnon-imputed      vs       imputed",
           show_rownames = F,
           gaps_col = ncol(data.s.wide.quant.log2.qnorm.noNa.minObs.imp),
           cluster_rows = F,
           cluster_cols = F)
  dev.copy(pdf, "Imputation_heatmap.pdf", width = ncol(data.s.wide.quant.log2.qnorm.noNa.minObs.imp))
  dev.off()

  ## Precursor level differential Abundance testing t.tests

  # Assemble differential expression testing result table
  res = copy(data.s.wide[,1:2])

  # Add raw data to the result table
  data.s.wide.quant.dt = as.data.table(data.s.wide.quant, keep.rownames = T)
  names(data.s.wide.quant.dt)[-1] = paste(names(data.s.wide.quant.dt)[-1], "raw_intensities")
  res = merge(res, data.s.wide.quant.dt,
              by.x = "Precursor.Id", by.y = "rn", all.x = F)

  # Add the filtered, normalized, imputed quantdata
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt = as.data.table(data.s.wide.quant.log2.qnorm.noNa.minObs.imp, keep.rownames = T)
  names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt)[-1] = paste(names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt)[-1], "log2_qnorm_imp_intensities")

  res = merge(res, data.s.wide.quant.log2.qnorm.noNa.minObs.imp.dt,
              by.x = "Precursor.Id", by.y = "rn", all.x = F)

  ## perform t.tests

  # get conditions
  if (is.null(condition_1)){
    condition_1 = unique(annotation$condition[1])
  }
  condition_1
  if (is.null(condition_1)){
    condition_2 = unique(annotation$condition[2])
  }
  condition_2

  # Keep only precursors that made it through the filtering
  res = res[Precursor.Id %in% row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp)]
  res[, comparison:=paste0(condition_1,"/",condition_2)]

  # preallocate result vectors
  prec = res$Precursor.Id
  n_prec = length(prec)
  pvals = numeric(length = n_prec)
  log2fcs = numeric(length = n_prec)

  # presplit matrix
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp =
    data.s.wide.quant.log2.qnorm.noNa.minObs.imp[order(row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp), method = "radix"),]

  # double-check ordering
  if (!(all(row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp) == res$Precursor.Id))){
    setorder(res, "Precursor.Id") # try to fix ordering
    if(!(all(row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp) == res$Precursor.Id))){
      stop("Index error, check ordering of precursors")
    }
  }
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond1 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp[, grep(condition_1, colnames(data.s.wide.quant.log2.qnorm.noNa.minObs.imp))]
  data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond2 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp[, grep(condition_2, colnames(data.s.wide.quant.log2.qnorm.noNa.minObs.imp))]

  # Start progress bar
  total = n_prec
  pb <- txtProgressBar(title = "Calculating precursor-level statistics", min = 0,
                       max = total, width = 300, style = 3)
  for (i in seq_along(prec)){
    setTxtProgressBar(pb, i)
    if (row.names(data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond1)[i] != prec[i]){
      stop("precursor mismatch")
    }
    v1 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond1[i,]
    v2 = data.s.wide.quant.log2.qnorm.noNa.minObs.imp.cond2[i,]

    if(all(v1 == v2)){
      pvals[i] = 1
      log2fcs[i] = 0
    } else{
      t_test_result =t.test(v1, v2)
      pvals[i] = t_test_result$p.value
      log2fcs[i] = mean(v1) - mean(v2)
    }
  }
  close(pb)
  res$p_value = pvals
  res$log2_fold_change = log2fcs

  # Visualize intermediate, precursor-level results in volcano plot
  res[, target_prot:=grepl(target_protein, Protein.Group), Precursor.Id]

  res_v = ggplot(res, aes(x = log2_fold_change, y = -log10(p_value), col = target_prot)) + geom_point() +
    scale_color_manual(values = c("darkgrey", "red")) +
    theme_bw() + ggtitle(paste0(condition_1, "/", condition_2, ", precursor-level statistics")) +
    geom_hline(yintercept = -log10(0.01), lty = 2) + geom_hline(yintercept = 0) +
    geom_vline(xintercept = c(-log2(2), log2(2)), lty = 2) + geom_vline(xintercept = 0)
  res_v
  if(plot_pdf){
    ggsave("Volcano_plot_precursorlevel.pdf", plot = res_v, height = 6, width = 6)
  }

  # Summarize to protein level
  res[, log2_fold_change_protein:=mean(log2_fold_change), Protein.Group]
  res[, n_precursors:=length(unique(Precursor.Id)), Protein.Group]
  res[, p_value_protein:=mean(p_value), Protein.Group]

  # do multiple testing correction
  pcorr = unique(res[, .(Protein.Group, p_value_protein)])
  pcorr[, p_value_BHadj_protein:=p.adjust(p_value_protein, method = "BH")]
  # Add corrected pvalues
  res = merge(res, pcorr[, .(Protein.Group, p_value_BHadj_protein)], by = "Protein.Group")

  # Add protein annotation from protein_group_annotation
  res = merge(res, protein_group_annotation, by = "Protein.Group")

  res_v_p = ggplot(res, aes(x = log2_fold_change_protein, y = -log10(p_value_BHadj_protein)))+
    geom_point() +
    scale_color_manual(values = c("darkgrey", "red")) +
    theme_bw() + ggtitle(paste0(condition_1, "/", condition_2, ", protein-level statistics")) +
    geom_hline(yintercept = -log10(0.01), lty = 2) +
    geom_vline(xintercept = c(-log2(2), log2(2)), lty = 2) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)
  res_v_p
  if(plot_pdf){
    ggsave("Volcano_plot_proteinlevel.pdf", plot = res_v_p, height = 6, width = 6)
  }

  # Plot labels for up- and down-regulated proteins
  res_v_p + geom_text_repel(data = res[p_value_BHadj_protein<=0.01 & log2_fold_change >= 1],
                            aes(label = Protein.Names), col = "red") +
    ggtitle(paste("Proteins higher abundant in condition_1: ", condition_1))
  ggsave(paste0("Volcano_plot_candidates_up_in_",condition_1,".pdf"))

  res_v_p + geom_text_repel(data = res[p_value_BHadj_protein<=0.01 & log2_fold_change <= -1],
                            aes(label = Protein.Names), col = "blue") +
    ggtitle(paste("Proteins higher abundant in condition_2: ", condition_2))
  ggsave(paste0("Volcano_plot_candidates_up_in_",condition_2,".pdf"))

  # Write out final result table
  write.table(res, "DiffTest_result.tsv", sep = "\t", quote = F, row.names = F)

  # return results
  names(data.s.long) = c("Precursor.Id", "filename", "Precursor.Quantity.log2", "normalization")
  data.s.long = merge(data.s.long, unique(data[, .(Protein.Group, Precursor.Id)]), by = "Precursor.Id", all.x = T)
  annotation[, colname_in_matrix:=paste(filename,condition,replicate)]

  res = list("data_source" = if(is.character(input_dt)){input_dt}else{deparse(substitute(input_dt))},
             "data_long" = data.s.long,
             "data_matrix_log2" = data.s.wide.quant.log2.qnorm.noNa.minObs,
             "data_matrix_log2_imp" = data.s.wide.quant.log2.qnorm.noNa.minObs.imp,
             "study_design" = annotation,
             "annotation_col" = annotation_col,
             "diffExpr_result_dt" = res,
             "candidates_condition1" = res[p_value_BHadj_protein<=0.01 & log2_fold_change >= 1],
             "candidates_condition2" = res[p_value_BHadj_protein<=0.01 & log2_fold_change <= -1])
  return(res)
}

