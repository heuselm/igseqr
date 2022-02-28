#' plotDifferentialAbundanceOverview
#'
#' @description Produce overview plots of multiple differential test results
#'
#' @param diffExpr_result_dt Input data table, typically appended result tables with differential expression testing results as produced by testDifferentialAbundance
#' @param significance_threshold_p Significance threshold on p-value to be indicated in plots (p_value_BHadj_protein <= x)
#' @param significance_threshold_fc Significance threshold on fold-change to be indicated in plots (|log2FC| >= log2(x)). A two-fold change is 2. The log2FC cutoff is then 1.
#' @param label_prefix Prefix label for the output plots. Defaults to "(input R object name)_".
#' @param remove_from_protein_names Regular expression of string components to remove from Protein.Names.
#' @param browsable_html Whether brosable volcanoe plot .htmls are to be generated. Default: TRUE.
#' @param protein_highlight_tag Substring of Protein names that get highlighted (color) in plots.
#' @param heatmap Whether an overview heatmap of FCs is to be drawn for the proteins passing the p-value cutoff at least once.
#'
#' @return A (list) containing (access by x$ or by "x[[name]]")
#' \itemize{
#' \item diffExpr_result_dt: Reduced diffExpr_result_dt processed for plots
#' \item volcanos: Volcano plots as ggplot object
#' \item heatmap: heatmap obj
#' \item FC_matrix: matrix of Log2FCs of proteins passing the p-value
#'  significance threshold in any one of the comparisons
#' }
#'
#' @author Moritz Heusel
#'
#' @import data.table ggplot2 pheatmap plotly crosstalk htmlwidgets
#'
#' @export

plotDifferentialAbundanceOverview <- function(diffExpr_result_dt,
                                      significance_threshold_p = 0.05,
                                      significance_threshold_fc = 2,
                                      label_prefix = NULL,
                                      remove_from_protein_names = "\\(Bos;|_MOUSE|_HUMAN",
                                      browsable_html = TRUE,
                                      protein_highlight_tag = "IGSeq",
                                      heatmap = TRUE){

  # Compare diffTest results in volcanos and heatmaps

  # collect protein-level results from differential table(s)
  cRes = unique(diffExpr_result_dt[, .(Protein.Group, Protein.Names, n_precursors,
                                       comparison, log2_fold_change_protein,
                                       p_value_protein, p_value_BHadj_protein)])
  # Highlight target_proteins
  cRes[, target_prot:=grepl(protein_highlight_tag, Protein.Group), Protein.Group]

  # get number of comparisons (volcano panels)
  ncomp = length(unique(cRes$comparison))

  # Get obj name as dataset label if not provided
  if(is.null(label_prefix)){
    label_prefix = deparse(substitute(diffExpr_result_dt))
  }

  # refine protein identifiers
  cRes[, n_proteins:=length(unlist(strsplit(Protein.Group, split = ";"))),Protein.Group]
  cRes[, Protein.Names.Short:=gsub(remove_from_protein_names, "", Protein.Names)]

  ## Volcano plots
  # Standard plot
  volcanos = ggplot(cRes, aes(log2_fold_change_protein, -log10(p_value_BHadj_protein),
                              text = Protein.Names.Short, col = target_prot)) +
    scale_color_manual(values = c("darkgrey", "red")) +
    geom_point(pch = ".") +
    facet_wrap(~comparison, nrow = 2) +
    theme_bw() +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(significance_threshold_p), col = "red", lty = 2) +
    geom_vline(xintercept = c(-log2(significance_threshold_fc), log2(significance_threshold_fc)), col = "red", lty = 2) +
    ggtitle(label_prefix)
  plot(volcanos)

  ggsave(paste0(label_prefix,".pdf"), height = ncomp*2+4, width = ncomp*2+4)

  if (browsable_html){
    # Save interactive html, standard
    htmlwidgets::saveWidget(ggplotly(volcanos), file = paste0(label_prefix,"_interactive_v1.html"))

    # Smarter plot with cross-panel highlighting
    volcanos_smart = ggplot(SharedData$new(cRes, key = ~Protein.Names.Short),
                            aes(log2_fold_change_protein, -log10(p_value_BHadj_protein),
                                text = Protein.Names.Short, col = target_prot)) +
      scale_color_manual(values = c("darkgrey", "red")) +
      geom_point(pch = ".") +
      facet_wrap(~comparison) +
      theme_bw() +
      theme(legend.position = "none") +
      geom_hline(yintercept = -log10(significance_threshold_p), col = "red", lty = 2) +
      geom_vline(xintercept = c(-log2(significance_threshold_fc), log2(significance_threshold_fc)),
                 col = "red", lty = 2) +
      ggtitle(label_prefix)

    htmlwidgets::saveWidget(ggplotly(volcanos_smart), file = paste0(label_prefix,"_interactive_v2.html"))
  }

  ## heatmap of Fold-change profiles
  diffprots = cRes[abs(log2_fold_change_protein)>=log2(significance_threshold_fc) &
                       p_value_BHadj_protein <= significance_threshold_p, unique(Protein.Names.Short)]
  diffprots.fcs = dcast(cRes[Protein.Names.Short %in% diffprots],
                        Protein.Names.Short~comparison, value.var = "log2_fold_change_protein")
  diffprots.fcs.m = as.matrix(diffprots.fcs[, 2:ncol(diffprots.fcs), with = F])
  row.names(diffprots.fcs.m) = diffprots.fcs$Protein.Names.Short
  diffprots.fcs.m[is.na(diffprots.fcs.m)] = 0

  hm = pheatmap(diffprots.fcs.m,
                cluster_cols = F,
                show_rownames = F,
                main = "differential proteins' global response profiles",
                # cutree_rows = 9,
                color = colorRampPalette(c("blue","white","red"))(100))
  dev.copy(pdf, paste0(label_prefix,"_heatmap_diffProtFcProfiles.pdf"), width = 4+ncol(diffprots.fcs.m)/10, height = 4+nrow(diffprots.fcs.m)/100)
  dev.off()

  res = list("diffExpr_result_dt" = cRes,
             "volcanos" = volcanos,
             "heatmap" = hm,
             "FC_matrix" = diffprots.fcs.m)
  return(res)
}
