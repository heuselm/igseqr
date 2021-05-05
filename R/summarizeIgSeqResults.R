#' @title summarizeIgSeqResults
#' @description summarize MiXCR IgSeq results as structured by https://github.com/heuselm/igseq_workflow.git
#' @param result_folder Location of the IgSeq-workflow result folder, not yet implemented (needs to be run IN that folder, for now)
#' @param write_tables Switch whether the collected summary information should be output as .csv tables
#' @import data.table ggplot2 ggrepel
#' @export

# library(data.table)
# library(ggplot2)
# library(ggrepel)

summarizeIgSeqResults <- function(result_folder = ".", write_tables = TRUE){

  ### collect & process / formalize results

  ## Step 1 - checkout
  step1_checkout = fread("01_checkout/checkout.log.txt")
  if ("SAMPLE" %in% names(step1_checkout)){
    samples = sort(step1_checkout$SAMPLE)
  } else {
    samples = sort(step1_checkout$V3[grep("S", step1_checkout$V3)])
  }

  ## Step 2 - histogram stats
  step2_histogram = fread("02_histogram/estimates.txt")
  step2_histogram_osq = fread("02_histogram/overseq.txt", header = T)

  ## Step 3 - assembly of UMI groups
  step3_assemble = fread("03_assemble/assemble.log.txt")

  ## Step 4 - Read Alignment to Germline Sequences
  step4_align = data.table()
  step4_align_files = list.files("04_aligned/")[grep("alignmentReport", list.files("04_aligned/"))]
  for (i in seq_along(step4_align_files)){
    dt.i = fread(paste0("04_aligned/", step4_align_files[i]), sep = ":", fill = T)[, SAMPLE:=gsub("alignmentReport_|.txt","",step4_align_files[i])]
    step4_align = rbind(step4_align, dt.i)
  }
  names(step4_align) = c("var", "value_merged", "V3", "V4", "SAMPLE")

  # massage to plottable data structure
  step4_align[, value_percent:=gsub("%)", "", unlist(strsplit(value_merged, split="\\("))[2]), .(var,value_merged)]
  step4_align[, value_abs:=gsub("%)", "", unlist(strsplit(value_merged, split="\\("))[1]), .(var,value_merged)]
  step4_align_num = step4_align[var %in% step4_align$var[7:30]]
  step4_align_num = step4_align_num[, .(SAMPLE, var, value_abs, value_percent)]
  step4_align_num[, value_percent := as.numeric(value_percent)]
  step4_align_num[, value_abs := as.numeric(value_abs)]

  ## Step 5 - clone assembly reports (part 1)
  step5_cloneAssemblyReports = data.table()
  step5_cloneAssemblyReports_files = list.files("05_finalclones/")[grep("assemblyReport", list.files("05_finalclones/"))]
  for (i in seq_along(step5_cloneAssemblyReports_files)){
    dt.i = fread(paste0("05_finalclones/", step5_cloneAssemblyReports_files[i]), sep = "X", header = F)[, SAMPLE:=gsub("cloneassemblyReport|.txt","",step5_cloneAssemblyReports_files[i])]
    dt.i[, var:=unlist(strsplit(V1, split = ": "))[1], V1]
    dt.i[, value_merged:=unlist(strsplit(V1, split = ": "))[2], V1]
    step5_cloneAssemblyReports = rbind(step5_cloneAssemblyReports, dt.i)
  }

  # massage to plottable data structure
  step5_cloneAssemblyReports[, value_percent:=gsub("%)", "", unlist(strsplit(value_merged, split="\\("))[2]), .(var,value_merged)]
  step5_cloneAssemblyReports[, value_abs:=gsub("%)", "", unlist(strsplit(value_merged, split="\\("))[1]), .(var,value_merged)]
  step5_cloneAssemblyReports_num = step5_cloneAssemblyReports[var %in% step5_cloneAssemblyReports$var[7:24]]
  step5_cloneAssemblyReports_num = step5_cloneAssemblyReports_num[, .(SAMPLE, var, value_abs, value_percent)]
  step5_cloneAssemblyReports_num[, value_percent := as.numeric(value_percent)]
  step5_cloneAssemblyReports_num[, value_abs := as.numeric(value_abs)]

  ## Step 5 - final clones (part 2)
  step5_finalclones = data.table()
  step5_finalclones_files = list.files("05_finalclones/")[grep("clones_all", list.files("05_finalclones/"))]
  for (i in seq_along(step5_finalclones_files)){
    dt.i = fread(paste0("05_finalclones/", step5_finalclones_files[i]))[, SAMPLE:=gsub("_clones_all.txt", "", step5_finalclones_files[i])]
    step5_finalclones = rbind(step5_finalclones, dt.i)
  }

  # Add histograms (or density plots) that display the distribution of how many times each clone was observed.
  head(step5_finalclones)
  head(step2_histogram)
  step2_histogram$SAMPLE = step2_histogram$`#SAMPLE_ID`
  step2_histogram[, overseq_ratio:=TOTAL_READS/TOTAL_MIGS]
  step5_finalclones = merge(step5_finalclones, step2_histogram[, .(SAMPLE, overseq_ratio)], by = "SAMPLE")
  step5_finalclones[, chain_type:=substr(allCHitsWithScore, 1, 3)]
  step5_finalclones[chain_type == "", chain_type:=substr(allVHitsWithScore, 1, 3)]

  # re-index of clones by chain type
  step5_finalclones[, cloneId_by_chain_type:=frank(-cloneCount,ties.method = "first"), .(SAMPLE, chain_type)]

  # Add Isotype information for the heavy chains
  step5_finalclones[, chain_isotype:="NA"]
  step5_finalclones[chain_type == "IGH", chain_isotype:=substr(allCHitsWithScore, 1, 4)]

  ### Visualize key info from each step in plot

  ## Step 1 - checkout
  names(step1_checkout)
  df = melt(step1_checkout[, .(SAMPLE, MASTER, `MASTER+SLAVE`)])

  ggplot(df, aes(SAMPLE, value, fill = variable)) + geom_bar(stat = "identity", position = position_dodge()) + theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("N assigned reads")
  ggsave("Step1_Nassigned_reads.pdf", height = 5, width = 5)

  ## Step 2 - overseq histogram
  head(step2_histogram)
  head(step2_histogram_osq)
  df2 = melt(step2_histogram_osq[, c(1,3:ncol(step2_histogram_osq)), with = F])

  ggplot(df2, aes(variable, value)) + geom_bar(stat = "identity", aes(fill = `#SAMPLE_ID`), position = position_dodge())+
    ggtitle("Oversequencing histogram per sample") + facet_wrap(~`#SAMPLE_ID`, ncol = 1)
  ggsave("Step2_Overseq_per_sample.pdf", height = 10)

  ## Step 3 - assemble
  head(step3_assemble)
  df3 = melt(step3_assemble)
  setnames(df3, '#SAMPLE_ID', 'SAMPLE')

  ggplot(df3, aes(SAMPLE, value)) + geom_bar(aes(fill = SAMPLE), stat = "identity") + facet_wrap(~variable, scales = "free") +
    ggtitle("Quality metrics from MIG assembly process") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("Step3_MIG_Consensus_assembly.pdf", width = 25, height = 14)

  ## Step 4 - alignment to germline sequences
  ggplot(step4_align_num) + geom_bar(aes(SAMPLE, value_abs, fill = SAMPLE), stat = "identity") + facet_wrap(~var, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Step4 alignment to germline sequences")
  ggsave("Step4_alignToGermline.pdf", width = 28, height = 14)

  # Step 5: Assembly/Clustering into final clones
  head(step5_cloneAssemblyReports)

  ggplot(step5_cloneAssemblyReports_num) + geom_bar(aes(SAMPLE, value_abs, fill = SAMPLE), stat = "identity") + facet_wrap(~var, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Step5 assemble clones")
  ggsave("Step5_assembleClones.pdf", width = 28, height = 14)

  ggplot(step5_cloneAssemblyReports_num[grep("^IG", var)]) + geom_bar(aes(SAMPLE, value_abs, fill = SAMPLE), stat = "identity") +
    geom_text(aes(SAMPLE, value_abs, label=value_abs), vjust = -0.2) + facet_wrap(~var) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Step5 assemble clones - IG types")
  ggsave("Step5_assembleClones_ChainTypes.pdf", width = 12, height = 6)

  names(step5_finalclones)
  ggplot(step5_finalclones) + geom_density(aes(nchar(aaSeqCDR3), col = SAMPLE), adjust = 4, size = 3, alpha = 0.5) + ggtitle("CDR3 length")
  ggsave("step5_finalClones_CDR3Length.pdf", height = 5, width = 5)

  ggplot(step5_finalclones) + geom_violin(aes(SAMPLE, y = minQualCDR3, fill = SAMPLE), bw = 1) + ggtitle("Minimum quality of CDR3 na Sequence region") +
    theme(axis.text.x = element_blank())
  ggsave("step5_finalClones_CDR3MinQual.pdf", height = 3, width = 5)

  # cloneCount per Sample (all chain types)
  ggplot(step5_finalclones, aes(cloneId, cloneCount)) +
    geom_line(aes(color = paste(SAMPLE, "osqR:", round(overseq_ratio, 2)))) + theme_bw() +
    ggtitle("LUIGSeq R01 cloneCount distribution, linear scale") +
    theme(legend.justification = c(1, 1.1), legend.position = c(1, 1))
  ggsave("step5_finalClones_cloneCount_lin.pdf")

  ggplot(step5_finalclones, aes(cloneId, cloneCount)) +
    geom_line(aes(color = paste(SAMPLE, "osqR:", round(overseq_ratio, 2))), size = 2) + theme_bw() +
    scale_y_log10() + ggtitle("LUIGSeq R01 cloneCount distribution, log scale") +
    theme(legend.justification = c(1, 1.1), legend.position = c(1, 1))
  ggsave("step5_finalClones_cloneCount_log.pdf", height = 7, width = 10)

  ggplot(step5_finalclones, aes(cloneId, cloneCount)) +
    geom_line(aes(color = paste(SAMPLE, "osqR:", round(overseq_ratio, 2))), size = 2) + theme_bw() +
    scale_y_log10() + facet_wrap(~SAMPLE) +
    ggtitle("LUIGSeq R01 cloneCount distribution, log scale")
  ggsave("step5_finalClones_cloneCount_log_split.pdf", height = 7, width = 10)

  # Clone count per sample split by chain
  ggplot(step5_finalclones, aes(cloneId_by_chain_type, cloneCount)) +
    geom_line(aes(color = paste(SAMPLE, "osqR:", round(overseq_ratio, 2)), linetype = chain_type), size = 1) +
    theme_bw() + scale_y_log10() + ggtitle("LUIGSeq R01 cloneCount distribution split by chain type, log scale") +
    theme(legend.justification = c(1, 1.1), legend.position = c(1, 1))
  ggsave("step5_finalClones_cloneCount_byChain_log.pdf", height = 7, width = 10)

  ggplot(step5_finalclones, aes(cloneId_by_chain_type, cloneCount)) +
    geom_line(aes(color = paste(SAMPLE, "osqR:", round(overseq_ratio, 2)), linetype = chain_type), size = 1) +
    theme_bw() + scale_y_log10() + facet_wrap(~SAMPLE) +
    ggtitle("LUIGSeq R01 cloneCount distribution, log scale")
  ggsave("step5_finalClones_cloneCount_byChain_log_split.pdf", height = 7, width = 14)

  ggplot(step5_finalclones, aes(cloneId_by_chain_type, cloneCount)) +
    geom_line(aes(color = paste(SAMPLE, "osqR:", round(overseq_ratio, 2)), linetype = chain_type), size = 1) +
    theme_bw() + scale_y_log10() + scale_x_log10() + facet_wrap(~SAMPLE) +
    ggtitle("LUIGSeq R01 cloneCount distribution, log scale")
  ggsave("step5_finalClones_cloneCount_byChain_logxy_split.pdf", height = 7, width = 14)

  # Count final number of sequenced chains
  ggplot(step5_finalclones[, length(unique(cloneId)), .(SAMPLE, chain_type)]) +
    geom_bar(aes(SAMPLE, V1, fill = SAMPLE, shape = chain_type), stat = "identity") +
    geom_text(aes(SAMPLE, V1, label = V1), vjust = -0.2) +
    facet_wrap(~chain_type) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("Number of unique chain sequences detected") +
    ggtitle("IGSeq number of chains sequenced")
  ggsave("step5_finalClones_NChains_sequenced_bar.pdf", height = 6, width = 12)

  ggplot(step5_finalclones[, length(unique(cloneId)), .(SAMPLE, chain_type)]) + geom_point(aes(SAMPLE, V1, color = SAMPLE, shape = chain_type), size = 3) + geom_text_repel(aes(SAMPLE, V1, label = paste(chain_type,"\n",V1)), nudge_x = 0.2) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Number of unique chain sequences detected") + ggtitle("IGSeq number of chains sequenced")
  ggsave("step5_finalClones_NChains_sequenced.pdf", height = 7, width = 7)

  # Plot isotype distribution
  ggplot(step5_finalclones[, length(unique(cloneId)), .(SAMPLE, chain_isotype)]) +
    geom_bar(aes(SAMPLE, V1, fill = chain_isotype), stat = "identity", position = "dodge") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Number of unique chain sequences detected") + ggtitle("IGSeq IGH chain isotypes observation")
  ggsave("step5_finalClones_IGH_isotype_obs.pdf", height = 7, width = 7)

  ## Write tables
  if (write_tables){
    fwrite(step1_checkout, paste0("summarizeIgSeq_", "step1_checkout.csv"))
    fwrite(step2_histogram, paste0("summarizeIgSeq_", "step2_histogram.csv"))
    fwrite(step3_assemble, paste0("summarizeIgSeq_", "step3_assemble.csv"))
    fwrite(step4_align, paste0("summarizeIgSeq_", "step4_align.csv"))
    fwrite(step5_cloneAssemblyReports_num, paste0("summarizeIgSeq_", "step5_cloneAssemblyReports_num.csv"))
    fwrite(step5_finalclones, paste0("summarizeIgSeq_", "step5_finalclones.csv"))
  }

}



