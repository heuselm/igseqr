#' @title summarizeIgSeqResults
#' @description summarize MiXCR IgSeq results as structured by https://github.com/heuselm/igseq_workflow.git
#' @param IGSeq_resultset IgSeq resultset imported by importIgSeq function.
#' @return creates ggplot graphs in out_folder
#' @import data.table ggplot2 ggrepel
#' @export

summarizeIgSeq = function(IGSeq_resultset, report_folder = "summarizeIgSeq_out"){

  # Create output folder if it does not exist yet
  if (!dir.exists(report_folder)){dir.create(report_folder)} else{
    warning("report_folder ", report_folder, " already exists, contents are potentially being overwritten")}

  ### Visualize key info from each step in plot
  ## Step 1 - checkout
  names(IGSeq_resultset$step1_checkout)
  df = melt(IGSeq_resultset$step1_checkout[, .(SAMPLE_ID, MASTER, `MASTER+SLAVE`)])

  ggplot(df, aes(SAMPLE_ID, value, fill = variable)) + geom_bar(stat = "identity", position = position_dodge()) + theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("N assigned reads")
  ggsave(paste0(report_folder,"/","Step1_Nassigned_reads.pdf"), height = 5, width = 5)

  ## Step 2 - overseq histogram
  head(IGSeq_resultset$step2_histogram)
  head(IGSeq_resultset$step2_histogram_osq)
  df2 = melt(IGSeq_resultset$step2_histogram_osq[, c(1,3:ncol(IGSeq_resultset$step2_histogram_osq)), with = F])

  ggplot(df2, aes(variable, value)) + geom_bar(stat = "identity", aes(fill = SAMPLE_ID), position = position_dodge())+
    ggtitle("Oversequencing histogram per sample") + facet_wrap(~SAMPLE_ID, ncol = 1)
  ggsave(paste0(report_folder,"/","Step2_Overseq_per_sample.pdf"), height = 10)

  ## Step 3 - assemble
  head(IGSeq_resultset$step3_assemble)
  df3 = melt(IGSeq_resultset$step3_assemble)

  ggplot(df3, aes(SAMPLE_ID, value)) + geom_bar(aes(fill = SAMPLE_ID), stat = "identity") + facet_wrap(~variable, scales = "free") +
    ggtitle("Quality metrics from MIG assembly process") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(report_folder,"/","Step3_MIG_Consensus_assembly.pdf"), width = 25, height = 14)

  ## Step 4 - alignment to germline sequences
  ggplot(IGSeq_resultset$step4_align_log) + geom_bar(aes(SAMPLE_ID, value_abs, fill = SAMPLE_ID), stat = "identity") + facet_wrap(~var, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Step4 alignment to germline sequences")
  ggsave(paste0(report_folder,"/","Step4_alignToGermline.pdf"), width = 28, height = 14)

  # Step 5: Assembly/Clustering into final clones
  head(IGSeq_resultset$step5_cloneAssembly_log)

  ggplot(IGSeq_resultset$step5_cloneAssembly_log) + geom_bar(aes(SAMPLE_ID, value_abs, fill = SAMPLE_ID), stat = "identity") +
    facet_wrap(~var, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Step5 assemble clones")
  ggsave(paste0(report_folder,"/","Step5_assembleClones.pdf"), width = 28, height = 14)

  ggplot(IGSeq_resultset$step5_cloneAssembly_log[grep("^IG", var)]) + geom_bar(aes(SAMPLE_ID, value_abs, fill = SAMPLE_ID), stat = "identity") +
    geom_text(aes(SAMPLE_ID, value_abs, label=value_abs), vjust = -0.2) + facet_wrap(~var) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Step5 assemble clones - IG types")
  ggsave(paste0(report_folder,"/","Step5_assembleClones_ChainTypes.pdf"), width = 12, height = 6)

  names(IGSeq_resultset$step5_finalclones)
  ggplot(IGSeq_resultset$step5_finalclones) + geom_density(aes(nchar(aaSeqCDR3), col = SAMPLE_ID), adjust = 4, size = 3, alpha = 0.5) + ggtitle("CDR3 length")
  ggsave(paste0(report_folder,"/","step5_finalClones_CDR3Length.pdf"), height = 5, width = 5)

  ggplot(IGSeq_resultset$step5_finalclones) + geom_violin(aes(SAMPLE_ID, y = minQualCDR3, fill = SAMPLE_ID), bw = 1) + ggtitle("Minimum quality of CDR3 na Sequence region") +
    theme(axis.text.x = element_blank())
  ggsave(paste0(report_folder,"/","step5_finalClones_CDR3MinQual.pdf"), height = 3, width = 5)

  # cloneCount per Sample (all chain types)
  ggplot(IGSeq_resultset$step5_finalclones, aes(cloneId, cloneCount)) +
    geom_line(aes(color = SAMPLE_ID)) + theme_bw() +
    ggtitle("cloneCount distribution, linear scale") +
    theme(legend.justification = c(1, 1.1), legend.position = c(1, 1))
  ggsave(paste0(report_folder,"/","step5_finalClones_cloneCount_lin.pdf"), height = 7, width = 10)

  ggplot(IGSeq_resultset$step5_finalclones, aes(cloneId, cloneCount)) +
    geom_line(aes(color = SAMPLE_ID), size = 2) + theme_bw() +
    scale_y_log10() + ggtitle("cloneCount distribution, log scale") +
    theme(legend.justification = c(1, 1.1), legend.position = c(1, 1))
  ggsave(paste0(report_folder,"/","step5_finalClones_cloneCount_log.pdf"), height = 7, width = 10)

  ggplot(IGSeq_resultset$step5_finalclones, aes(cloneId, cloneCount)) +
    geom_line(aes(color = SAMPLE_ID), size = 2) + theme_bw() +
    scale_y_log10() + facet_wrap(~SAMPLE_ID) +
    ggtitle("cloneCount distribution, log scale")
  ggsave(paste0(report_folder,"/","step5_finalClones_cloneCount_log_split.pdf"), height = 7, width = 10)

  # Clone count per sample split by chain
  ggplot(IGSeq_resultset$step5_finalclones, aes(cloneId_by_chain_type, cloneCount)) +
    geom_line(aes(color = SAMPLE_ID, linetype = chain_type), size = 1) +
    theme_bw() + scale_y_log10() + ggtitle("cloneCount distribution split by chain type, log scale") +
    theme(legend.justification = c(1, 1.1), legend.position = c(1, 1))
  ggsave(paste0(report_folder,"/","step5_finalClones_cloneCount_byChain_log.pdf"), height = 7, width = 10)

  ggplot(IGSeq_resultset$step5_finalclones, aes(cloneId_by_chain_type, cloneCount)) +
    geom_line(aes(color = SAMPLE_ID, linetype = chain_type), size = 1) +
    theme_bw() + scale_y_log10() + facet_wrap(~SAMPLE_ID) +
    ggtitle("cloneCount distribution, log scale")
  ggsave(paste0(report_folder,"/","step5_finalClones_cloneCount_byChain_log_split.pdf"), height = 7, width = 14)

  ggplot(IGSeq_resultset$step5_finalclones, aes(cloneId_by_chain_type, cloneCount)) +
    geom_line(aes(color = SAMPLE_ID, linetype = chain_type), size = 1) +
    theme_bw() + scale_y_log10() + scale_x_log10() + facet_wrap(~SAMPLE_ID) +
    ggtitle("cloneCount distribution, log scale")
  ggsave(paste0(report_folder,"/","step5_finalClones_cloneCount_byChain_logxy_split.pdf"), height = 7, width = 14)

  # Count final number of sequenced chains
  ggplot(IGSeq_resultset$step5_finalclones[, length(unique(cloneId)), .(SAMPLE_ID, chain_type)]) +
    geom_bar(aes(SAMPLE_ID, V1, fill = SAMPLE_ID, col = chain_type), stat = "identity") +
    geom_text(aes(SAMPLE_ID, V1, label = V1), vjust = -0.2) +
    facet_wrap(~chain_type) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylab("Number of unique chain sequences detected") +
    ggtitle("number of Ig chains sequenced")
  ggsave(paste0(report_folder,"/","step5_finalClones_NChains_sequenced_bar.pdf"), height = 6, width = 12)

  ggplot(IGSeq_resultset$step5_finalclones[, length(unique(cloneId)), .(SAMPLE_ID, chain_type)]) +
    geom_point(aes(SAMPLE_ID, V1, color = SAMPLE_ID, shape = chain_type), size = 3) +
    geom_text_repel(aes(SAMPLE_ID, V1, label = paste(chain_type,"\n",V1)), nudge_x = 0.2) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Number of unique chain sequences detected") +
    ggtitle("IGSeq number of chains sequenced")
  ggsave(paste0(report_folder,"/","step5_finalClones_NChains_sequenced.pdf"), height = 7, width = 7)

  # Plot isotype distribution
  ggplot(IGSeq_resultset$step5_finalclones[, length(unique(cloneId)), .(SAMPLE_ID, chain_isotype)]) +
    geom_bar(aes(SAMPLE_ID, V1, fill = chain_isotype), stat = "identity", position = "dodge") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Number of unique chain sequences detected") +
    ggtitle("IGSeq IGH chain isotypes observation")
  ggsave(paste0(report_folder,"/","step5_finalClones_IGH_isotype_obs.pdf"), height = 7, width = 7)

}
