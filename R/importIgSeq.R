#' @title summarizeIgSeqResults
#' @description summarize MiXCR IgSeq results as structured by https://github.com/heuselm/igseq_workflow.git
#' @param result_folder Location of the IgSeq-workflow result folder.
#' @param write_tables Switch whether the collected summary information should be output as .csv tables
#' @value IGSeq result object (list) with tables
#' * step1_checkout_log
#' * step2_histogram_log
#' * step2_histogram_osq
#' * step3_assemble_log
#' * step4_align_log
#' * step5_cloneAssembly_log"
#' * step5_finalclones
#' @import data.table ggplot2 ggrepel
#' @export

importIgSeq <- function(result_folder = ".", write_tables = TRUE, summarize_in_pdfs = TRUE){

  ### collect & process / formalize results

  ## Step 1 - checkout
  step1_checkout = fread(paste0(result_folder,  "01_checkout/checkout.log.txt"))
  if ("SAMPLE" %in% names(step1_checkout)){
    samples = sort(step1_checkout$SAMPLE)
  } else {
    samples = sort(step1_checkout$V3[grep("S", step1_checkout$V3)])
  }

  ## Step 2 - histogram stats
  step2_histogram = fread(paste0(result_folder,  "02_histogram/estimates.txt"))
  step2_histogram_osq = fread(paste0(result_folder,  "02_histogram/overseq.txt"), header = T)

  ## Step 3 - assembly of UMI groups
  step3_assemble = fread(paste0(result_folder,  "03_assemble/assemble.log.txt"))

  ## Step 4 - Read Alignment to Germline Sequences
  step4_align = data.table()
  step4_align_files = list.files(paste0(result_folder,  "04_aligned/"))[grep("alignmentReport", list.files(paste0(result_folder,  "04_aligned/")))]
  for (i in seq_along(step4_align_files)){
    dt.i = fread(paste0(paste0(result_folder,  "04_aligned/"), step4_align_files[i]), sep = ":",
                 fill = T)[, SAMPLE:=gsub("alignmentReport_|.txt","",step4_align_files[i])]
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
  step5_cloneAssemblyReports_files = list.files(paste0(result_folder,  "05_finalclones/"))[grep("assemblyReport", list.files(paste0(result_folder,  "05_finalclones/")))]
  for (i in seq_along(step5_cloneAssemblyReports_files)){
    dt.i = fread(paste0(paste0(result_folder,  "05_finalclones/"), step5_cloneAssemblyReports_files[i]),
                 sep = "X", header = F)[, SAMPLE:=gsub("cloneassemblyReport|.txt","",step5_cloneAssemblyReports_files[i])]
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
  step5_finalclones_files = list.files(paste0(result_folder,  "05_finalclones/"))[grep("clones_all", list.files(paste0(result_folder,  "05_finalclones/")))]
  for (i in seq_along(step5_finalclones_files)){
    dt.i = fread(paste0(paste0(result_folder,  "05_finalclones/"),
                        step5_finalclones_files[i]))[, SAMPLE:=gsub("_clones_all.txt", "", step5_finalclones_files[i])]
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

  ## Write tables
  if (write_tables){
    fwrite(step1_checkout, paste0("summarizeIgSeq_", "step1_checkout.csv"))
    fwrite(step2_histogram, paste0("summarizeIgSeq_", "step2_histogram.csv"))
    fwrite(step3_assemble, paste0("summarizeIgSeq_", "step3_assemble.csv"))
    fwrite(step4_align, paste0("summarizeIgSeq_", "step4_align.csv"))
    fwrite(step5_cloneAssemblyReports_num, paste0("summarizeIgSeq_", "step5_cloneAssemblyReports_num.csv"))
    fwrite(step5_finalclones, paste0("summarizeIgSeq_", "step5_finalclones.csv"))
  }

  res = list("step1_checkout_log" = step1_checkout,
             "step2_histogram_log" = step2_histogram,
             "step2_histogram_osq" = step2_histogram_osq,
             "step3_assemble_log" = step3_assemble,
             "step4_align_log" = step4_align_num,
             "step5_cloneAssembly_log" = step5_cloneAssemblyReports_num,
             "step5_finalclones" = step5_finalclones)

  if(summarize_in_pdfs){
    summarizeIgSeq(res)
  }

  return(res)
}



