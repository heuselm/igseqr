#' @title assembleFastaDb
#'
#' @description assemble fasta database of variable domains from MiXCR IgSeq results
#'
#' @param IGSeq_resultset IgSeq resultset imported by importIgSeq function.
#' @param dataset_tag Dataset tag introduced to cloneIdGlobal and output files. Default: Input IGSeq R object name
#' @param write_fasta Whether to write .fasta files to working directory. Default: TRUE
#' @param write_fasta_per_sample Whether to write separate .fasta files per sample. Default: FALSE
#' @param write_clone_tables Whether to write clone .tsv tables next to .fasta files. Default: TRUE
#' @param subset_cloneIdGlobal Whether to subset clones to be included in fasta file. Character vector of global clone IDs.
#' @param gsub_sample_id character vector of length 2, Whether to replace certain string in sample IDs (element 1) with (element 2). Default: Null (no editing of SAMPLE_ID)
# @param prepend_chaintype whether chain_type should be prepended to(at the beginning of) clonIdGlobal
#'
#' Default: Null, All clones included.
#'
#' @return IGSeq_resultset with added element step6_clonestofasta
#'
#' @import data.table ggplot2 ggrepel seqinr
#'
#' @export

assembleFastaDb = function(IGSeq_resultset, dataset_tag = NULL, write_fasta = TRUE, write_fasta_per_sample = FALSE, write_clone_tables = TRUE, write_fasta_subset_cloneIdGlobal = NULL, gsub_sample_id = NULL){

  # Get Dataset tag if not provided
  if (is.null(dataset_tag)){
    dataset_tag = deparse(substitute(IGSeq_resultset))
  }

  # Clean missing values /rows where one sequence is missing
  # & Assemble longest-possible Var-domain sequences
  ##########################################################

  if (!is.null(gsub_sample_id)){
    IGSeq_resultset$step5_finalclones[, SAMPLE_ID:=gsub(gsub_sample_id[1],gsub_sample_id[2], SAMPLE_ID)]
  }

  IGSeq_resultset$step5_finalclones[, cloneIdGlobal:=paste0(dataset_tag, SAMPLE_ID, "_", cloneId), .(cloneId, SAMPLE_ID)]
  IGSeq_resultset$step6_clonestofasta = copy(IGSeq_resultset$step5_finalclones)

  # Assemble metainfo for Protein Description
  IGSeq_resultset$step6_clonestofasta[, aaSeq_lengths:=paste(nchar(aaSeqFR1), nchar(aaSeqCDR1), nchar(aaSeqFR2), nchar(aaSeqCDR2), nchar(aaSeqFR3), nchar(aaSeqCDR3), nchar(aaSeqFR4), collapse = ";"), cloneIdGlobal]
  IGSeq_resultset$step6_clonestofasta[, aaSeq_AllComplete:=!grepl(" 0 ", aaSeq_lengths), cloneIdGlobal]
  IGSeq_resultset$step6_clonestofasta[, aaSeq_fullVarDomain:=gsub("_", "", paste0(aaSeqFR1,aaSeqCDR1,aaSeqFR2,aaSeqCDR2,aaSeqFR3,aaSeqCDR3,aaSeqFR4)), cloneIdGlobal]

  # kick out incomplete/ambiguous sequences containing *
  IGSeq_resultset$step6_clonestofasta[, aaSeq_AllComplete :=!grepl("\\*", aaSeq_fullVarDomain), row.names(IGSeq_resultset$step6_clonestofasta)]
  IGSeq_resultset$step6_clonestofasta[, aaSeq_AllComplete :=!grepl("_", aaSeq_fullVarDomain), row.names(IGSeq_resultset$step6_clonestofasta)]

  # Make fasta headers
  IGSeq_resultset$step6_clonestofasta[, fasta_header:=gsub(",", ";", paste0("igseq|", cloneIdGlobal, "|", chain_type, "_", cloneIdGlobal, " ", " primaryCHitWithScore: ",
                                                          substr(allCHitsWithScore, 1,15), " vdj_genes:", allVHitsWithScore, allDHitsWithScore, allJHitsWithScore,
                                                          " domain_lengths(fr1-cdr1-fr2-cdr2-fr3-cdr3-fr4):", aaSeq_lengths,
                                                          " countInSeq:", cloneCount)), cloneIdGlobal]
  IGSeq_resultset$step6_clonestofasta[, subset:=TRUE]
  if(!(is.null(write_fasta_subset_cloneIdGlobal))){
    message("Clones to fasta are being subsetted")
    IGSeq_resultset$step6_clonestofasta[, subset:= cloneIdGlobal %in% write_fasta_subset_cloneIdGlobal, cloneIdGlobal]
    message("Clones per sample total:", IGSeq_resultset$step6_clonestofasta[, length(unique(cloneId)), .(SAMPLE,chain_type)])
    message("Clones per sample subset:", IGSeq_resultset$step6_clonestofasta[subset == TRUE, length(unique(cloneId)), .(SAMPLE,chain_type)])
  }

  if (write_fasta){
    # Write out final clones in tsv and fasta format
    if(write_clone_tables){
      write.table(IGSeq_resultset$step6_clonestofasta[subset == TRUE], paste0(dataset_tag, "final_clones_annotated.tsv"), sep = "\t", quote = F, row.names = F)
    }
    write.fasta(as.list(IGSeq_resultset$step6_clonestofasta[subset == TRUE]$aaSeq_fullVarDomain), as.string = TRUE,
                names = as.list(IGSeq_resultset$step6_clonestofasta[subset == TRUE]$fasta_header), file.out = paste0(dataset_tag, "final_clones_annotated.fasta"))
    if (write_fasta_per_sample){
      samples = unique(IGSeq_resultset$step6_clonestofasta[subset == TRUE]$SAMPLE)
      # Write out individual, filtered tables and fasta file(s)
      for (i in seq_along(samples)){
        if(write_clone_tables){
        write.table(IGSeq_resultset$step6_clonestofasta[subset == TRUE][SAMPLE == samples[i]],
                    file = paste0(dataset_tag, samples[i], "final_clones_annotated.tsv"),
                    sep = "\t", quote = F, row.names = F)
        }
        write.fasta(as.list(IGSeq_resultset$step6_clonestofasta[subset == TRUE][SAMPLE == samples[i]]$aaSeq_fullVarDomain),
                    as.string = TRUE,
                    names = as.list(IGSeq_resultset$step6_clonestofasta[subset == TRUE][SAMPLE == samples[i]]$fasta_header),
                    file.out = paste0(dataset_tag, samples[i], "final_clones_annotated.fasta"))
      }
      }
    }

    return(IGSeq_resultset)
    }
