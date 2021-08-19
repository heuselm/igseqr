summarizeMSdetection = function(IGSeq_resultset_MSdetection){

    # Assemble & write combined table from MSannotation results
    table2 = IGSeq_resultset_MSdetection$MS_abundance_allprecursors_mean
    colnames(table2) = paste("MS_abundance_allprecursors_mean",colnames(table2))

    table3 = IGSeq_resultset_MSdetection$MS_abundance_uniqueprecursors_mean
    colnames(table3) = paste("MS_abundance_uniqueprecursors_mean", colnames(table3))

    table4 = IGSeq_resultset_MSdetection$MS_abundance_uniqueprecursors_mean_cond
    colnames(table4) = paste("MS_abundance_uniqueprecursors_mean_cond", colnames(table4))

    table5 = IGSeq_resultset_MSdetection$MS_abundance_uniqueprecursors_sd_cond
    colnames(table5) = paste("MS_abundance_uniqueprecursors_sd_cond", colnames(table5))

    table6 = IGSeq_resultset_MSdetection$MS_detection_uniqueprecursors_nrep_cond
    colnames(table6) = paste("MS_detection_uniqueprecursors_nrep_cond", colnames(table6))

    clone_detection_and_quant = cbind(IGSeq_resultset_MSdetection$clones_detected,
                                      table2,table3,table4,table5,table6)
    fwrite(clone_detection_and_quant,
           file = paste0(deparse(substitute(IGSeq_resultset_MSdetection)), ".csv"))

}


