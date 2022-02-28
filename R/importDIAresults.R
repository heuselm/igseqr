#' @title importDIAresults
#' @description import DIA-MS quantitative results from Spectronaut
#' @param report_location path to the Spectronaut report (Factory default)
#' @param study_design alternative study design to use, data.table of filename, condition, replicate.
#' If not provided, an attempt will be made to extract the study design from the imported DIA results.
#' @return DIA_resultset object (list) with
#' data_long - data.table of the report
#' data_wide - data.table of Protein.Group+Precursor.Id over the Ms runs/R.FileName
#' study_design - study design table of filename, condition, replicate
#' @import data.table
#' @export

importDIAresults = function(report_location, study_design = NULL){

  # get the path
  path = dirname(report_location)
  print(path)

  # read full long format report
  data_long = fread(report_location)

  # extract study design if available
  if(is.null(study_design)){
    study_design = unique(data_long[, .(R.FileName, R.Condition, R.Replicate)])
    names(study_design) = c("filename", "condition", "replicate")
  } else if (all(c("filename", "condition", "replicate") %in% names(study_design))){
    study_design = study_design
  } else{
      warning("study_design does not contain necessary columns filename, condition, replicate")
    }

  # reformat to wide:
  data_long_s = data_long[, .(PG.ProteinGroups, FG.Charge, FG.LabeledSequence, FG.Quantity, R.FileName)]
  data_long_s[, Precursor.Id:=paste0(gsub("_", "", FG.LabeledSequence), FG.Charge), row.names(data_long_s)]
  data_long_s = merge(data_long_s, study_design, by.x = "R.FileName", by.y = "filename", all.y = F)
  setnames(data_long_s, "PG.ProteinGroups", "Protein.Group")
  data_wide = dcast(data_long_s, Protein.Group+Precursor.Id~R.FileName,
                    value.var = "FG.Quantity")

  # Assemble result
  DIAresult = list("data_long" = data_long,
                   "data_wide" = data_wide,
                   "study_design" = study_design)

  return(DIAresult)
}
