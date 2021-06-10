#' @title importDIAresults
#' @description import MS results from Spectronaut
#' @param report_location path to the Spectronaut report (Factory default)
#' @return DIAresult object (list) with
#' data_long - data.table of the report
#' data_wide - data.table of Protein.Group+Precursor.Id over the Ms runs/R.FileName
#' overview - analysis overview from SN
#' settings - settings file
#' @import data.table
#' @export

importDIAresults = function(report_location){

  # get the path
  path = dirname(report_location)
  print(path)

  # read full long format report
  data_long = fread(report_location)

  # read additional metainfo of SN analysis, if available
  overview = fread(paste0(path, "/AnalyisOverview.txt"))
  if (length(list.files(path)[grep("ExperimentSetup", list.files(path))])>0){
    settings = fread(paste0(path, "/", list.files(path)[grep("ExperimentSetup", list.files(path))]))
  } else {
    settings = "N/A"
  }

  # extract study design if available
  study_design = unique(data_long[, .(R.FileName, R.Condition, R.Replicate)])
  names(study_design) = c("filename", "condition", "replicate")

  # reformat to wide:
  data_long_s = data_long[, .(PG.ProteinGroups, FG.Charge, FG.LabeledSequence, FG.Quantity, R.FileName)]
  data_long_s[, Precursor.Id:=paste0(gsub("_", "", FG.LabeledSequence), FG.Charge), row.names(data_long_s)]
  setnames(data_long_s, "PG.ProteinGroups", "Protein.Group")
  data_wide = dcast(data_long_s, Protein.Group+Precursor.Id~R.FileName,
                    value.var = "FG.Quantity")

  # Assemble result
  DIAresult = list("data_long" = data_long,
                   "data_wide" = data_wide,
                   "study_design" = study_design,
                   "overview" = overview,
                   "settings" = settings)

  return(DIAresult)
}
