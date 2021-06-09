#' @title importDIAresults
#' @description import MS results from Spectronaut
#' @param report_location path to the Spectronaut report (Factory default)
#' @return DIAresult object (list) with
#' * data_long - data.table of the report
#' * overview - analysis overview from SN
#' * settings - settings file
#' @import data.table
#' @export

importDIAresults = function(report_location){

  path = dirname(report_location)
  print(path)

  data_long = fread(report_location)
  overview = fread(paste0(path, "/AnalyisOverview.txt"))
  message(list.files(path)[grep("ExperimentSetup", list.files(path))])

  if (length(list.files(path)[grep("ExperimentSetup", list.files(path))])>0){
    settings = fread(paste0(path, "/", list.files(path)[grep("ExperimentSetup", list.files(path))]))
  } else {
    settings = "N/A"
  }


  DIAresult = list("data_long" = data_long,
                   "overview" = overview,
                   "settings" = settings)

  return(DIAresult)
}
