# helper functions

matrix_please = function(datatable, rn_colidx = 1, transform_colnames = NULL){
  m = as.matrix(datatable[, 2:ncol(datatable)])
  row.names(m) = datatable[[names(datatable[, rn_colidx])]]

  if(!is.null(transform_colnames)){
    colnames_new = sapply(colnames(m), function(x){transform_colnames})
    colnames(m) = colnames_new
  }
  return(m)
}
