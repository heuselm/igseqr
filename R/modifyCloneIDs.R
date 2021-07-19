#' @title modifyCloneIDs
#'
#' @description modify Clone IDs (or protein identifiers) within vector of protein groups.
#' Useful when non-trivial modifications are required, such as replacement with different information
#'  from cloneID table.
#'
#' @param cloneIDs_to_modify Character vector of protein groups including cloneIDs (with internal | char.)
#' @param modification_table modification mapping data table with columns cloneID_old & cloneID_new
#'
#' @return mapping table with input protein groups in column cloneIDs_to_simplify and modified
#' @import data.table
#'
#' @export

modifyCloneIDs = function(cloneIDs_to_modify, modification_table){

  if (!all(c("cloneID_old", "cloneID_new") %in% names(modification_table))){
    stop("Modification table must contain columns cloneID_old and cloneID_new")
  }

  cloneIDs_modified = copy(cloneIDs_to_modify)

  # subset modification table to those truly present in cloneIDs_to_modify
  modification_table[, present:=length(grep(gsub("\\|","\\\\|", cloneID_old),
                                            cloneIDs_to_modify)), cloneID_old]
  # Modify each protein group
  for (i in 1:nrow(modification_table[present>0])){
    cloneIDs_modified = gsub(gsub("\\|","\\\\|", modification_table[present>0]$cloneID_old[i]),
                               modification_table[present>0]$cloneID_new[i], cloneIDs_modified)
  }

  return(data.table(cloneIDs_to_modify, cloneIDs_modified))
}

