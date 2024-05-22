

#' Remove normalization assays from a SummarizedExperiment object
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param assays_to_remove Character vector of assay names to remove from the SummarizedExperiment object
#'
#' @return SummarizedExperiment object with the normalization assays removed
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- remove_assays_from_SE(tuberculosis_TMT_se,
#'                                  assays_to_remove = c("IRS_on_RobNorm"))
#'
remove_assays_from_SE <- function(se, assays_to_remove){
  for(ain in assays_to_remove){
    # check if ain in se
    if(ain %in% names(SummarizedExperiment::assays(se))){
      # remove ain from se
      SummarizedExperiment::assays(se)[[ain]] <- NULL
    } else {
      warning("Assay ", ain, " not found in se!")
    }
  }
  return(se)
}

#' Subset SummarizedExperiment object by normalization assays
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Character vector of assay names to keep in the SummarizedExperiment object
#'
#' @return SummarizedExperiment object with only the selected normalization assays
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- subset_SE_by_norm(tuberculosis_TMT_se,
#'                           ain = c("raw", "log2", "IRS_on_RobNorm"))
#'
subset_SE_by_norm <- function(se, ain){
  not_available_ains <- ain[!ain %in% names(SummarizedExperiment::assays(se))]
  if(length(not_available_ains) > 0){
    warning("Assay(s) ", paste(not_available_ains, collapse = ", "), " not found in se!")
  }
  assays_to_remove <- names(SummarizedExperiment::assays(se))[!names(SummarizedExperiment::assays(se)) %in% ain]
  se <- remove_assays_from_SE(se, assays_to_remove)
  return(se)
}
