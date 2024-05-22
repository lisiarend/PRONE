
#' Generate SummarizedExperiment object containing the samples of all assays in a single object (for overall PCA)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input (default NULL)
#'            If NULL then all normalization of the SummarizedExperiment object are plotted next to each other (except raw).
#'
#' @return SummarizedExperiment object with single assay containing the normalized samples of all assays
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se_all <- generate_complete_SE(tuberculosis_TMT_se,
#'                                                 ain = NULL)
generate_complete_SE <- function(se, ain = NULL){
  # create one big data table
  big_dt <- NULL
  big_cd <- NULL
  # check if assays all in se
  ain <- check_input_assays(se, ain)
  for(a in ain){
    if(!a %in% c("raw")){
      dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[a]])
      colnames(dt) <- paste0(a, "_", colnames(dt))
      cd <- data.table::as.data.table(SummarizedExperiment::colData(se))
      cd$Column <- paste0(a,"_",cd$Column)
      cd$Normalization <- a
      if(is.null(big_dt)){
        big_dt <- dt
        big_cd <- cd
      } else {
        big_dt <- cbind(big_dt, dt)
        big_cd <- rbind(big_cd, cd)
      }
    }
  }
  se_big <- SummarizedExperiment::SummarizedExperiment(assays = list(all = big_dt), colData = big_cd, rowData = data.table::as.data.table(SummarizedExperiment::rowData(se)))
  S4Vectors::metadata(se_big) <- S4Vectors::metadata(se)
  return(se_big)
}
