#' Export the SummarizedExperiment object, the meta data, and the normalized data.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param out_dir Path of output directory
#' @param ain Vector of strings which assay should be downloaded (default NULL).
#'            If NULL then all assays of the se object are saved.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' \dontrun{export_data(tuberculosis_TMT_se, out_dir = "data/",
#'             ain = c("IRS_on_RobNorm", "IRS_on_Median"))}
#'
export_data <- function(se, out_dir, ain = NULL) {
  ain <- check_input_assays(se, ain)

  # Write column data
  message(".. saving meta_data.csv")
  utils::write.csv(data.frame(SummarizedExperiment::colData(se), check.names=FALSE), file.path(out_dir, "meta_data.csv"), row.names = FALSE)

  # Extract Assay Data: The assay data is the core data matrix (or matrices) in a SummarizedExperiment object.
  message(".. saving assays as .csv")
  row_data <- data.frame(SummarizedExperiment::rowData(se), check.names=FALSE)
  assays <- SummarizedExperiment::assays(se)
  for (assay in ain) {
    cols <- names(SummarizedExperiment::assays(se)[[assay]])
    assay_data <- data.frame(SummarizedExperiment::assays(se)[[assay]])
    colnames(assay_data) <- cols
    combined_data <- cbind(row_data, assay_data)
    combined_data$IDs <- NULL
    utils::write.csv(combined_data, file.path(out_dir, paste0(assay, "_normalized_data.csv")), row.names = FALSE)
  }
  message(".. saving SummarizedExperiment as RDS file")
  # Save the SummarizedExperiment object as an RDS file
  saveRDS(se, file.path(out_dir, "se.rds"))
}
