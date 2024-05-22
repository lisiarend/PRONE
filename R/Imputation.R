
#' Method to impute SummarizedExperiment. 
#' This method performs a mixed imputation on the proteins. It uses a k-nearest neighbor imputation for proteins with missing values at random (MAR) and imputes missing values by random draws from a left-shifted Gaussian distribution for proteins with missing values not at random (MNAR).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics dataset
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition name of column of colData(se) representing the conditions of the data
#'
#' @return SummarizedExperiment with imputed intensities
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- remove_samples_manually(tuberculosis_TMT_se,
#'                           column = "Label", values = c("1.HC_Pool1"))
#' tuberculosis_TMT_se <- impute_se(tuberculosis_TMT_se, ain = NULL,
#'                           condition = NULL)
#'
#' @importFrom MSnbase MSnSet
impute_se <- function(se, ain = NULL, condition = NULL){
  se_imputed <- se
  # check input
  condition <- get_condition_value(se, condition)
  assays <- check_input_assays(se, ain)

  for(ain in assays){
    dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
    dt$ID <- rownames(dt)
    dt <- data.table::melt(dt, id.vars = "ID", variable.name = "Column", value.name = "Intensity")
    coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
    dt <- merge(dt, coldata, by = "Column")

    proteins_MNAR <- dt %>%
      dplyr::group_by(ID, get(condition)) %>%
      dplyr::summarize(NAs = all(is.na(Intensity))) %>%
      dplyr::filter(NAs) %>%
      dplyr::pull(ID) %>%
      unique()

    # Get a logical vector
    original_dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
    MNAR <- rownames(original_dt) %in% proteins_MNAR
    MAR <- !rownames(original_dt) %in% proteins_MNAR

    # Impute MAR
    # Perform knn imputation using the MSnBase impute method
    MAR_dt <- original_dt[MAR,,drop=FALSE] # use Msnbase impute method with fun="knn
    # Use MsnSet::impute function
    MSnSet_dt <- methods::new("MSnSet",exprs=as.matrix(MAR_dt))
    MsnSet_imputed <-  MSnbase::impute(MSnSet_dt,method="knn")
    MAR_dt <- MSnbase::exprs(MsnSet_imputed)

    # Impute MNAR
    # Rebuild DEP method (impute_manual)
    MNAR_dt <- original_dt[MNAR,,drop=FALSE]

    MNAR_prot_id <- row.names(MNAR_dt)
    MNAR_dt <- as.data.frame(MNAR_dt)
    MNAR_dt <- cbind(ID = MNAR_prot_id , MNAR_dt)

    scale <- 0.3
    shift <- 2.0

    # Get descriptive parameters of the current sample distributions
    stat <- MNAR_dt %>%
      tidyr::gather(Column, Intensity, -ID) %>%
      dplyr::filter(!is.na(Intensity)) %>%
      dplyr::group_by(Column) %>%
      dplyr::summarise(mean = mean(Intensity), median = stats::median(Intensity), sd = stats::sd(Intensity), n = dplyr::n(), infin = nrow(MNAR_dt) - n)

    names(stat) <- c("Column", "mean", "median", "sd", "n", "infin")
    # Impute missing values by random draws from a distribution
    # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
    for (a in seq_len(nrow(stat))){
      MNAR_dt[is.na(MNAR_dt[, stat[, "Column"][[1]][a]]), stat[, "Column"][[1]][a]] <-
        stats::rnorm(stat$infin[a],
              mean = stat$median[a] - shift * stat$sd[a],
              sd = stat$sd[a] * scale)
    }
    # Remove ID
    MNAR_dt$ID <- NULL

    # Combine MNAR imputed counts and MAR imputed counts
    imputed_dt <- rbind(MAR_dt, MNAR_dt)
    SummarizedExperiment::assay(se_imputed, paste0(ain, "_imputed"), FALSE) <- data.table::as.data.table(imputed_dt)
  }
  return(se_imputed)
}
