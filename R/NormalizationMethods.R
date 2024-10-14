## ----- Normalization Methods ----- ##

#' Total Intensity Normalization
#'
#' Intensities of each variable in a sample are divided with the sum of intensities
#' of all variables in the sample and multiplied with the median or mean of sum of intensities
#' of all variables in all samples. Raw data should be taken as input (on_raw = TRUE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param type String whether to use median or mean to calculate the scaling factor
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#' @importFrom SummarizedExperiment assay
#'
#' @return SummarizedExperiment containing the total intensity normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- globalIntNorm(tuberculosis_TMT_se, ain = "raw",
#'                                            aout = "GlobalMedian",
#'                                            type = "median",
#'                                            on_raw = TRUE)
#'
globalIntNorm <- function(se, ain = "raw", aout="GlobalMedian", type = "median", on_raw = TRUE){
  dt <- SummarizedExperiment::assays(se)[[ain]]
  dt <- as.data.frame(dt)
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  colSums <- colSums(dt, na.rm = TRUE)
  if(type == "median"){
    colSumsM <- stats::median(colSums)
  } else if (type == "mean"){
    colSumsM <- mean(colSums)
  }
  norm_dt <- matrix(nrow = nrow(dt), ncol = ncol(dt),
                    byrow = TRUE)
  normFunc <- function(colIndex) {
    (dt[rowIndex, colIndex]/colSums[colIndex]) *
      colSumsM
  }
  for (rowIndex in seq_len(nrow(dt))) {
    norm_dt[rowIndex, ] <- vapply(seq_len(ncol(dt)),
                                  normFunc, 0)
  }
  # transform data to log2-scale if necessary
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  colnames(norm_dt) <- colnames(dt)
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Total Intensity Normalization Using the Mean for the Calculation of Scaling Factors
#'
#' Intensities of each variable in a sample are divided with the sum of intensities
#' of all variables in the sample and multiplied with the mean of sum of intensities
#' of all variables in all samples. Raw data should be taken as input (on_raw = TRUE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the total intensity normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- globalMeanNorm(tuberculosis_TMT_se, ain = "raw",
#'                               aout = "GlobalMean", on_raw = TRUE)
#'
globalMeanNorm <- function(se, ain = "raw", aout = "GlobalMean", on_raw = TRUE){
  se <- globalIntNorm(se, ain, aout, type = "mean", on_raw = on_raw)
  return(se)
}

#' Total Intensity Normalization Using the Median for the Calculation of Scaling Factors
#'
#' Intensities of each variable in a sample are divided with the sum of intensities
#' of all variables in the sample and multiplied with the median of sum of intensities
#' of all variables in all samples. Raw data should be taken as input (on_raw = TRUE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the total intensity normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- globalMedianNorm(tuberculosis_TMT_se, ain = "raw",
#'                               aout = "GlobalMedian", on_raw = TRUE)
#'
globalMedianNorm <- function(se, ain = "raw", aout = "GlobalMedian", on_raw = TRUE){
  se <- globalIntNorm(se, ain, aout, type = "median")
  return(se)
}


#' Median Normalization
#'
#' The intensity of each protein group in a given sample is divided by the median of the
#' intensities of all protein groups in that sample and then multiplied by the mean of
#' median of sum of intensities of all protein groups in all samples.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the median normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- medianNorm(tuberculosis_TMT_se, ain = "raw",
#'                           aout = "Median", on_raw = TRUE)
#'
medianNorm <- function(se, ain = "raw", aout = "Median", on_raw = TRUE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  # find median of each sample
  sample_med <- apply(dt, 2, stats::median, na.rm=TRUE) # columns
  # find mean of medians
  mean_med <- mean(sample_med, na.rm=TRUE)
  # divide data by median
  norm_dt <- t(t(dt)/sample_med)
  # multiply data by mean of medians
  norm_dt <- norm_dt * mean_med
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  # transform data to log2-scale if necessary
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Mean Normalization
#'
#' The intensity of each protein group in a given sample is divided by the mean of the
#' intensities of all protein groups in that sample and then multiplied by the mean of
#' mean of sum of intensities of all protein groups in all samples.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean whether normalized should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the mean normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- meanNorm(tuberculosis_TMT_se, ain = "raw",
#'                           aout = "Mean", on_raw = TRUE)
#'
meanNorm <- function(se, ain = "raw", aout="Mean", on_raw = TRUE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  # find means of each sample
  sample_mean <- apply(dt, 2, mean, na.rm=TRUE)
  # find mean of means of each sample
  mean_mean <- mean(sample_mean, na.rm=TRUE)
  # divide data by mean
  norm_dt <- t(t(dt)/sample_mean)
  # multiply data by mean of means
  norm_dt <- norm_dt * mean_mean
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  # transform data to log2-scale if necessary
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Internal Reference Scaling Normalization
#'
#' IRS makes different measurements of the same thing all exactly the same and puts
#' all of the intensities on the same scale. Raw data should be taken as input (on_raw = TRUE)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the IRS normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- irsNorm(tuberculosis_TMT_se, ain = "raw",
#'                                 aout = "IRS", on_raw = TRUE)
irsNorm <- function(se, ain="raw", aout="IRS", on_raw = TRUE){
  # extract necessary info of SummarizedExperiment
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  batch <- S4Vectors::metadata(se)$batch
  refs <- S4Vectors::metadata(se)$refs
  md <- data.table::as.data.table(SummarizedExperiment::colData(se))

  # get md of reference samples
  refs_md <- md[md$Column %in% refs,]
  # separate data by batch
  dt_list <- lapply(unique(md[[batch]]), function(b){
    md_chunk <- md[md[[batch]] == b,]
    dt_chunk <- subset(dt, select = md_chunk$Column)
    return(dt_chunk)
  })
  names(dt_list) <- unique(md[[batch]])
  # check if ref_samples have been specified
  if (is.null(refs)){
    # create mock channel with rowsums from each frame
    irs <- lapply(dt_list, function(dt_chunk){
      return(rowSums(dt_chunk, na.rm=TRUE))
    })
    irs <- do.call(cbind, irs)
  } else {
    # take reference sample intensities
    irs <- subset(dt, select = refs_md$Column)
    colnames(irs) <- as.character(refs_md[refs_md$Column %in% refs,][[batch]])
  }
  # get the geometric average intensity for each protein
  irs <- tibble::as_tibble(irs)
  irs$average <- apply(irs, 1, function(x) exp(mean(log(x), na.rm=TRUE)))
  # normalize data
  dt_irs_list <- lapply(names(dt_list), function(b){
    # compute scaling factor vectors
    fac <- irs$average / irs[,b]
    # normalize
    dt_irs_chunk <- dt_list[[b]] * fac[,1]
    return(dt_irs_chunk)
  })
  # reconstruct data after irs normalization
  dt_irs <- do.call(cbind, dt_irs_list)
  dt_irs <- data.table::as.data.table(dt_irs)
  dt_irs <- subset(dt_irs, select = colnames(dt))
  # transform data to log2-scale if necessary
  if(on_raw){
    dt_irs <- log2(dt_irs)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- dt_irs
  return(se)
}


#' Quantile Normalization of preprocessCore package.
#'
#' Forces distributions of the samples to be the same on the basis of the quantiles of the samples by replacing
#' each protein of a sample with the mean of the corresponding quantile. Log2-scaled data should be taken as input (on_raw = FALSE)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the quantile normalized data as assay (on log2 scale)
#' @export
#' @seealso \code{\link[preprocessCore]{normalize.quantiles}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- quantileNorm(tuberculosis_TMT_se, ain = "log2",
#'                               aout = "Quantile", on_raw = FALSE)
#'
quantileNorm <- function(se, ain="log2", aout="Quantile", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  if(on_raw & ain != "raw"){
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  if(ain == "log2" & on_raw == TRUE){
    warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    dt <- 2^dt
  }
  norm_dt <- preprocessCore::normalize.quantiles(as.matrix(dt), copy=TRUE)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  # transform data to log2-scale if necessary
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Variance Stabilization Normalization of limma package.
#'
#' Raw data should be taken as input (on_raw = TRUE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#' @param VSN_quantile Numeric of length 1. The quantile that is used for the resistant least trimmed sum of squares regression (see vsn2 lts.quantile)
#'
#' @return SummarizedExperiment containing the vsn normalized data as assay (on log2-scale)
#' @export
#' @seealso \code{\link[limma]{normalizeVSN}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- vsnNorm(tuberculosis_TMT_se, ain = "raw",
#'                             aout = "VSN", on_raw = TRUE, VSN_quantile = 0.9)
#'
vsnNorm <- function(se, ain="raw", aout="VSN", on_raw = TRUE, VSN_quantile = 0.9){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  norm_dt <- suppressMessages(limma::normalizeVSN(dt, lts.quantile = VSN_quantile))
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  # no log transformation because limma already returns log2 data
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}


#' Weighted Trimmed Mean of M Values (TMM) Normalization of edgeR package.
#'
#' Raw data should be taken as input (on_raw = TRUE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the TMM normalized data as assay (on log2 scale)
#' @export
#' @seealso \code{\link[edgeR]{calcNormFactors}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- tmmNorm(tuberculosis_TMT_se, ain = "raw",
#'                                   aout = "TMM", on_raw = TRUE)
#'
tmmNorm <- function(se, ain="raw", aout="TMM", on_raw = TRUE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  tmm <- edgeR::calcNormFactors(stats::na.omit(dt))
  dt_norm <- sweep(dt, 2, tmm, FUN="/")
  dt_norm <- data.table::as.data.table(dt_norm)
  if(on_raw){
    dt_norm <- log2(dt_norm)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- dt_norm
  return(se)
}

#' Robust Linear Regression Normalization of NormalyzerDE.
#'
#' Uses median values over all samples as reference
#' sample to which all the other samples in the data are normalized to. Log2 data should be taken as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the rlr normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- rlrNorm(tuberculosis_TMT_se, ain = "log2",
#'                                             aout = "Rlr", on_raw = FALSE)
#'
rlrNorm <- function(se, ain="log2", aout="Rlr", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt <- as.data.frame(dt)
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  sample_median <- matrixStats::rowMedians(as.matrix(dt), na.rm = TRUE, useNames = TRUE)
  for (i in seq_len(ncol(dt))){ # iterate over samples
    sampleA <- dt[,i] # sample to normalize
    lrFit <- MASS::rlm(as.matrix(sampleA) ~
                         sample_median, na.action = stats::na.exclude)
    coeffs <- lrFit$coefficients
    coefIntercept <- coeffs[1]
    coefSlope <- coeffs[2]
    globalFittedRLRCol <- (sampleA - coefIntercept)/coefSlope
    dt[,i] <- globalFittedRLRCol
  }
  colnames(dt) <- colnames(dt)
  rownames(dt) <- rownames(dt)
  if(on_raw){
    dt <- log2(dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(dt)
  return(se)
}

#' Linear Regression Normalization on MA Transformed Data
#'
#' Similar to Rlr, but data are MA transformed before normalization,
#' (A = median sample, M = difference of that sample to A). Log2 data should be taken as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the RlrMA normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- rlrMANorm(tuberculosis_TMT_se, ain = "log2",
#'                                   aout = "RlrMA", on_raw = FALSE)
#'
rlrMANorm <- function(se, ain="log2",aout="RlrMA", on_raw = FALSE){
  # extract intensities
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  dt_matrix <- as.matrix(dt)
  # MA transformation
  ref_a <- apply(dt_matrix, 1, stats::median, na.rm = TRUE)
  for (i in seq_len(ncol(dt))){ # iterate over samples
    comp_sample <- dt_matrix[, i] # sample to normalize
    m <- comp_sample - ref_a  # M = sample to normalize
    # rlm fit
    fit <- MASS::rlm(m ~ ref_a, na.action = stats::na.exclude)
    fit_values <- stats::predict(fit)
    dt_matrix[, i] <- comp_sample - fit_values
  }
  norm_dt <- data.table::as.data.table(dt_matrix)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Cyclic Linear Regression Normalization on MA Transformed Data
#'
#' No reference, but MA transformation and normalization
#' of samples is done pairwise between two samples with A = average of two samples and M = difference. The process is iterated through all samples pairs.
#' Log2 data should be taken as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#' @param iterations Number of cyclic iterations to be performed
#'
#' @return SummarizedExperiment containing the RlrMACyc normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- rlrMACycNorm(tuberculosis_TMT_se, ain = "log2",
#'                             aout = "RlrMACyc", on_raw = FALSE, iterations=3)
#'
rlrMACycNorm <- function(se, ain="log2", aout="RlrMACyc", on_raw = FALSE, iterations=3){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  dt_matrix <- as.matrix(dt)
  n <- ncol(dt)
  # iteration over all pairs
  for (k in seq_len(iterations)){
    for (j in seq_len(n)){
      for (i in seq(2,5)){
        # MA transformation
        ref_sample <- dt_matrix[, j]
        comp_sample <- dt_matrix[, i]
        m <- comp_sample - ref_sample
        a <- (ref_sample + comp_sample)/2
        # rlm fit
        fit <- MASS::rlm(m ~ a, na.action = stats::na.exclude)
        fit_values <- stats::predict(fit)
        # reorder values
        comp_norm <- comp_sample - fit_values/2
        ref_norm <- ref_sample + fit_values/2
        dt_matrix[, j] <- ref_norm
        dt_matrix[, i] <- comp_norm
      }
    }
  }
  norm_dt <- data.table::as.data.table(dt_matrix)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}

#' Cyclic Loess Normalization of limma
#'
#' Two samples of the data are MA transformed and normalized at a time, and all pairs of samples are iterated through. Log2-scaled data should be taken as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the loessCyc normalized data as assay (on log2 scale)
#' @export
#' @seealso \code{\link[limma]{normalizeCyclicLoess}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- loessCycNorm(tuberculosis_TMT_se, ain = "log2",
#'                                       aout = "LoessCyc", on_raw = FALSE)
#'
loessCycNorm <- function(se, ain="log2", aout="LoessCyc", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  norm_dt <- limma::normalizeCyclicLoess(as.matrix(dt), method="pairs")
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Fast Loess Normalization of limma
#'
#' Using mean intensities over all the samples as its reference A sample. Log2-scaled data should be used as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the LoessF normalized data as assay (on log2 scale)
#' @export
#' @seealso \code{\link[limma]{normalizeCyclicLoess}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- loessFNorm(tuberculosis_TMT_se, ain = "log2",
#'                                     aout = "LoessCyc", on_raw = FALSE)
#'
loessFNorm <- function(se, ain="log2", aout="LoessF", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  norm_dt <- limma::normalizeCyclicLoess(as.matrix(dt), method="fast")
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' EigenMS Normalization
#'
#' EigenMS fits an analysis of variance model to estimate the effects of the experimental factors on the data using the
#' knowledge about the experimental design, and then applies singular value decomposition to identify systematic trends 
#' contributing to significant variation not explained by the experimental factors Log2-scaled data should be used as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the EigenMS normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- eigenMSNorm(tuberculosis_TMT_se, ain = "log2",
#'                                     aout = "EigenMS", on_raw = FALSE)
#'
eigenMSNorm <- function(se, ain="log2", aout="EigenMS", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  prot.info <- cbind(data.frame(SummarizedExperiment::rowData(se)$Protein.IDs), data.frame(SummarizedExperiment::rowData(se)$Protein.IDs))
  colnames(prot.info) <- c("pepIDs", "prID")
  condition <- S4Vectors::metadata(se)$condition
  grps <- factor(data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]], levels = unique(data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]]))
  ints_eig1 <- eig_norm1(m=dt, treatment=grps, prot.info = prot.info)
  ints_norm <- eig_norm2(rv=ints_eig1)
  # present
  prot_present <- ints_eig1$present$pepIDs
  # missing
  prot_missing <- prot.info[!prot.info$pepIDs %in% prot_present,]
  prot_missing$prID <- NULL
  # complete norm_dt
  norm_dt <- data.table::as.data.table(ints_norm$normalized)
  norm_dt$prID <- NULL
  if(nrow(prot_missing)>0){
    # add these proteins with complete NAs
    for(col in names(norm_dt)){ # -1 because of pepIDs columns
      if (col!="pepIDs") prot_missing[,col] <- NA
    }
    norm_dt <- rbind(norm_dt, prot_missing)
    # take ordering of rowData
    ordering <- SummarizedExperiment::rowData(se)$Protein.IDs
    norm_dt <- norm_dt %>% dplyr::arrange(factor(pepIDs, levels = ordering))
    norm_dt$pepIDs <- NULL
  } else {
    norm_dt$pepIDs <- NULL

  }
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' Median Absolute Deviation Normalization
#'
#' Substracts the median and divides the data by the median absolute deviation (MAD). Log2-scaled data should be used as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scale data
#'
#' @return SummarizedExperiment containing the MAD normalized data as assay (on log2 scale)
#' @export
#' @seealso \code{\link[NormalyzerDE]{performSMADNormalization}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- medianAbsDevNorm(tuberculosis_TMT_se, ain = "log2",
#'                                               aout = "MAD", on_raw = FALSE)
#'
medianAbsDevNorm <- function(se, ain="log2", aout="MAD", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt <- as.matrix(dt)
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  norm_dt <- NormalyzerDE::performSMADNormalization(dt, noLogTransform = TRUE)
  norm_dt <- data.table::as.data.table(norm_dt)
  colnames(norm_dt) <- colnames(dt)
  rownames(norm_dt) <- rownames(dt)
  if(on_raw){
    norm_dt <- log2(norm_dt)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
  return(se)
}

#' RobNorm Normalization
#'
#' Log2-scaled data should be used as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#' @param gamma.0 Numeric representing the exponent of the weighted density. When the sample size
#'                is small, the fitted population of some proteins could be locally trapped such
#'                that the variance of those proteins was very small under a large gamma. To avoid
#'                this, a small gamma is recommended. When sample size smaller than 40, then set
#'                gamma to 0.5 or 0.1.
#'
#' @return SummarizedExperiment containing the RobNorm normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- robnormNorm(tuberculosis_TMT_se, ain = "log2",
#'                             aout = "RobNorm", on_raw = FALSE, gamma.0 = 0.1)
#'
robnormNorm <- function(se, ain="log2", aout="RobNorm", on_raw = FALSE, gamma.0 = 0.1){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  X.0 <- as.matrix(dt)
  rownames(X.0) <- rownames(dt)
  tryCatch(
    {
      robnorm_res <- RobNorm(X.0, gamma.0=gamma.0, tol=10^(-4), step=200) # default parameter
      norm_dt <- data.table::as.data.table(robnorm_res$norm.data)
      colnames(norm_dt) <- colnames(dt)
      rownames(norm_dt) <- rownames(dt)
      if(on_raw){
        norm_dt <- log2(norm_dt)
      }
      SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(norm_dt)
    },
    error = function(e){
      message(paste0("RobNorm: ", e$message))
    }
  )
  return(se)
}

#' limma::removeBatchEffects (limBE)
#'
#' Log2-scaled data should be used as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#'
#' @return SummarizedExperiment containing the limBE normalized data as assay (on log2 scale)
#' @export
#' @seealso \code{\link[limma]{removeBatchEffect}()}
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- limmaNorm(tuberculosis_TMT_se, ain = "log2",
#'                                       aout = "limBE", on_raw = FALSE)
#'
limmaNorm <- function(se, ain = "log2", aout = "limBE", on_raw = FALSE){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  batch <- S4Vectors::metadata(se)$batch
  batch_vector <- coldata[[batch]]
  dt_batch <- limma::removeBatchEffect(dt, batch = batch_vector)
  if(on_raw){
    dt_batch <- log2(dt_batch)
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- data.table::as.data.table(dt_batch)
  return(se)
}

#' Normics Normalization (Normics using VSN or using Median)
#'
#' Log2-scaled data should be used as input (on_raw = FALSE).
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomic dataset
#' @param ain String which assay should be used as input
#' @param aout String which assay should be used to save normalized data
#' @param method String specifying the method to use (NORMICS or NORMICSmedian)
#' @param on_raw Boolean specifying whether normalization should be performed on raw or log2-scaled data
#' @param reduce_correlation_by If the data is too big for the computation of the params, increase this parameter by 2,3,4.... The whole data will still be normalized, but the params are calculated on every second row etc.
#' @param NormicsVSN_quantile The quantile that is used for the resistant least trimmed sum of squares regression. A value of 0.8 means focusing on the central 80\% of the data, reducing the influence of outliers.
#' @param TMT_ratio Indicates if the data involves Tandem Mass Tag (TMT) ratio-based measurements (common in proteomics). If TRUE, the method may handle the data differently.
#' @param top_x Number of reference proteins extracted for the calculation of parameters
#'
#' @return SummarizedExperiment containing the NormicsVSN/NormicsMedian normalized data as assay (on log2 scale)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- normicsNorm(tuberculosis_TMT_se, ain = "raw",
#'                                 aout = "NormicsVSN", method = "NormicsVSN",
#'                                 on_raw = TRUE)
#'
#'
normicsNorm <- function(se, ain = "raw", aout = "NormicsVSN", method = "NormicsVSN", on_raw = TRUE, reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, TMT_ratio = FALSE, top_x = 50){
  # check method
  stopifnot(method %in% c("NormicsVSN", "NormicsMedian"))
  stopifnot(is.numeric(top_x))
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  # check if top_x is not bigger than the number of proteins
  stopifnot(top_x <= nrow(dt))
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  # check if ain is not raw but normalization should be performed on raw data
  if(on_raw & ain != "raw"){
    if(ain == "log2"){
      warning("Log2 data specified as ain but on_raw is set to TRUE. On_raw=TRUE forces data to be transformed to raw-scale.")
    }
    dt <- 2^dt
  }
  if(ain == "raw" & on_raw == FALSE){
    warning("Raw data specified as ain but on_raw is set to FALSE. On_raw=FALSE forces data to be transformed to log2-scale.")
    dt <- log2(dt)
  }
  dt_reduced <- dt[seq(1, nrow(dt), by = reduce_correlation_by), ]
  rowdata_reduced <- rowdata[seq(1, nrow(rowdata), by = reduce_correlation_by), ]
  # create longlist
  cols <- c("Protein_ID", "Rank_sum_RS", "Mean_correlation_MC", "Coefficient_of_variation_CV", "Variance_correlation_VC")
  # use matrices for faster computation
  dt_matrix <- as.matrix(dt_reduced)
  # precompute means and variances to avoid recalculating them
  means <- apply(dt_matrix, 1, mean, na.rm = TRUE)
  variances <- apply(dt_matrix, 1, stats::var, na.rm = TRUE)
  # initialize empty dataframe for longlist
  longlist <- data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(longlist) <- cols
  # compute pairwise correlations
  cor_matrix <- stats::cor(t(dt_matrix), method = "spearman", use = "pairwise.complete.obs")
  cor_df <- as.data.frame(cor_matrix)
  # populate longlist
  for (i in seq_len(nrow(cor_df))){
    # exclude self-correlation by removing the i-th column from the matrix
    temp <- as.numeric(cor_df[i, -i])
    CV <- ifelse(!TMT_ratio, sqrt(variances[i]) / means[i], sqrt(variances[i]))
    # add row to longlist
    longlist <- rbind(longlist,
                      data.frame( Protein_ID = rowdata_reduced$Protein.IDs[i],
                                  Rank_sum_RS = NA,
                                  Mean_correlation_MC = mean(temp, na.rm = TRUE),
                                  Coefficient_of_variation_CV = CV,
                                  Variance_correlation_VC = stats::var(temp, na.rm = TRUE)))
  }

  # sort and determine rank sum
  longlist <- longlist[order(longlist[[cols[4]]]), ]
  longlist$Rank_sum_RS <- seq_len(nrow(longlist)) - 1
  longlist <- longlist[order(-longlist[[cols[3]]]), ]
  longlist$Rank_sum_RS <- longlist$Rank_sum_RS + seq_len(nrow(longlist)) - 1
  longlist <- longlist[order(longlist[[cols[2]]]), ]

  # extract top
  shortlist <- longlist[seq_len(top_x), ]

  if(method == "NormicsVSN"){
    # NORMICS VSN normalization
    # generate subset of counts including just the top proteins
    dt_ids <- cbind(dt, rowdata[,"Protein.IDs"])
    colnames(dt_ids)[length(dt)+1] <- "Protein.IDs"
    dt_shortlist <- dt_ids[dt_ids$Protein.IDs %in% shortlist$Protein_ID,]
    dt_shortlist$Protein.IDs <- NULL
    # get parameters
    dt_to_norm <- Biobase::ExpressionSet(assayData = as.matrix(dt_shortlist))
    fit <- vsn::vsn2(dt_to_norm, lts.quantile = NormicsVSN_quantile)
    #params <- coef(fit)[1,,]
    # normalize on complete proteins data
    dt_to_norm <- Biobase::ExpressionSet(assayData = as.matrix(dt))
    norm_data <- vsn::predict(fit, dt_to_norm, log2scale = TRUE)
    norm_dt <- express_to_DT(expr_data = norm_data,
                                 column_names = colnames(dt),
                                 row_names = rownames(dt))
  } else {
    # NORMICSmedian normalization
    # calculate median ratios for shortlist proteins per column
    dt_ids <- cbind(dt, rowdata[["Protein.IDs"]])
    colnames(dt_ids)[length(dt)+1] <- "Protein.IDs"
    dt_shortlist <- dt_ids[dt_ids$Protein.IDs %in% shortlist$Protein_ID,]
    dt_shortlist$Protein.IDs <- NULL
    ratios <- as.data.frame(lapply(dt_shortlist, stats::median, na.rm = TRUE), col.names = colnames(dt_shortlist))
    # normalize on complete proteins data
    dt <- as.data.frame(dt)
    norm_dt <- dt
    # iterate over the columns and perform calculations
    for (j in seq_along(dt)) {
      ri <- ratios[1, j]
      norm_dt[, j] <- dt[, j] / ri * mean(as.numeric(ratios[1, ]), na.rm = TRUE)
    }
    # apply log2 transformation
    if(on_raw){
      norm_dt <- log2(norm_dt)
    }
    norm_dt <- tib_to_DF(norm_dt, colnames(norm_dt), rownames(norm_dt))
  }
  SummarizedExperiment::assay(se, aout, FALSE) <- norm_dt
  return(se)
}


## ----- Main Normalization Method Functions ----- ##

#' Function to return available normalization methods' identifier names
#'
#' @return Vector of normalization methods
#' @export
#'
#' @examples
#' get_normalization_methods()
#'   
#'
get_normalization_methods <- function(){
  norm_names <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                  "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "limBE", "NormicsVSN", "NormicsMedian")
  return(norm_names)
}

#' Normalize SummarizedExperiment object using different normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param methods Vector of normalization methods to apply for normalizing the proteomics data of the SummarizedExperiment object (identifier of normalization methods can be retrieved using get_normalization_methods())
#' @param on_raw Logical indicating if the normalization should be performed on the raw data or on log2-transformed data. If on_raw = NULL(default), the normalization is performed on the default method specific on_raw setting (suggestion based on publications).
#' @param gamma.0 Numeric representing the exponent of the weighted density of RobNorm normalization. When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.
#' @param reduce_correlation_by If the data is too big for the computation of the params, increase this parameter by 2,3,4.... The whole data will still be normalized, but the params are calculated on every second row etc.
#' @param NormicsVSN_quantile The quantile that is used for the resistant least trimmed sum of squares regression. A value of 0.8 means focusing on the central 80\% of the data, reducing the influence of outliers.
#' @param top_x Number of reference proteins extracted for the calculation of parameters
#' @param VSN_quantile Numeric of length 1. The quantile that is used for the resistant least trimmed sum of squares regression. (see vsn2 lts.quantile)
#'
#' @return SummarizedExperiment object with normalized data saved as assays
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- normalize_se_single(tuberculosis_TMT_se,
#'           methods = c("RobNorm", "Median", "NormicsVSN","VSN"),
#'           on_raw = NULL, gamma.0 = 0.1, reduce_correlation_by = 1,
#'           NormicsVSN_quantile = 0.8, top_x = 50, VSN_quantile = 0.9)
#'
normalize_se_single <- function(se, methods = NULL, on_raw = NULL, gamma.0 = 0.1, reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, top_x = 50, VSN_quantile = 0.9){
  # vector with available normalization methods
  norm_functions <- norm_functions <- list(globalMeanNorm, globalMedianNorm, medianNorm, meanNorm, irsNorm,
                                           quantileNorm, vsnNorm, loessFNorm, loessCycNorm, rlrNorm,
                                           rlrMANorm, rlrMACycNorm, eigenMSNorm, medianAbsDevNorm, robnormNorm, tmmNorm, limmaNorm, normicsNorm, normicsNorm)
  norm_names <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                  "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "limBE", "NormicsVSN", "NormicsMedian")
  names(norm_functions) <- norm_names

  # retrieve normalization methods & check if all methods available
  if(!is.null(methods)){
    not_available_methods <- methods[!methods %in% norm_names]
    if(length(not_available_methods) > 0){
      warning(paste0(paste0(not_available_methods, collapse = ", "), " normalization methods not available!"))
    }
  } else {
    message("All available normalization methods will be performed.")
    methods <- norm_names
  }

  # check if IRS or limBE can be executed
  batch <-S4Vectors::metadata(se)$batch
  if(sum(c("limBE", "IRS") %in% methods) > 0){
    if(is.null(batch)){
      stop("No batch specified! Batch need to be specified for limBE, IRS in the SummarizedExperiment under metadata(se)$batch!")
    }
  }
  refs <- S4Vectors::metadata(se)$refs
  if("IRS" %in% methods){
    if(is.null(refs)){
      stop("No reference samples specified! Reference samples need to be specified for IRS in the SummarizedExperiment under metadata(se)$refs!")
    }
  }

  # normalization
  for(method in methods){
    func <- norm_functions[[method]]
    if(is.null(on_raw)){
      if(method == "RobNorm"){
        se <- func(se, aout = method, gamma.0 = gamma.0)
      } else if(method %in% c("NormicsVSN", "NormicsMedian")){
        se <- func(se, aout = method, method = method, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x)
      } else if(method == "VSN"){
        se <- func(se, aout = method, VSN_quantile = VSN_quantile)
      } else {
        se <- func(se, aout = method)
      }
    } else {
      if(method == "RobNorm"){
        se <- func(se, aout = method, on_raw = on_raw, gamma.0 = gamma.0)
      } else if(method %in% c("NormicsVSN", "NormicsMedian")){
        se <- func(se, aout = method, on_raw = on_raw, method = method, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x)
      } else if(method == "VSN"){
        se <- func(se, aout = method, on_raw = on_raw, VSN_quantile = VSN_quantile)
      } else {
        se <- func(se, aout = method, on_raw = on_raw)
      }
    }


    # TODO: error handling
    message(paste0(method, " completed."))
  }
  return(se)
}


#' Normalize SummarizedExperiment object using combinations of normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param methods Vector of normalization methods to apply for normalizing the proteomics data of the SummarizedExperiment object (identifier of normalization methods can be retrieved using get_normalization_methods())
#' @param ains Vector of assays of SummarizedExperiment object to apply the normalization methods (e.g. if you want to perform Median normalization on IRS-normalized data)
#' @param on_raw Logical indicating if the normalization should be performed on the raw data or on log2-transformed data. If on_raw = NULL(default), the normalization is performed on the default method specific on_raw setting (suggestion based on publications).
#' @param combination_pattern String to give name to combination of methods (e.g. IRS_on_Median --> "_on_")
#' @param gamma.0 Numeric representing the exponent of the weighted density of RobNorm normalization. When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.
#' @param reduce_correlation_by If the data is too big for the computation of the params, increase this parameter by 2,3,4.... The whole data will still be normalized, but the params are calculated on every second row etc.
#' @param NormicsVSN_quantile The quantile that is used for the resistant least trimmed sum of squares regression. A value of 0.8 means focusing on the central 80\% of the data, reducing the influence of outliers.
#' @param top_x Number of reference proteins extracted for the calculation of parameters
#' @param VSN_quantile Numeric of length 1. The quantile that is used for the resistant least trimmed sum of squares regression. (see vsn2 lts.quantile)
#'
#' @return SummarizedExperiment object with normalized data saved as assays
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- normalize_se_combination(tuberculosis_TMT_se,
#'           methods = c("Median","NormicsVSN"), ains = c("IRS"), on_raw = NULL,
#'           combination_pattern = "_on_", gamma.0 = 0.1,
#'           reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, top_x = 50,
#'           VSN_quantile = 0.9)
#'
normalize_se_combination <- function(se, methods, ains, on_raw = NULL, combination_pattern = "_on_", gamma.0 = 0.1, reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, top_x = 50, VSN_quantile = 0.9){

  # vector with available normalization methods
  norm_functions <- norm_functions <- list(globalMeanNorm, globalMedianNorm, medianNorm, meanNorm, irsNorm,
                                           quantileNorm, vsnNorm, loessFNorm, loessCycNorm, rlrNorm,
                                           rlrMANorm, rlrMACycNorm, eigenMSNorm, medianAbsDevNorm, robnormNorm, tmmNorm, normicsNorm, normicsNorm, limmaNorm)
  norm_names <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                  "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsVSN", "NormicsMedian", "limBE")
  names(norm_functions) <- norm_names

  # retrieve normalization methods & check if all methods available
  if(!is.null(methods)){
    # check if all methods available
    not_available_methods <- methods[!methods %in% norm_names]
    if(length(not_available_methods) > 0){
      warning(paste0(paste0(not_available_methods, collapse = ", "), " normalization methods not available!"))
    }
  } else {
    message("All available normalization methods will be performed.")
    methods <- norm_names
  }

  # perform combination of normalization
  for(ain in ains){
    # check if ain already in se --> if not: perform now
    if(! ain %in% names(SummarizedExperiment::assays(se))){
      message(paste0(ain, " normalization not yet performed. Single ", ain, " normalization performed now."))
      se <- normalize_se_single(se, methods = c(ain), on_raw = on_raw, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x)
    }

    for(method in methods){
      aout <- paste0(method, combination_pattern, ain)
      func <- norm_functions[[method]]
      if(is.null(on_raw)){
        if(method == "RobNorm"){
          se <- func(se, ain = ain, aout = aout, gamma.0 = gamma.0)
        } else if(method %in% c("NormicsVSN", "NormicsMedian")){
          se <- func(se, ain = ain, aout = aout, method = method, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x)
        } else if(method == "VSN"){
          se <- func(se, ain = ain,aout = aout, VSN_quantile = VSN_quantile)
        } else {
          se <- func(se, ain = ain, aout = aout)
        }
      } else {
        if(method == "RobNorm"){
          se <- func(se, ain = ain, aout = aout, on_raw = on_raw, gamma.0 = gamma.0)
        } else if(method %in% c("NormicsVSN", "NormicsMedian")){
          se <- func(se, ain = ain, aout = aout, on_raw = on_raw, method = method, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x)
        } else if(method == "VSN"){
          se <- func(se, ain = ain,aout = aout, on_raw = on_raw, VSN_quantile = VSN_quantile)
        } else {
          se <- func(se, ain = ain, aout = aout, on_raw = on_raw)
        }
      }

      message(paste0(method, " normalization performed on ", ain, "-normalized data completed."))
    }
  }
  return(se)
}


#' Normalize SummarizedExperiment object using single normalization methods or specified combinations of normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param methods Vector of normalization methods to apply for normalizing the proteomics data of the SummarizedExperiment object (identifier of normalization methods can be retrieved using get_normalization_methods())
#' @param combination_pattern String specifying how normalization methods are combined. For instance, methods = c("IRS", "Median_on_IRS"), combination_pattern = "_on_".
#' @param on_raw Logical indicating if the normalization should be performed on the raw data or on log2-transformed data. If on_raw = NULL(default), the normalization is performed on the default method specific on_raw setting (suggestion based on publications).
#' @param gamma.0 Numeric representing the exponent of the weighted density of RobNorm normalization. When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.
#' @param reduce_correlation_by If the data is too big for the computation of the params, increase this parameter by 2,3,4.... The whole data will still be normalized, but the params are calculated on every second row etc.
#' @param NormicsVSN_quantile The quantile that is used for the resistant least trimmed sum of squares regression. A value of 0.8 means focusing on the central 80\% of the data, reducing the influence of outliers.
#' @param top_x Number of reference proteins extracted for the calculation of parameters
#' @param VSN_quantile Numeric of length 1. The quantile that is used for the resistant least trimmed sum of squares regression (see vsn2 lts.quantile)
#'
#' @return SummarizedExperiment object with normalized data saved as assays
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- normalize_se(tuberculosis_TMT_se,
#'            methods = c("IRS_on_GlobalMedian", "IRS_on_Median",
#'            "limBE_on_NormicsVSN"), on_raw = NULL,
#'            combination_pattern = "_on_", gamma.0 = 0.1,
#'            reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, top_x = 50,
#'            VSN_quantile = 0.9)
#'
normalize_se <- function(se, methods, combination_pattern = "_on_", on_raw = NULL, gamma.0 = 0.1, reduce_correlation_by = 1, NormicsVSN_quantile = 0.8, top_x = 50, VSN_quantile = 0.9){
  # extract combination of methods
  if(!is.null(combination_pattern)){
    comb_methods <- methods[stringr::str_detect(methods, combination_pattern)] #  combined methods
    sing_methods <- methods[!methods %in% comb_methods] # single methods
  } else {
    sing_methods <- methods
  }

  # single normalization
  if(length(sing_methods) > 0){
    se <- normalize_se_single(se, sing_methods, on_raw = on_raw, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x, VSN_quantile = VSN_quantile)
  }
  # combined normalization
  if(!is.null(combination_pattern)){
    if(length(comb_methods) > 0){
      for(m in comb_methods){
        methods_split <- strsplit(m, combination_pattern)[[1]]

        # Initialize method and ain with the last two elements
        length_split <- length(methods_split)
        ain <- methods_split[length_split]
        method <- methods_split[length_split - 1]
        se <- normalize_se_combination(se, c(method), c(ain), on_raw = on_raw, combination_pattern, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x, VSN_quantile = VSN_quantile)

        # If there are more than two methods in the combination, process the rest
        if (length(methods_split) > 2) {
          for (i in (length_split - 2):1) {
            ain <- paste(method, ain, sep = combination_pattern)
            method <- methods_split[i]

            # Apply the subsequent combination
            se <- normalize_se_combination(se, c(method), c(ain), on_raw = on_raw, combination_pattern, gamma.0 = gamma.0, reduce_correlation_by = reduce_correlation_by, NormicsVSN_quantile = NormicsVSN_quantile, top_x = top_x, VSN_quantile = VSN_quantile)
          }
        }
      }
    }
  }
  return(se)
}




