% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FilterSamples.R
\name{detect_outliers_POMA}
\alias{detect_outliers_POMA}
\title{Outlier detection via POMA R Package}
\usage{
detect_outliers_POMA(
  se,
  ain = "log2",
  condition = NULL,
  method = "euclidean",
  type = "median",
  group = TRUE,
  coeff = 1.5
)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{ain}{String which data type should be used (default raw)}

\item{condition}{Column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)}

\item{method}{String specifying the method that should be used to calculate the distance matrix}

\item{type}{String specifying the type of distance calculation to centroid or spatial median}

\item{group}{String specifying if the outlier detection should be performed multi-variate
(with conditions) or on the complete data set}

\item{coeff}{This value corresponds to the classical 1.5 in Q3 + 1.5 * IQR formula to detect outliers. By changing this value, the permissiveness in outlier detection will change.}
}
\value{
list of two ggplot objects and a data.table with outlier samples
}
\description{
Outlier detection via POMA R Package
}
\examples{
data(tuberculosis_TMT_se)
poma_res <- detect_outliers_POMA(tuberculosis_TMT_se, ain="raw",
                                 condition = NULL, method="euclidean",
                                 type="median", group=TRUE, coeff = 1.5)
}
