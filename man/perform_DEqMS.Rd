% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEAnalysis.R
\name{perform_DEqMS}
\alias{perform_DEqMS}
\title{Perform DEqMS}
\usage{
perform_DEqMS(
  fit,
  se,
  DEqMS_PSMs_column = NULL,
  logFC = TRUE,
  logFC_up = 1,
  logFC_down = -1,
  p_adj = TRUE,
  alpha = 0.05
)
}
\arguments{
\item{fit}{eBayes object resulting from perform_limma method}

\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{DEqMS_PSMs_column}{String specifying which column name to use for DEqMS (default NULL). Any column of the rowData(se) is accepted.}

\item{logFC}{Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)}

\item{logFC_up}{Upper log2 fold change threshold (dividing into up regulated)}

\item{logFC_down}{Lower log2 fold change threshold (dividing into down regulated)}

\item{p_adj}{Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)}

\item{alpha}{Threshold for adjusted p-values or p-values}
}
\value{
data.table of DE results
}
\description{
Perform DEqMS
}
