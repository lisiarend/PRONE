% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpikeDEPlots.R
\name{plot_logFC_thresholds_spiked}
\alias{plot_logFC_thresholds_spiked}
\title{Line plot of number of true and false positives when applying different logFC thresholds}
\usage{
plot_logFC_thresholds_spiked(
  se,
  de_res,
  condition,
  ain = NULL,
  comparisons = NULL,
  nrow = 2,
  alpha = 0.05
)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{de_res}{data table resulting of run_DE}

\item{condition}{column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)}

\item{ain}{Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)}

\item{comparisons}{Vector of comparisons (must be valid comparisons saved in stats)}

\item{nrow}{number of rows for facet wrap}

\item{alpha}{threshold for adjusted p-values}
}
\value{
list of ggplot objects
}
\description{
Line plot of number of true and false positives when applying different logFC thresholds
}
\examples{
data(spike_in_se)
data(spike_in_de_res)
plot_logFC_thresholds_spiked(spike_in_se, spike_in_de_res,
                             condition = "Condition", ain = NULL,
                             comparisons = NULL, nrow = 2, alpha = 0.05)

}
