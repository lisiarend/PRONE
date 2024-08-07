% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpikeDEPlots.R
\name{plot_ROC_AUC_spiked}
\alias{plot_ROC_AUC_spiked}
\title{Plot ROC curve and barplot of AUC values for each method for a specific comparion or for all comparisons}
\usage{
plot_ROC_AUC_spiked(se, de_res, ain = NULL, comparisons = NULL)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{de_res}{data table resulting of run_DE}

\item{ain}{Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)}

\item{comparisons}{Vector of comparisons (must be valid comparisons saved in stats)}
}
\value{
list of ggplot objects
}
\description{
Plot ROC curve and barplot of AUC values for each method for a specific comparion or for all comparisons
}
\examples{
data(spike_in_se)
data(spike_in_de_res)
plot_ROC_AUC_spiked(spike_in_se, spike_in_de_res)

}
