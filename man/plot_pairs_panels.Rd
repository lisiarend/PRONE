% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetOverview.R
\name{plot_pairs_panels}
\alias{plot_pairs_panels}
\title{Plot correlation, histogram, and scatterplot of samples of a specific group for a selection of normalization methods}
\usage{
plot_pairs_panels(se, ain = NULL, condition = NULL, label_by = NULL)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{ain}{Vector of strings which assay should be used as input (default NULL).
If NULL then all normalization of the se object are plotted next to each other.}

\item{condition}{Column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)}

\item{label_by}{String specifying the column in the metadata used to label the samples for the UpSet plot}
}
\value{
list of ggplot objects
}
\description{
Plot correlation, histogram, and scatterplot of samples of a specific group for a selection of normalization methods
}
\examples{
data(tuberculosis_TMT_se)
plot_pairs_panels(tuberculosis_TMT_se, ain = c("log2", "IRS_on_RobNorm"),
                  condition = NULL, label_by = NULL)

}
