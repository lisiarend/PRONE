% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FilterProteins.R
\name{plot_NA_frequency}
\alias{plot_NA_frequency}
\title{Plot Protein identification overlap (x = Identified in Number of Samples, y=Number of Proteins)}
\usage{
plot_NA_frequency(se)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}
}
\value{
ggplot object
}
\description{
Plot Protein identification overlap (x = Identified in Number of Samples, y=Number of Proteins)
}
\examples{
data(tuberculosis_TMT_se)
plot_NA_frequency(tuberculosis_TMT_se)

}
