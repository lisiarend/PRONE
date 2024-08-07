% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEPlots.R
\name{plot_coverage_DE_markers}
\alias{plot_coverage_DE_markers}
\title{Barplot of coverage of DE markers per normalization method in any comparison. (If you want to have a look at a specific comparison, just subset the de_res data table before plotting.)}
\usage{
plot_coverage_DE_markers(se, de_res, ain, markers, id_column = "Protein.IDs")
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{de_res}{data table resulting of run_DE}

\item{ain}{Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)}

\item{markers}{vector of the IDs of the markers to plot}

\item{id_column}{String specifying the column of the rowData of the SummarizedExperiment object which includes the IDs of the markers}
}
\value{
ggplot object
}
\description{
Barplot of coverage of DE markers per normalization method in any comparison. (If you want to have a look at a specific comparison, just subset the de_res data table before plotting.)
}
\examples{
data(tuberculosis_TMT_se)
data(tuberculosis_TMT_de_res)
plot_coverage_DE_markers(tuberculosis_TMT_se, tuberculosis_TMT_de_res,
                          ain = NULL, markers = c("Q7Z7F0", "Q13790"),
                          id_column = "Protein.IDs")

}
