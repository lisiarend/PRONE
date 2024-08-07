% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEPlots.R
\name{plot_overview_DE_tile}
\alias{plot_overview_DE_tile}
\title{Overview heatmap plot of DE results}
\usage{
plot_overview_DE_tile(de_res, ain = NULL, comparisons = NULL)
}
\arguments{
\item{de_res}{data table resulting of run_DE}

\item{ain}{Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)}

\item{comparisons}{Vector of comparisons (must be valid comparisons saved in de_res)}
}
\value{
ggplot object
}
\description{
Overview heatmap plot of DE results
}
\examples{
data(tuberculosis_TMT_de_res)
plot_overview_DE_tile(tuberculosis_TMT_de_res, ain = NULL,
                      comparisons = NULL)

}
