% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEPlots.R
\name{plot_heatmap_DE}
\alias{plot_heatmap_DE}
\title{Heatmap of DE results}
\usage{
plot_heatmap_DE(
  se,
  de_res,
  ain,
  comparison,
  condition = NULL,
  label_by = NULL,
  pvalue_column = "adj.P.Val",
  col_vector = NULL
)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set (including the normalized intensities)}

\item{de_res}{data table resulting of run_DE}

\item{ain}{Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)}

\item{comparison}{String of comparison (must be a valid comparison saved in de_res)}

\item{condition}{column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)}

\item{label_by}{String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)}

\item{pvalue_column}{column name of p-values in de_res}

\item{col_vector}{Vector of colors to use for the heatmap. If NULL, default colors are used.}
}
\value{
list of ComplexHeatmaps for each method
}
\description{
Heatmap of DE results
}
\examples{
data(tuberculosis_TMT_se)
data(tuberculosis_TMT_de_res)
plot_heatmap_DE(tuberculosis_TMT_se, tuberculosis_TMT_de_res, ain = NULL,
                comparison = "PTB-HC",
                condition = NULL, label_by = NULL,
                pvalue_column = "adj.P.Val", col_vector = NULL)

}
