% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Enrichment.R
\name{plot_intersection_enrichment}
\alias{plot_intersection_enrichment}
\title{Intersect top N enrichment terms per normalization method}
\usage{
plot_intersection_enrichment(
  se,
  de_res,
  ain = NULL,
  comparisons = NULL,
  id_column = "Gene.Names",
  organism = "hsapiens",
  per_comparison = TRUE,
  sources = c("GO:BP", "GO:MF", "GO:CC"),
  top = 10
)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomics data set}

\item{de_res}{data table resulting of run_DE}

\item{ain}{Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)}

\item{comparisons}{Vector of comparisons (must be valid comparisons saved in de_res)}

\item{id_column}{String specifying the column of the rowData of the SummarizedExperiment object which includes the gene names}

\item{organism}{Organism name (gprofiler parameter)}

\item{per_comparison}{Boolean specifying whether the enrichment analysis should be performed per comparison (TRUE) or on all given comparisons together (FALSE)}

\item{sources}{Vector of data sources to use (gprofiler parameter)}

\item{top}{Number of enrichment terms to extract for each normalization method}
}
\value{
list of ggplot objects or single ggplot object
}
\description{
Intersect top N enrichment terms per normalization method
}
\examples{
data(tuberculosis_TMT_se)
data(tuberculosis_TMT_de_res)
plot_intersection_enrichment(tuberculosis_TMT_se, tuberculosis_TMT_de_res,
                ain = c("IRS_on_RobNorm", "IRS_on_Median"),
                comparisons = NULL, id_column = "Gene.Names",
                organism = "hsapiens", per_comparison = TRUE,
                sources = c("GO:BP", "GO:MF", "GO:CC"), top = 10)

}
