% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NormalizationMethods.R
\name{meanNorm}
\alias{meanNorm}
\title{Mean Normalization}
\usage{
meanNorm(se, ain = "raw", aout = "Mean", on_raw = TRUE)
}
\arguments{
\item{se}{SummarizedExperiment containing all necessary information of the proteomic dataset}

\item{ain}{String which assay should be used as input}

\item{aout}{String which assay should be used to save normalized data}

\item{on_raw}{Boolean whether normalized should be performed on raw or log2-scaled data}
}
\value{
SummarizedExperiment containing the mean normalized data as assay (on log2 scale)
}
\description{
The intensity of each protein group in a given sample is divided by the mean of the
intensities of all protein groups in that sample and then multiplied by the mean of
mean of sum of intensities of all protein groups in all samples.
}
\examples{
data(tuberculosis_TMT_se)
tuberculosis_TMT_se <- meanNorm(tuberculosis_TMT_se, ain = "raw",
                          aout = "Mean", on_raw = TRUE)

}
