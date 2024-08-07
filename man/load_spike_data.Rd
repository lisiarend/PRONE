% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadData.R
\name{load_spike_data}
\alias{load_spike_data}
\title{Load spike-in proteomics data into a SummarizedExperiment}
\usage{
load_spike_data(
  data,
  md,
  spike_column,
  spike_value,
  spike_concentration,
  protein_column = "Protein.IDs",
  gene_column = "Gene.Names",
  ref_samples = NULL,
  batch_column = NULL,
  condition_column = NULL,
  label_column = NULL
)
}
\arguments{
\item{data}{tabular data table with rows = proteins and columns = samples (such as proteinGroups.txt of MaxQuant)}

\item{md}{experimental design table (requires a column named "Column" for the column names of the sample intensities in data)}

\item{spike_column}{name of the column specifying which proteins are the spike-ins}

\item{spike_value}{String value specifying the spike-in proteins in the spike-in column}

\item{spike_concentration}{name of the column in md defining the spike-in concentration per sample}

\item{protein_column}{name of the column in data containing the protein IDs}

\item{gene_column}{name of the column in data containing the gene names}

\item{ref_samples}{reference samples if TMT experiment provided (names of samples)}

\item{batch_column}{name of the column in md defining the batches}

\item{condition_column}{name of the column in md defining the condition (can still be changed afterwards)}

\item{label_column}{name of the column in md containing simple sample names (for visualization)}
}
\value{
SummarizedExperiment object
}
\description{
Load spike-in proteomics data into a SummarizedExperiment
}
\examples{
data_path <- readPRONE_example("Ecoli_human_MaxLFQ_protein_intensities.csv")
md_path <- readPRONE_example("Ecoli_human_MaxLFQ_metadata.csv")
data <- read.csv(data_path)
md <- read.csv(md_path)
mixed <- grepl("Homo sapiens.*Escherichia|Escherichia.*Homo sapiens", data$Fasta.headers)
data <- data[!mixed,]
data$Spiked <- rep("HUMAN", nrow(data))
data$Spiked[grepl("ECOLI", data$Fasta.headers)] <- "ECOLI"
se <- load_spike_data(data, md, spike_column = "Spiked", spike_value = "ECOLI", 
        spike_concentration = "Concentration", protein_column = "Protein.IDs", 
        gene_column = "Gene.names", ref_samples = NULL, batch_column = NULL, 
        condition_column = "Condition", label_column = "Label")

}
