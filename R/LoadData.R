
#' Load real-world proteomics data into a SummarizedExperiment
#'
#' @param data tabular data table with rows = proteins and columns = samples (such as proteinGroups.txt of MaxQuant)
#' @param md experimental design table (requires a column named "Column" for the column names of the sample intensities in data)
#' @param protein_column name of the column in data containing the protein IDs
#' @param gene_column name of the column in data containing the gene names
#' @param ref_samples reference samples if TMT experiment provided (names of samples)
#' @param batch_column name of the column in md defining the batches
#' @param condition_column name of the column in md defining the condition (can still be changed afterwards)
#' @param label_column name of the column in md containing simple sample names (for visualization)
#' @importFrom magrittr %>%
#'
#' @return SummarizedExperiment object
#' @export
#'
#' @examples
#' data_path <- readPRONE_example("tuberculosis_protein_intensities.csv")
#' md_path <- readPRONE_example("tuberculosis_metadata.csv")
#' data <- read.csv(data_path)
#' md <- read.csv(md_path)
#' md$Column <- stringr::str_replace_all(md$Column, " ", ".")
#' ref_samples <- md[md$Group == "ref",]$Column
#' se <- load_data(data, md, protein_column = "Protein.IDs",
#'                gene_column = "Gene.names", ref_samples = ref_samples,
#'                batch_column = "Pool", condition_column = "Group",
#'                label_column = "Label")
#'
load_data <- function(data, md, protein_column = "Protein.IDs", gene_column = "Gene.Names", ref_samples = NULL, batch_column = NULL, condition_column = NULL, label_column = NULL){
  # convert to data.table
  data <- data.table::as.data.table(data)
  md <- data.table::as.data.table(md)
  # assays
  raw <- data[,md$Column, with=FALSE]
  raw[raw == 0] <- NA
  log2 <- log2(raw)
  # rowData
  rowData <- data[, !colnames(data) %in% md$Column, with = FALSE]
  rowData$IDs <- rownames(raw)
  rowData <- rowData %>% dplyr::rename("Protein.IDs" = protein_column, "Gene.Names" = gene_column)
  # colData
  colData <- md

  # check if ref samples in data
  if(!is.null(ref_samples)){
    if(length(c(ref_samples) %in% colData$Column) != length(c(ref_samples))){
      stop("Reference samples not all included in data. Aborted!")
    }
  }
  # check if batch column in data
  if(!is.null(batch_column)){
    if(!batch_column %in% colnames(colData)){
      stop("Batch column not in the data. Aborted!")
    }
  }
  # check if label column in the data
  if(!is.null(label_column)){
    if(!label_column %in% colnames(colData)){
      stop("Label column not in the data. Aborted!")
    }
  }
  # check if condition column in the data
  if(!is.null(condition_column)){
    if(!condition_column %in% colnames(colData)){
      stop("Condition column not in the data. Aborted!")
    }
  }

  # metadata
  metadata <- list(condition = condition_column, batch = batch_column, refs = ref_samples, label = label_column)
  # SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = raw, log2 = log2), colData = colData, rowData = rowData, metadata = metadata)
  return(se)
}


#' Load spike-in proteomics data into a SummarizedExperiment
#'
#' @param data tabular data table with rows = proteins and columns = samples (such as proteinGroups.txt of MaxQuant)
#' @param md experimental design table (requires a column named "Column" for the column names of the sample intensities in data)
#' @param spike_column name of the column specifying which proteins are the spike-ins
#' @param spike_value String value specifying the spike-in proteins in the spike-in column
#' @param spike_concentration name of the column in md defining the spike-in concentration per sample
#' @param protein_column name of the column in data containing the protein IDs
#' @param gene_column name of the column in data containing the gene names
#' @param ref_samples reference samples if TMT experiment provided (names of samples)
#' @param batch_column name of the column in md defining the batches
#' @param condition_column name of the column in md defining the condition (can still be changed afterwards)
#' @param label_column name of the column in md containing simple sample names (for visualization)
#' @importFrom magrittr %>%
#'
#' @return SummarizedExperiment object
#' @export
#'
#' @examples
#' data_path <- readPRONE_example("Ecoli_human_MaxLFQ_protein_intensities.csv")
#' md_path <- readPRONE_example("Ecoli_human_MaxLFQ_metadata.csv")
#' data <- read.csv(data_path)
#' md <- read.csv(md_path)
#' mixed <- grepl("Homo sapiens.*Escherichia|Escherichia.*Homo sapiens", data$Fasta.headers)
#' data <- data[!mixed,]
#' data$Spiked <- rep("HUMAN", nrow(data))
#' data$Spiked[grepl("ECOLI", data$Fasta.headers)] <- "ECOLI"
#' se <- load_spike_data(data, md, spike_column = "Spiked", spike_value = "ECOLI", 
#'         spike_concentration = "Concentration", protein_column = "Protein.IDs", 
#'         gene_column = "Gene.names", ref_samples = NULL, batch_column = NULL, 
#'         condition_column = "Condition", label_column = "Label")
#'
load_spike_data <- function(data, md, spike_column, spike_value, spike_concentration, protein_column = "Protein.IDs", gene_column = "Gene.Names", ref_samples = NULL, batch_column = NULL, condition_column = NULL, label_column = NULL){
  # convert to data.table
  data <- data.table::as.data.table(data)
  md <- data.table::as.data.table(md)
  # assays
  raw <- data[,md$Column, with=FALSE]
  raw[raw == 0] <- NA
  log2 <- log2(raw)
  # rowData
  rowData <- data[, !colnames(data) %in% md$Column, with = FALSE]
  rowData$IDs <- rownames(raw)
  rowData <- rowData %>% dplyr::rename("Protein.IDs" = protein_column, "Gene.Names" = gene_column)
  # colData
  colData <- md

  # check if ref samples in data
  if(!is.null(ref_samples)){
    if(length(c(ref_samples) %in% colData$Column) != length(c(ref_samples))){
      stop("Reference samples not all included in data. Aborted!")
    }
  }
  # check if batch column in data
  if(!is.null(batch_column)){
    if(!batch_column %in% colnames(colData)){
      stop("Batch column not in the data. Aborted!")
    }
  }
  # check if label column in the data
  if(!is.null(label_column)){
    if(!label_column %in% colnames(colData)){
      stop("Label column not in the data. Aborted!")
    }
  }
  # check if condition column in the data
  if(!is.null(condition_column)){
    if(!condition_column %in% colnames(colData)){
      stop("Condition column not in the data. Aborted!")
    }
  }

  # check if spike column in the data
  if(is.null(spike_column)){
    stop("Spike column need to be specified!")
  } else {
    if(!spike_column %in% colnames(rowData)){
      stop("Spike column not in the data. Aborted!")
    } else {
      # check if spike_value in data
      sc <- rowData[[spike_column]]
      if(!spike_value %in% sc){
        stop("Spike value not in spike column. Aborted!")
      }
    }
  }

  # check if spike_concentration column in data
  if(!is.null(spike_concentration)){
    if(!spike_concentration %in% colnames(colData)){
      stop("Spike concentration column not in the data. Aborted!")
    }
  }

  # metadata
  metadata <- list(condition = condition_column, batch = batch_column, refs = ref_samples, label = label_column, spike_column = spike_column, spike_value = spike_value, spike_concentration = spike_concentration)
  # SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = raw, log2 = log2), colData = colData, rowData = rowData, metadata = metadata)
  return(se)
}
