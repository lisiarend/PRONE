################# ----------- Filter Proteins Functions ----------- ###################

#' Remove proteins with NAs in all samples
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberulosis_TMT_se <- filter_out_complete_NA_proteins(tuberculosis_TMT_se)
#'
filter_out_complete_NA_proteins <- function(se){
  # find proteins that only have NAs
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])
  row_idx_no_NA <- rowSums(is.na(dt)) != ncol(dt)
  # subset SummarizedExperiment
  se_subset <- se[row_idx_no_NA,]
  if(nrow(se_subset) == 0){
    stop("All proteins would be removed. Aborted!")
  } else if(nrow(se_subset) == nrow(se)){
    message("No proteins were removed.")
  } else{
    rm <- nrow(se) - nrow(se_subset)
    message(paste0(rm, " proteins were removed."))
  }
  return(se_subset)
}


#' Remove proteins by value in specific column
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param column_name name of column of which proteins with a specific value should be removed
#' @param values value of the column defining the proteins that should be removed
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- filter_out_proteins_by_value(tuberculosis_TMT_se,
#'                                  column_name = "Reverse", values = c("+"))
#'
filter_out_proteins_by_value <- function(se, column_name = "Reverse", values = c("+")){
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  if(column_name %in% colnames(rd)){
    rd_subset <- rd[! rd[, column_name] %in% c(values),]
    se_subset <- se[rd_subset$ID, ]
    if(nrow(se_subset) == 0){
      stop("All proteins would be removed. Aborted!")
    } else if(nrow(se_subset) == nrow(se)){
      message("No proteins were removed.")
    } else{
      rm <- nrow(se) - nrow(se_subset)
      message(paste0(rm, " proteins were removed."))
    }
  } else {
    warning(paste0("Column ", column_name, " not in SummarizedExperiment object!"))
    se_subset <- se
  }
  return(se_subset)
}

#' Get proteins by value in specific column
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param column_name name of column of which proteins with a specific value should be identified
#' @param values value of the column defining the proteins that should be identified
#'
#' @return vector of protein IDs
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' proteins <- get_proteins_by_value(tuberculosis_TMT_se,
#'                       column_name = "Potential.contaminant", values = c("+"))
#'
get_proteins_by_value <- function(se, column_name = "Reverse", values = c("+")){
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  if(column_name %in% colnames(rd)){
    rd_subset <- rd[rd[[column_name]] %in% c(values),]
    if(nrow(rd_subset) == nrow(rd)){
      message("All proteins were identified.")
    } else if(nrow(rd_subset) == 0){
      message("No proteins were identified.")
    } else{
      message(paste0(nrow(rd_subset), " proteins were identified."))
    }
    return(rd_subset$Protein.IDs)
  } else {
    warning(paste0("Column ", column_name, " not in SummarizedExperiment object!"))
    return(c())
  }
}

#' Remove proteins by their ID
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param protein_ids Vector of protein IDs that should be kept
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- filter_out_proteins_by_ID(tuberculosis_TMT_se,
#'                                 protein_ids = c("P0A8V2", "P0A8V2"))
#'
filter_out_proteins_by_ID <- function(se, protein_ids){
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  rd_subset <- rd[!rd$Protein.IDs %in% protein_ids,]
  se_subset <- se[rd_subset$ID, ]
  if(nrow(se_subset) == 0){
    stop("No proteins would be remaining. Aborted!")
  } else if(nrow(se_subset) == nrow(se)){
    message("No proteins were removed.")
  } else{
    rm <- nrow(se) - nrow(se_subset)
    message(paste0(rm, " proteins were removed."))
  }
  return(se_subset)
}

#' Filter proteins based on their NA pattern using a specific threshold
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param thr Threshold for the minimum fraction of valid values allowed for any protein
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- filter_out_NA_proteins_by_threshold(tuberculosis_TMT_se,
#'                                                        thr = 0.8)
#'
filter_out_NA_proteins_by_threshold <- function(se, thr = 0.8){
  # prepare data
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])
  dt$ID <- SummarizedExperiment::rowData(se)$Protein.IDs
  melted_dt <- data.table::melt(dt, variable.name = "Column", value.name = "Intensity", id.vars = "ID")
  melted_dt <- merge(melted_dt, data.table::as.data.table(SummarizedExperiment::colData(se)),by="Column")

  # make values binary (1 = valid value, 0 = missing value)
  bin_dt <- melted_dt
  bin_dt$Intensity[is.na(bin_dt$Intensity)] <- 0
  bin_dt$Intensity[bin_dt$Intensity != 0] <- 1

  # filter data
  nr_samples <- dim(se)[2]
  keep <- bin_dt %>% dplyr::group_by(ID) %>% dplyr::summarise(miss_val = dplyr::n() - sum(Intensity), valid_val = sum(Intensity)) %>% dplyr::filter(valid_val >= thr * nr_samples)

  # get indices of keep proteins
  rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  rowdata <- rowdata[rowdata$Protein.IDs %in% keep$ID, ]

  se_subset <- se[rowdata$IDs,]
  rm <- nrow(se) - nrow(se_subset)
  if(nrow(se_subset) == 0){
    stop("All proteins would be removed. Aborted!")
  } else if (rm == 0){
    message("No proteins were removed.")
  } else{
    message(paste0(rm, " proteins were removed."))
  }
  return(se_subset)
}


################# ----------- Plotting Functions ----------- ###################

#' Plot heatmap of the NA pattern
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#' @param cluster_samples Boolean. TRUE if samples should be clustered, else FALSE.
#' @param cluster_proteins Boolean. TRUE if proteins should be clustered, else FALSE.
#' @param show_row_dend Boolean. TRUE if row dendrogram should be shown.
#' @param show_column_dend Boolean. TRUE if column dendrogram should be shown.
#' @importFrom magrittr %>%
#'
#' @return ComplexHeatmap plot (only showing proteins with at least one missing value)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_NA_heatmap(tuberculosis_TMT_se, color_by = NULL,
#'                 label_by = NULL, cluster_samples = TRUE,
#'                 cluster_proteins = TRUE, show_row_dend = TRUE,
#'                 show_column_dend = FALSE)
#'
plot_NA_heatmap <- function(se, color_by = NULL, label_by = NULL, cluster_samples = TRUE, cluster_proteins = TRUE, show_row_dend = TRUE, show_column_dend = FALSE){
  # get color and label values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # prepare data --> make values binary (1 = valid value, 0 = missing value)
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])
  dt[is.na(dt)] <- 0
  dt[dt != 0] <- 1

  # only take those proteins with at least one missing value
  dt <- dt[apply(dt, 1, function(x) any(x==0)),]

  # preserve order
  dt <- dt[, SummarizedExperiment::colData(se)$Column, with=FALSE]
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  colnames(dt) <- coldata$Column

  if(show_sample_names){
    colnames(dt) <- as.character(coldata[, label_by])
  } else {
    colnames(dt) <- coldata$Column
  }

  # color vector
  if(!is.null(color_by)){
    # prepare data for color bar
    annotation_col <- data.table::data.table(color_by = coldata[, color_by])
    colnames(annotation_col) <- c(color_by)

    if(show_sample_names){
      rownames(annotation_col) <- coldata[, label_by]
    } else {
      rownames(annotation_col) <- coldata$Column
    }

    qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- rev(col_vector)

    condition_colors <- col_vector[seq_len(length(unique(annotation_col[, color_by])))]
    names(condition_colors) <- unique(annotation_col[, color_by])
    annotation_colors <- list(color_by = condition_colors)
    ha <- ComplexHeatmap::HeatmapAnnotation(color_by = annotation_col[, color_by], col = annotation_colors, annotation_label = c(color_by))
    # Plot Heatmap
    p <- ComplexHeatmap::Heatmap(as.matrix(dt),
                                 name = "Value Type",
                                 cluster_rows = cluster_proteins,
                                 cluster_columns = cluster_samples,
                                 col = c("black", "grey90"),
                                 show_row_dend = show_row_dend,
                                 show_column_dend = show_column_dend,
                                 show_column_names = show_sample_names,
                                 top_annotation = ha,
                                 heatmap_legend_param = list(at = c(0,1),labels = c("Missing", "Valid"))
    )
  } else {
    # no color bar
    p <- ComplexHeatmap::Heatmap(as.matrix(dt),
                                 name = "Value Type",
                                 cluster_rows = cluster_proteins,
                                 cluster_columns = cluster_samples,
                                 col = c("black", "grey90"),
                                 show_row_dend = show_row_dend,
                                 show_column_dend = show_column_dend,
                                 show_column_names = show_sample_names,
                                 heatmap_legend_param = list(at = c(0,1),labels = c("Missing", "Valid"))
    )
  }
  return(p)
}


#' Plot the intensity distribution of proteins with and without NAs
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_NA_density(tuberculosis_TMT_se)
#'
plot_NA_density <- function(se){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["log2"]])
  dt$ID <- SummarizedExperiment::rowData(se)$Protein.IDs
  melted_dt <- data.table::melt(dt, variable.name = "Column", value.name = "Intensity", id.vars = "ID")
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  melted_dt <- merge(melted_dt, coldata ,by="Column")
  stat <- melted_dt %>% dplyr::group_by(get("ID")) %>% dplyr::summarize(mean = mean(get("Intensity"), na.rm = TRUE), Type = any(is.na(get("Intensity"))))
  stat$Type[stat$Type == TRUE] <- "Missing Value"
  stat$Type[stat$Type == FALSE] <- "Valid Value"
  p <- ggplot2::ggplot(stat, ggplot2::aes(mean, col=get("Type"))) +
    ggplot2::geom_density(na.rm=TRUE) +
    ggplot2::labs(x = "Log2 Intensity", y="Density") +
    ggplot2::scale_colour_manual(name="Value Type", values = c("#A92C23", "#345995"))
  return(p)
}

#' Plot Protein identification overlap (x = Identified in Number of Samples, y=Number of Proteins)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_NA_frequency(tuberculosis_TMT_se)
#'
plot_NA_frequency <- function(se){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])
  dt$ID <- SummarizedExperiment::rowData(se)$Protein.IDs
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  melted_dt <- data.table::melt(dt, variable.name = "Column", value.name = "Intensity", id.vars = "ID")
  melted_dt <- merge(melted_dt, coldata,by="Column")

  nr_samples <- nrow(coldata)
  stat <- melted_dt %>% dplyr::mutate(bin = ifelse(is.na(get("Intensity")),0,1)) %>% dplyr::group_by(get("ID")) %>% dplyr::summarize(sum = sum(get("bin")))
  table <- table(stat$sum) %>% data.table::as.data.table()
  table <- table %>% dplyr::mutate(V1 = factor(table$V1, levels=seq(0,nr_samples)))
  p <- ggplot2::ggplot(table, ggplot2::aes(x=get("V1"), y=get("N"))) +
    ggplot2::geom_col() +
    ggplot2::labs(title="", x = "Identified in Number of Samples", y="Number of Proteins")
  return(p)
}
