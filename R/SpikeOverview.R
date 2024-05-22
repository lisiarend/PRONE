
#' Plot number of identified spike-in proteins per sample.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)#'
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' plot_identified_spiked_proteins(spike_in_se, color_by = NULL,
#'                                 label_by = NULL)
#'
plot_identified_spiked_proteins <- function(se, color_by = NULL, label_by = NULL){
  # check parameters
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # extract data
  dt <- SummarizedExperiment::assays(se)$log2
  spike_column <- S4Vectors::metadata(se)$spike_column
  spike_val <- S4Vectors::metadata(se)$spike_value
  rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))

  # Make Data Binary (NA --> 0, all other values --> 1)
  dt[is.na(dt)] <- 0
  dt[dt != 0] <- 1
  dt$Spiked <- rowdata[[spike_column]]

  # How Many spiked and BG Proteins per Sample
  stats <- dt %>% dplyr::group_by(Spiked) %>% dplyr::summarise_each(dplyr::funs(sum))
  stats <- t(stats)
  colnames(stats) <- stats[1,]
  stats <- stats[2:nrow(stats),]
  column <- row.names(stats)
  stats <- data.table::as.data.table(stats)
  stats$Column <- column

  # Add coldata
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  stats <- merge(stats, coldata, by="Column")
  stats$Spiked <- as.numeric(stats[,spike_val])

  # Ordering
  if(!is.null(color_by)){
    stats <- stats[order(stats[[color_by]]),]

    # color vector
    qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- rev(col_vector)

    # Plot
    p <- ggplot2::ggplot(stats, ggplot2::aes(x=get("Label"), y=get("Spiked"), group = get(color_by), fill=get(color_by), label=get("Spiked"))) +
      ggplot2::geom_bar(stat = "identity", width=0.8, position=ggplot2::position_dodge(width=0.2)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::labs(x = "Samples", y=paste0("Number of Identified ", spike_val," Proteins")) +
      ggplot2::scale_fill_manual(name=color_by, values=col_vector) +
      ggplot2::geom_label(fill="white", show.legend = FALSE, color="black")

  } else {
    # Plot
    p <- ggplot2::ggplot(stats, ggplot2::aes(x=get("Label"), y=get("Spiked"), label=get("Spiked"))) +
      ggplot2::geom_bar(stat = "identity", width=0.8, position=ggplot2::position_dodge(width=0.2)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::labs(x = "Samples", y=paste0("Number of Identified ", spike_val," Proteins")) +
      ggplot2::geom_label(fill="white", show.legend = FALSE, color="black")
  }
  if(!show_sample_names){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  return(p)
}

#' Plot histogram of the spike-in and background protein intensities per condition.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' plot_histogram_spiked(spike_in_se, condition = NULL)
#'
plot_histogram_spiked <- function(se, condition = NULL){
  # get condition
  condition <- get_condition_value(se, condition)

  dt <- SummarizedExperiment::assays(se)$log2
  spike_column <- S4Vectors::metadata(se)$spike_column
  spike_val <- S4Vectors::metadata(se)$spike_value

  rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))

  dt[[spike_column]] <- rowdata[[spike_column]]

  dt <- data.table::melt(dt, id.vars = spike_column, variable.name = "Column", value.name = "Intensity")
  dt <- merge(dt, coldata, by = "Column")

  types <- unique(dt[[spike_column]])
  levels <- c(spike_val, types[!types %in% spike_val])
  dt[[spike_column]] <- factor(dt[[spike_column]], levels = levels)


  colors <- c("red", "grey80")

  p <- ggplot2::ggplot(dt, ggplot2::aes(x=get("Intensity"), fill=get(spike_column))) +
    ggplot2::geom_histogram(bins = 50, position="identity", alpha=0.7) +
    ggplot2::facet_wrap(~get(condition)) +
    ggplot2::scale_fill_manual(name = spike_column, values = colors) +
    ggplot2::labs(x="log2 Intensity", y="Protein Number")
  return(p)
}

#' Plot profiles of the spike-in and background proteins using the log2 average protein intensities as a function of the different concentrations.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param xlab String for the x-label of the plot
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' plot_profiles_spiked(spike_in_se, xlab = "Concentration")
#'
plot_profiles_spiked <- function(se, xlab = "Concentration"){
  # prepare data
  dt <- SummarizedExperiment::assays(se)[["log2"]]
  spike_column <- S4Vectors::metadata(se)$spike_column
  spike_val <- S4Vectors::metadata(se)$spike_value
  spike_concentration <- S4Vectors::metadata(se)$spike_concentration

  rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  samples <- colnames(dt)

  dt$Protein.IDs <- rowdata[,"Protein.IDs"]

  dt <- data.table::melt(data.table::as.data.table(dt), id.vars = "Protein.IDs", measure.vars = colnames(dt), variable.name = "Column", value.name = "Intensity")

  dt <- merge(dt, coldata, by = "Column")
  dt <- merge(dt, rowdata[, c("Protein.IDs", spike_column), with = FALSE], by = "Protein.IDs")
  bg_value <- unique(dt[[spike_column]])[unique(dt[[spike_column]]) != spike_val]


  col_vector <- c("red", "grey70")
  dt[[spike_column]] <- factor(dt[[spike_column]], levels = c(spike_val, bg_value))

  group <- c("Protein.IDs", spike_concentration, spike_column)
  dt$Intensity <- as.numeric(dt$Intensity)
  dt[[spike_concentration]] <- as.numeric(dt[[spike_concentration]])

  average <- dt %>% dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>% dplyr::summarize(Mean = mean(Intensity, na.rm=TRUE)) %>% data.table::as.data.table()

  p <- ggplot2::ggplot(average, ggplot2::aes(x= get(spike_concentration), y = get("Mean"), group = get("Protein.IDs"), color = get(spike_column))) +
  ggplot2::geom_point() +
  ggplot2::geom_line(alpha = 0.7) +
  ggplot2::labs(x = xlab, y = "Log2 Average Intensity") +
  ggplot2::scale_x_continuous(breaks = unique(average[, spike_concentration])) +
  ggplot2::theme_bw() +
  ggplot2::scale_color_manual(name = spike_column, values = col_vector)  +
  ggplot2::facet_wrap(~get(spike_column), scales = "free_y")
  return(p)
}
