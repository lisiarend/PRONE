

#' Function returning some values on the numbers of NA in the data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which data type should be used (default raw)
#'
#' @return list with total amount of values in the data, amount of NA values, and the percentage of NAs
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' get_NA_overview(tuberculosis_TMT_se, ain="log2")
get_NA_overview <- function(se, ain="log2"){
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  na_nr <- sum(is.na(dt))
  tot_nr <- dim(dt)[1] * dim(dt)[2]
  na_perc <- na_nr/tot_nr * 100
  return(data.table::data.table("Total.Values"=c(tot_nr), "NA.Values" = c(na_nr), "NA.Percentage"= c(na_perc)))
}


#' Barplot showing the number of samples per condition
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_condition_overview(tuberculosis_TMT_se, condition = NULL)
plot_condition_overview <- function(se, condition = NULL){
  # get condition
  condition <- get_condition_value(se, condition)

  # prepare data
  md <- data.table::as.data.table(SummarizedExperiment::colData(se))
  dt <- md %>% dplyr::count(get(condition))
  colnames(dt) <- c(condition, "n")

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual",]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(condition), y=get("n"), fill=get(condition))) +
    ggplot2::geom_col() +
    ggplot2::labs(x=condition, y="Number of Samples") +
    ggplot2::scale_fill_manual(name=condition, values=col_vector) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=0.5, vjust=0.5))
  return(p)
}


#' Boxplots of intensities of specific markers
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param markers Vector of the IDs of the markers to plot
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param id_column String specifying the column of the rowData of the SummarizedExperiment object which includes the IDs of the markers
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param shape_by String specifying the column to shape the samples (If NULL or "No", no shaping of samples is done.)
#' @param facet_norm Boolean indicating whether to facet by normalization method (TRUE) or not (FALSE)
#' @param facet_marker Boolean indicating whether to facet by comparison (TRUE) or not (FALSE). Only valid if facet_norm = FALSE.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_markers_boxplots(tuberculosis_TMT_se, markers = c("Q7Z7F0", "Q13790"),
#'                      ain = c("log2"), id_column = "Protein.IDs",
#'                      color_by = NULL,
#'                      shape_by = "Pool",
#'                      facet_norm = FALSE,
#'                      facet_marker = TRUE)
#'
plot_markers_boxplots <- function(se, markers, ain = NULL, id_column = "Protein.IDs", color_by = NULL, shape_by = NULL, facet_norm = TRUE, facet_marker = FALSE){
  # check input parameters
  color_by <- get_color_value(se, color_by)
  shape_by <- get_shape_value(se, shape_by)
  ain <- check_input_assays(se, ain)

  # check id_column
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  if(!id_column %in% colnames(rd)){
    # if id_column not in rowData
    stop(paste0(id_column, " not in rowData of the SummarizedExperiment object!"))
  }

  # prepare data
  overall_dt <- NULL
  for(assay in c(ain)){
    dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[assay]])
    rd <- subset(data.table::as.data.table(SummarizedExperiment::rowData(se)), select = id_column)
    dt <- cbind(rd, dt)
    melted_dt <- data.table::melt(dt, variable.name = "Column", value.name = "Intensity", id.vars = id_column)
    melted_dt$Method <- assay
    if(is.null(overall_dt)){
      overall_dt <- melted_dt
    } else {
      overall_dt <- rbind(overall_dt, melted_dt)
    }
  }

  # merge colData
  cd <- data.table::as.data.table(SummarizedExperiment::colData(se))
  overall_dt <- merge(overall_dt, cd, by = "Column")
  # subset intensities by markers
  found_markers <- stringr::str_detect(overall_dt[[id_column]], paste(markers, collapse = "|") )
  overall_dt <- overall_dt[found_markers,]

  # color vector
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # facet by normalization method
  plots <- list()
  if(facet_norm){
    for(marker in markers){
      found_markers <- stringr::str_detect(overall_dt[[id_column]], marker )
      tmp <- overall_dt[found_markers,]
      p <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(color_by), y = get("Intensity"), fill = get(color_by))) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(x = color_by, fill = color_by, shape = shape_by, y = "Intensity") +
        ggplot2::scale_fill_manual(values = col_vector) +
        ggplot2::facet_wrap(~Method)
      if(is.null(shape_by)){
        p <- p + ggplot2::geom_point(alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(shape = get(shape_by)),alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      }
      plots[[marker]] <- p
    }
  # facet by marker
  } else if(facet_marker){
    for(method in ain){
      tmp <- overall_dt[overall_dt$Method == method,]
      p <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(color_by), y = get("Intensity"), color = get(color_by))) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(x = color_by, color = color_by, shape = shape_by, y = "Intensity") +
        ggplot2::scale_color_manual(values = col_vector) +
        ggplot2::facet_wrap(~Protein.IDs)
      if(is.null(shape_by)){
        p <- p + ggplot2::geom_point(alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(shape = get(shape_by)),alpha =0.8, size = 3, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
      }
      plots[[method]] <- p
    }
  # do not facet by anything
  } else {
    for(method in ain){
      for(marker in markers){
        tmp <- overall_dt[overall_dt$Method == method,]
        found_markers <- stringr::str_detect(tmp[[id_column]], marker )
        tmp <- tmp[found_markers,]
        p <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(color_by), y = get("Intensity"), color = get(color_by))) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(x = color_by, color = color_by, shape = shape_by, y = "Intensity") +
          ggplot2::scale_color_manual(values = col_vector)
        if(is.null(shape_by)){
          p <- p + ggplot2::geom_point(alpha =0.8, size = 1.5, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
        } else {
          p <- p + ggplot2::geom_point(ggplot2::aes(shape = get(shape_by)),alpha =0.8, size = 1.5, position = ggplot2::position_jitterdodge(jitter.width = 0.2, jitter.height = 0))
        }
        plots[[paste0(method, "_", marker)]] <- p
      }
    }
  }
  return(plots)
}

#' Create an UpSet Plot from SummarizedExperiment Data
#'
#' This function generates an UpSet plot from a given SummarizedExperiment object.
#' It allows for the visualization of overlaps between sets defined by a specific
#' column in the metadata. The function supports subsetting to reference samples
#' and customizable color mapping.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used.)
#' @param label_by String specifying the column in the metadata used to label the samples for the UpSet plot
#' @param mb.ratio A numeric vector of length 2, specifying the barplot and matrix area ratios
#' @param only_refs Logical, if TRUE, only reference samples (ComRef) are included in the plot
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_upset(tuberculosis_TMT_se, color_by = NULL, label_by = NULL,
#'            mb.ratio = c(0.7, 0.3), only_refs = FALSE)
#'
plot_upset <- function(se, color_by = NULL, label_by = NULL, mb.ratio = c(0.7,0.3), only_refs = FALSE){

  # get color and label values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  if(!show_sample_names){
    warning("This type of plot requires a labeling. Hence, the sample names are used!")
    label_by <- "Column"
  }

  # Extract metadata as a DataFrame
  md <- data.table::as.data.table(SummarizedExperiment::colData(se))

  # Create color map and assign colors to md$Color
  md$Color <- stats::setNames(RColorBrewer::brewer.pal(n = length(unique(md[[color_by]])), name = "Set3"),
                       unique(md[[color_by]]))[md[[color_by]]]


  # Extract data from the SummarizedExperiment object
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[["raw"]])
  dt[!is.na(dt)] <- 1
  dt[is.na(dt)] <- 0
  # order md according to dt and color_by (safe)
  md <- md[order(md[[color_by]]),]
  dt <- dt[, md$Column, with = FALSE]
  colnames(dt) <- md[[label_by]]

  # If 'only_refs' is TRUE, subset 'md' to include only reference samples
  if (only_refs) {
    refs <- S4Vectors::metadata(se)$refs
    if(is.null(refs)){
      stop("No reference samples in SummarizedExperiment!")
    } else {
      md <- md[md$Column %in% refs,]
    }
  }

  p <- UpSetR::upset(dt, sets = md[[label_by]], mb.ratio = mb.ratio, keep.order = TRUE,
                     order.by = "degree", sets.bar.color = md$Color,
                     main.bar.color = "#a4bccc", matrix.color = "black",
                     shade.color = "wheat4")
  return(p)
}

#' Plot a heatmap of the sample intensities with optional column annotations for a selection of normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param color_by Vector of strings specifying the columns to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bars added.)
#' @param label_by String specifying the column in the metadata used to label the samples for the UpSet plot
#' @param only_refs Logical, if TRUE, only reference samples (ComRef) are included in the plot
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_heatmap(tuberculosis_TMT_se, ain = c("log2"), color_by = NULL,
#'              label_by = NULL, only_refs = FALSE)
#'
plot_heatmap <- function(se, ain = NULL, color_by = c("Group", "Pool"), label_by = NULL, only_refs = FALSE){
  plots <- list()

    # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }
  # get color_by and label_by values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  # If 'only_refs' is TRUE, subset 'md' to include only reference samples
  if (only_refs) {
    refs <- S4Vectors::metadata(se)$refs
    if(is.null(refs)){
      stop("No reference samples in SummarizedExperiment!")
    } else {
      coldata <- coldata[coldata$Column %in% refs, ]
    }
  }

  # build heatmap annotation
  df_anno <- subset(coldata, select = c(color_by))
  # color vector
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  if(!is.null(color_by)){
    # color generation for annotation
    i <- 1
    color_list <- list()
    for(col in c(color_by)){
      vals <- unique(df_anno[[col]])
      vals <- vals[!is.na(vals)]
      colors <- col_vector[i:(i+length(vals)-1)]
      names(colors) <- vals
      i <- i+length(vals)
      color_list[[col]] <- colors
    }

    # heatmap annotations
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      df = df_anno,
      col = color_list,
      gp = grid::gpar(col = "black")
    )
  } else {
    top_anno <- NULL
  }



  for(method in ain){
    data <- data.table::as.data.table(SummarizedExperiment::assays(se)[[method]])
    if(only_refs){
      data <- data[,c(refs)]
    }
    # order clustering
    clustering <- dendsort::dendsort(stats::hclust(stats::dist(t(data))))
    p <- ComplexHeatmap::Heatmap(stats::na.omit(as.matrix(data)),
                                 name = "Intensity",
                                 column_labels = coldata[[label_by]],
                                 show_column_names = show_sample_names,
                                 top_annotation = top_anno,
                                 show_row_names = FALSE,
                                 cluster_columns = clustering,
                                 use_raster = FALSE)
    plots[[ain]] <- p
  }
  return(plots)
}

#' Plot correlation, histogram, and scatterplot of samples of a specific group for a selection of normalization methods
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition Column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param label_by String specifying the column in the metadata used to label the samples for the UpSet plot
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_pairs_panels(tuberculosis_TMT_se, ain = c("log2", "IRS_on_RobNorm"),
#'                   condition = NULL, label_by = NULL)
#'
plot_pairs_panels <- function(se, ain = NULL, condition = NULL, label_by = NULL){
  plots <- list()
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }

  # get label values
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  if(!show_sample_names){
    warning("This type of plot requires a labeling. Hence, the sample names are used!")
    label_by <- "Column"
  }

  # get condition
  condition <- get_condition_value(se, condition)

  # condition vector
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  condition_vector <- coldata[[condition]]

  # get group
  group <- unique(condition_vector)

  # function for scatter plots with smoothed trend line
  upper_plots <- function(data, mapping, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_point(color = "black", size = 1, alpha = 1) +
      ggplot2::geom_smooth(method = "lm", color = "red", size = 1)
  }

  # working with a trick of global assignment
  diag_plots <- function(data, mapping, ...) {
    # increase counter each run globally so outside the function as well and this does the trick!
    ggplot2::ggplot(data = data, mapping = mapping) +
      # choose color by counter and send bin width argument in
      ggplot2::geom_histogram(fill = "#79CAFF") +
      ggplot2::geom_density()
  }

  # loop over data
  for(method in c(ain)){
    plots[[method]] <- list()
    data <- data.table::as.data.table(SummarizedExperiment::assays(se)[[method]])
    coldata <- coldata[match(colnames(data), coldata$Column),]
    colnames(data) <- coldata[[label_by]]
    for(g in c(group)){
      dt <- data[,coldata[condition_vector == g,][[label_by]], with = FALSE]
      p <- GGally::ggpairs(dt,
              title = paste0(g, " (", method, " data)"),
              lower = list(continuous = "cor"),
              diag = list(continuous = GGally::wrap(diag_plots, color = "red", se=FALSE)),
              upper = list(continuous = GGally::wrap(upper_plots), bins = 10))
      plots[[method]][[g]] <- p
    }
  }
  return(plots)
}


