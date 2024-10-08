################# ----------- Filter Samples Functions ----------- ###################

#' Remove outliers samples detected by the detect_outliers_POMA function
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param poma_res_outliers Outliers data.table returned by the detect_outliers_POMA function
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' poma_res <- detect_outliers_POMA(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- remove_POMA_outliers(tuberculosis_TMT_se, poma_res$outliers)
#'
remove_POMA_outliers <- function(se, poma_res_outliers){
  if(nrow(poma_res_outliers) == 0){
    message("No outlier samples detected! Hence, no samples were removed.")
    return(se)
  } else {
    se_subset <- se[, !se$Column %in% poma_res_outliers$sample]
    rm <- ncol(se) - ncol(se_subset)
    if(rm == 0){
      message("No outlier samples detected! Hence, no samples removed.")
      return(se)
    } else {
      message(paste0(rm, " outlier samples removed."))
      return(se_subset)
    }
  }
}


#' Remove samples with specific value in column manually
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param column String specifying the column of the meta data (samples with the specified value in this column will be removed)
#' @param values Vector of Strings specifying the value for the removal of samples (samples with this value in the specified column will be removed)
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- remove_samples_manually(tuberculosis_TMT_se,
#'                                 column = "Label", values = c("1.HC_Pool1"))
#'
remove_samples_manually <- function(se, column, values){
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  # check if column in object
  if(! column %in% colnames(coldata)){
    stop(paste0(column, " not in data!"))
  } else {
    col <- coldata[[column]]
    # check if there is a sample with the specific value in the corresponding column
    if(sum(!c(values) %in% col) == length(c(values))){
      stop(paste0("No sample(s) found with ", values, " in ", column, "column!"))
    } else {
      se_subset <- se[, !se[[column]] %in% values]
      # check if all samples would be removed
      if(ncol(se_subset) == 0){
        stop("All samples would be removed. Aborted!")
      } else {
        rm <- ncol(se) - ncol(se_subset)
        message(paste0(rm, " samples removed."))
        return(se_subset)
      }
    }
  }
}


#' Remove reference samples of SummarizedExperiment object (reference samples specified during loading)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#'
#' @return filtered SummarizedExperiment object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' tuberculosis_TMT_se <- remove_reference_samples(tuberculosis_TMT_se)
#'
remove_reference_samples <- function(se){
  refs <- S4Vectors::metadata(se)$refs
  condition <- S4Vectors::metadata(se)$condition
  if(is.null(refs)){
    stop("No reference samples specified during data loading!")
  } else {
    se_subset <- se[,! se$Column %in% c(refs)]
    rm <- ncol(se) - ncol(se_subset)
    message(paste0(rm, " reference samples removed from the SummarizedExperiment object."))
    if(!is.null(levels(se[[condition]]))){
      se_subset[[condition]] <- droplevels(se_subset[[condition]])
    }
    return(se_subset)
  }
}


################# ----------- Plotting Functions ----------- ###################

#' Plot number of non-zero proteins per sample
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which data type should be used (default raw)
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_nr_prot_samples(tuberculosis_TMT_se, ain="raw", color_by = "Group",
#'                      label_by = "Label")
#'
plot_nr_prot_samples <- function(se, ain="raw", color_by = NULL, label_by = NULL){
  # get color and label values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # prepare data
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt$ID <- rownames(dt)
  dt <- data.table::melt(dt, id.vars = "ID", variable.name = "Column", value.name = "Intensity")
  dt$Intensity[is.na(dt$Intensity)] <- 0
  dt <- dt[dt$Intensity != 0,]
  dt <- dt %>% dplyr::count(Column)
  dt <- merge(dt, data.table::as.data.table(SummarizedExperiment::colData(se)), by="Column")
  dt <- data.table::as.data.table(dt)

  # order data
  if(!is.null(color_by)){
    dt[, color_by] <- factor(dt[,color_by], levels = unique(dt[,color_by]))
    dt <- dt[order(dt[,color_by]),]
    if(show_sample_names){
      dt[, label_by] <- factor(dt[,label_by], levels = unique(dt[,label_by]))
    } else {
      dt$Column <- factor(dt$Column, levels = dt$Column)
    }
  }

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  if(is.null(color_by)){
    p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(label_by), y=get("n"))) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      ggplot2::labs(x = "Samples", y ="Number of Proteins \n (non-zero protein abundance)")
  } else {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(label_by), y=get("n"), fill=get(color_by))) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      ggplot2::scale_fill_manual(name = color_by, values = col_vector ) +
      ggplot2::labs(x = "Samples", y ="Number of Proteins \n (non-zero protein abundance)")
  }
  if(!show_sample_names){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  }
  return(p)
}

#' Plot total protein intensity per sample
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which data type should be used (default raw)
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#'
#' @return list of a ggplot object and the dataframe of outliers
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_tot_int_samples(tuberculosis_TMT_se, ain="raw", color_by = NULL,
#'                      label_by = NULL)
#'
plot_tot_int_samples <- function(se, ain="raw", color_by = NULL, label_by = NULL){
  # get color and label values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # prepare data
  dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
  dt$ID <- rownames(dt)
  dt <- data.table::melt(dt, id.vars = "ID", variable.name = "Column", value.name = "Intensity")
  dt$Intensity[is.na(dt$Intensity)] <- 0
  dt <- dt %>% dplyr::group_by(Column) %>% dplyr::summarise(total = sum(Intensity)) %>% data.table::as.data.table()
  dt <- merge(dt, data.table::as.data.table(SummarizedExperiment::colData(se)), by="Column")
  dt$total <- as.numeric(dt$total)

  # order data
  if(!is.null(color_by)){
    dt[, color_by] <- factor(dt[,color_by], levels = unique(dt[,color_by]))
    dt <- dt[order(dt[,color_by]),]
    if(show_sample_names){
      dt[, label_by] <- factor(dt[,label_by], levels = unique(dt[,label_by]))
    } else {
      dt$Column <- factor(dt$Column, levels = dt$Column)
    }
  }

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  if(is.null(color_by)){
    p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(label_by), y=get("total"))) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      ggplot2::labs(x = "Samples", y ="Total Protein Intensity")
  } else {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x=get(label_by), y=get("total"), fill=get(color_by))) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      ggplot2::scale_fill_manual(name = color_by, values = col_vector ) +
      ggplot2::labs(x = "Samples", y ="Total Protein Intensity")
  }
  if(!show_sample_names){
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  }
  return(p)
}


#' Outlier detection via POMA R Package
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which data type should be used (default raw)
#' @param condition Column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param method String specifying the method that should be used to calculate the distance matrix
#' @param type String specifying the type of distance calculation to centroid or spatial median
#' @param group String specifying if the outlier detection should be performed multi-variate
#'             (with conditions) or on the complete data set
#' @param coeff This value corresponds to the classical 1.5 in Q3 + 1.5 * IQR formula to detect outliers. By changing this value, the permissiveness in outlier detection will change.
#' @importFrom magrittr %>%
#' @importFrom ggtext element_textbox_simple
#'
#' @return list of two ggplot objects and a data.table with outlier samples
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' poma_res <- detect_outliers_POMA(tuberculosis_TMT_se, ain="raw",
#'                                  condition = NULL, method="euclidean",
#'                                  type="median", group=TRUE, coeff = 1.5)
detect_outliers_POMA <- function(se, ain="log2", condition = NULL, method="euclidean", type="median", group=TRUE, coeff = 1.5){

  # get condition
  condition <- get_condition_value(se, condition)
  
  stopifnot(length(ain) == 1)

  # prepare data
  md <- data.table::as.data.table(SummarizedExperiment::colData(se))
  if(group){
    condition <- S4Vectors::metadata(se)$condition
    condition_vector <- md[[condition]]
  } else {
    condition <- "complete"
    condition_vector <- rep("complete",nrow(md))
  }
  new_md <- data.table::data.table("Condition" = as.factor(condition_vector), "Column" = md$Column)
  colnames(new_md)[1] <- condition
  data <- SummarizedExperiment::SummarizedExperiment(assays = SummarizedExperiment::assays(se)[[ain]], colData = new_md, rowData=SummarizedExperiment::rowData(se))
  
  # perform POMA outlier detection
  browser()
  poma_res <- POMA::PomaOutliers(data, method=method, type=type, coef = coeff)

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # modify plots
  poma_res$polygon_plot <- poma_res$polygon_plot +
    ggplot2::scale_fill_manual(name = condition, values = col_vector) +
    ggplot2::scale_color_manual(condition, values = c(rep("black", length(unique(md[[condition]]))))) +
    ggplot2::scale_shape_manual(name = condition, values = seq_len(length(unique(md[[condition]]))))

  poma_res$distance_boxplot <- poma_res$distance_boxplot +
    ggplot2::scale_fill_manual(name = condition, values = col_vector) +
    ggplot2::labs(x = condition, y = "Distance to Group Centroid", fill = condition)

  poma_res$outliers <- data.table::as.data.table(poma_res$outliers)

  return(poma_res)
}
