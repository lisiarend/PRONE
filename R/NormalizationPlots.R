## ----- Helper Functions ----- ##

#' Function to get a long data table of all intensities of all kind of normalization
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input (default NULL)
#'            If NULL then all normalization of the SummarizedExperiment object are plotted next to each other (except raw).
#'
#' @return data table
#'
get_complete_dt <- function(se, ain=NULL){
  complete_dt <- NULL
  if(is.null(ain)){
    assays <- names(SummarizedExperiment::assays(se))
  } else {
    assays <- ain
  }
  for(ain in assays){
    if (ain != "raw"){
      dt <- data.table::as.data.table(SummarizedExperiment::assays(se)[[ain]])
      dt$ID <- data.table::as.data.table(SummarizedExperiment::rowData(se))$Protein.IDs
      melted_dt <- data.table::melt(dt, variable.name = "Column", value.name = "Intensity", id.vars = "ID")
      melted_dt <- merge(melted_dt, data.table::as.data.table(SummarizedExperiment::colData(se)),by="Column")
      melted_dt$Assay <- rep(ain,nrow(melted_dt))
      if (is.null(complete_dt)){
        complete_dt <- melted_dt
      } else {
        complete_dt <- rbind(complete_dt, melted_dt)
      }
    }
  }
  return(complete_dt)
}

#' Function to get a long data table of all PCA1 and PCA2 values of all kind of normalization
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain String which assay should be used as input (default NULL)
#'            If NULL then all normalization of the SummarizedExperiment object are plotted next to each other (except raw).
#'
#' @return data table
#'
get_complete_pca_dt <- function(se, ain=NULL){
  complete_pca_dt <- NULL
  if(is.null(ain)){
    assays <- names(SummarizedExperiment::assays(se))
  } else {
    assays <- c(ain)
  }
  for(ain in assays){
    if (ain != "raw"){
      dt <- SummarizedExperiment::assays(se)[[ain]]
      dt <- dt[!apply(is.na(dt), 1, any),]
      # check if protein group is constant
      rm <- which(apply(t(dt), 2, stats::var)== 0)
      if (!identical(rm, integer(0))){
        dt <- dt[-rm,]
      }
      pca <- stats::prcomp(t(dt), scale.=TRUE)
      pc1 <- round(summary(pca)$importance["Proportion of Variance","PC1"] * 100, 2)
      pc2 <- round(summary(pca)$importance["Proportion of Variance","PC2"] * 100, 2)
      pca_dt <- data.table::data.table(PC1 = pca$x[, "PC1"],
                                       PC2 = pca$x[, "PC2"],
                                       Column = row.names(pca$x))
      pca_dt$Assay <- ain
      pca_dt$PC1_Perc <- pc1
      pca_dt$PC2_Perc <- pc2
      if (is.null(complete_pca_dt)){
        complete_pca_dt <- pca_dt
      } else {
        complete_pca_dt <- rbind(complete_pca_dt, pca_dt)
      }
    }
  }
  return(complete_pca_dt)
}

## ----- Qualitative Evaluation Functions ----- ##

#' Plot the distributions of the normalized data as boxplots
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#' @param facet_norm Boolean specifying whether to facet by normalization methods (default TRUE). If FALSE, list of plots of the different normalized data is returned.
#' @param ncol Number of columns in plot (for faceting)
#'
#'
#' @return if facet_norm = TRUE, ggplot object, else list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_boxplots(tuberculosis_TMT_se, ain = NULL, color_by = NULL, label_by = NULL,
#'               facet_norm = TRUE, ncol = 3)
#' plot_boxplots(tuberculosis_TMT_se, ain = c("log2", "IRS_on_RobNorm"), color_by = "Pool",
#'               label_by = "Label", facet_norm = FALSE)
#'
plot_boxplots <- function(se, ain = NULL, color_by = NULL, label_by = NULL, facet_norm = TRUE, ncol = 3){
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }

  # get color and label values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # get data for plot
  melted_dt <- get_complete_dt(se, ain=ain)
  melted_dt <- melted_dt[order(melted_dt[,color_by]),]

  if(show_sample_names){
    melted_dt$Label <- factor(melted_dt[,label_by], levels = unique(melted_dt[,label_by]))
  } else {
    melted_dt$Label <- factor(melted_dt$Column, levels = unique(melted_dt$Column))
  }

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  if(facet_norm){
    p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x=get("Intensity"), y=get(label_by), fill=get(color_by))) +
      ggplot2::geom_boxplot() +
      ggplot2::stat_boxplot(geom ='errorbar', width = 0.4) +
      ggplot2::scale_fill_manual(name = color_by, values = col_vector) +
      ggplot2::labs(x="Intensity", y="Samples") +
      ggplot2::facet_wrap(~Assay, scales="free_x", ncol=ncol)
    if(!show_sample_names){
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }
  } else {
    p <- list()
    for(method in c(ain)){
      dt <- melted_dt[melted_dt$Assay == method,]
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("Intensity"), y=get(label_by), fill=get(color_by))) +
        ggplot2::geom_boxplot() +
        ggplot2::stat_boxplot(geom ='errorbar', width = 0.4) +
        ggplot2::scale_fill_manual(name = color_by, values = col_vector) +
        ggplot2::labs(x="Intensity", y="Samples")
      if(!show_sample_names){
        tmp <- tmp + ggplot2::theme(axis.text.y = ggplot2::element_blank())
      }
      p[[method]] <- tmp
    }
  }
  return(p)
}

#' Plot the densities of the normalized data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param facet_norm Boolean specifying whether to facet by normalization methods (default TRUE). If FALSE, list of plots of the different normalized data is returned.
#' @param ncol Number of columns in plot (for faceting)
#'
#'
#' @return if facet_norm = TRUE, ggplot object, else list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_densities(tuberculosis_TMT_se, ain = NULL, color_by = NULL,
#'                           facet_norm = TRUE, ncol = 3)
#' plot_densities(tuberculosis_TMT_se, ain = c("log2", "IRS_on_RobNorm"),
#'                           color_by = "Label",
#'               facet_norm = FALSE)
#'
plot_densities <- function(se, ain = NULL, color_by = NULL, facet_norm = TRUE, ncol = 3){
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }

  # get color and label values
  color_by <- get_color_value(se, color_by)

  # get data for plot
  melted_dt <- get_complete_dt(se, ain=ain)

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  if(facet_norm){
    p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x=get("Intensity"), color=get(color_by))) +
      ggplot2::geom_density(na.rm=TRUE) +
      ggplot2::labs(x="Intensity", y="Density") +
      ggplot2::facet_wrap(~Assay, scales="free", ncol = ncol) +
      ggplot2::scale_color_manual(name = color_by, values = col_vector)
  } else {
    p <- list()
    for(method in c(ain)){
      dt <- melted_dt[melted_dt$Assay == method,]
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("Intensity"), color=get(color_by))) +
        ggplot2::geom_density(na.rm=TRUE) +
        ggplot2::labs(x="Intensity", y="Density") +
        ggplot2::scale_color_manual(name = color_by, values = col_vector)
      p[[method]] <- tmp
    }
  }
  return(p)
}


#' PCA plot of the normalized data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#' @param shape_by String specifying the column to shape the samples (If NULL or "No", no shaping of samples is done.)
#' @param facet_norm Boolean specifying whether to facet by normalization methods (default TRUE). If FALSE, list of plots of the different normalized data is returned. However, then you can also facet by any column of the metadata.
#' @param facet_by String specifying the column to facet the samples (If facet = FALSE, the plot will not be faceted by the normalization methods, but instead a list of plots of each normalization method is returned. Then, the PCA plot
#'                 can be faceted by any column of the metadata, for instance by "Batch". If facet_by is NULL or "No", no faceting is performed.)
#' @param ellipse Boolean to indicate if ellipses should be drawn
#' @param ncol Number of columns in plot (for faceting)
#'
#' @return if facet_norm = TRUE, ggplot object, else list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_PCA(tuberculosis_TMT_se, ain = NULL, color_by = NULL, label_by = NULL,
#'          shape_by = "Pool",
#'          facet_norm = TRUE, ncol = 3)
#' plot_PCA(tuberculosis_TMT_se, ain = c("IRS_on_RobNorm"), color_by = "Group",
#'               label_by = "Label", facet_norm = FALSE, facet_by = "Pool")
#'
plot_PCA <- function(se, ain=NULL, color_by = NULL, label_by = NULL, shape_by = NULL, facet_norm = TRUE, facet_by = NULL, ellipse = FALSE, ncol = 3){
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }

  # get color and label values
  color_by <- get_color_value(se, color_by)
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]

  # get shape value
  shape_by <- get_shape_value(se, shape_by)

  # get facet by value
  if(!facet_norm){
    facet_by <- get_facet_value(se, facet_by)
  }

  # get data for plot
  pca_dt <- get_complete_pca_dt(se, ain=ain)
  pca_dt <- merge(pca_dt, data.table::as.data.table(SummarizedExperiment::colData(se)), by="Column")

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  if(facet_norm){
    # facet by normalization method
    pca_dt$AssayTitle <- paste0(pca_dt$Assay, " (PC1: ", pca_dt$PC1_Perc, "% & PC2: ", pca_dt$PC2_Perc, "%)")
    if(!show_sample_names){
      # no label
      if(is.null(color_by)){
        # no coloring
        if(is.null(shape_by)){
          # no shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol)
        } else {
          # shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), shape = get(shape_by))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol) +
            ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
        }
      } else {
        # coloring
        if(is.null(shape_by)){
          # no shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_color_manual(name = color_by, values = col_vector) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol)
        } else {
          # shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by), shape = get(shape_by))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_color_manual(name = color_by, values = col_vector) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol) +
            ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
        }
      }
    } else {
      # label
      if(is.null(color_by)){
        # no coloring
        if(is.null(shape_by)){
          # no shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), label = get(label_by))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::geom_text(hjust = 0.5, vjust = 1) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol)
        } else {
          # shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), label = get(label_by)), shape = get(shape_by)) +
            ggplot2::geom_point(size = 3) +
            ggplot2::geom_text(hjust = 0.5, vjust = 1) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol) +
            ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
        }
      } else {
        # coloring
        if(is.null(shape_by)){
          # no shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by), label = get(label_by))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_color_manual(name = color_by, values = col_vector) +
            ggplot2::geom_text(hjust = 0.5, vjust = 1) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol)
        } else {
          # shaping
          p <- ggplot2::ggplot(pca_dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by), label = get(label_by), shape = get(shape_by))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_color_manual(name = color_by, values = col_vector) +
            ggplot2::geom_text(hjust = 0.5, vjust = 1) +
            ggplot2::xlab("PC1") +
            ggplot2::ylab("PC2") +
            ggplot2::facet_wrap(~AssayTitle, scales = "free", ncol = ncol) +
            ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
        }
      }
    }

    if(ellipse){
      p <- p + ggplot2::stat_ellipse()
    }
  } else {
    # no facet by normalization method --> individual plots for each normalization method --> but facet_by maybe any other column
    p <- list()
    for(method in c(ain)){
      dt <- pca_dt[pca_dt$Assay == method,]
      if(!show_sample_names){
        # no label
        if(is.null(color_by)){
          # no coloring
          if(is.null(shape_by)){
            # no shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)"))
          } else {
            # shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), shape = get(shape_by))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)")) +
              ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
          }
        } else {
          # coloring
          if(is.null(shape_by)){
            # no shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::scale_color_manual(name = color_by, values = col_vector) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)"))
          } else {
            # shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by), shape = get(shape_by))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::scale_color_manual(name = color_by, values = col_vector) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)")) +
              ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
          }
        }
      } else {
        # label
        if(is.null(color_by)){
          # no coloring
          if(is.null(shape_by)){
            # no shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), label = get(label_by))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::geom_text(hjust = 0.5, vjust = 1) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)"))
          } else {
            # shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), label = get(label_by)), shape = get(shape_by)) +
              ggplot2::geom_point(size = 3) +
              ggplot2::geom_text(hjust = 0.5, vjust = 1) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)")) +
              ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
          }
        } else {
          # coloring
          if(is.null(shape_by)){
            # no shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by), label = get(label_by))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::scale_color_manual(name = color_by, values = col_vector) +
              ggplot2::geom_text(hjust = 0.5, vjust = 1) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)"))
          } else {
            # shaping
            tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("PC1"), y = get("PC2"), color = get(color_by), label = get(label_by), shape = get(shape_by))) +
              ggplot2::geom_point(size = 3) +
              ggplot2::scale_color_manual(name = color_by, values = col_vector) +
              ggplot2::geom_text(hjust = 0.5, vjust = 1) +
              ggplot2::xlab(paste0("PC1 (", unique(dt$PC1_Perc), " %)")) +
              ggplot2::ylab(paste0("PC2 (", unique(dt$PC2_Perc), " %)")) +
              ggplot2::guides(shape=ggplot2::guide_legend(shape_by))
          }
        }
      }
      if(!is.null(facet_by)){
        tmp <- tmp + ggplot2::facet_wrap(~get(facet_by), ncol = ncol)
      }
      if(ellipse){
        tmp <- tmp + ggplot2::stat_ellipse()
      }
      p[[method]] <- tmp
    }
  }
  return(p)
}

## ----- Quantitative Evaluation Functions ----- ##

#' Plot intragroup correlation of the normalized data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param method String specifying the method for correlation calculation (pearson, spearman or kendall)
#'
#' @return ggplot object (boxplot)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_intragroup_correlation(tuberculosis_TMT_se, ain = NULL,
#'                             condition = NULL, method = "pearson")
#'
plot_intragroup_correlation <- function(se, ain = NULL, condition = NULL, method = "pearson"){
  # check method
  stopifnot(method %in% c("pearson", "spearman", "kendall"))

  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }
  assays <- SummarizedExperiment::assays(se)
  exclude <- names(assays)[!names(assays) %in% ain]
  for(e in exclude){
    assays[e] <- NULL
  }

  # get condition
  condition <- get_condition_value(se, condition)

  # condition vector
  condition_vector <- data.table::as.data.table(SummarizedExperiment::colData(se)[, condition])

  # prepare data
  sampleGroupsWithReplicates <- names(table(condition_vector)[table(condition_vector) > 1])
  repCor <- NormalyzerDE:::calculateSummarizedCorrelationVector(
    assays,
    condition_vector,
    sampleGroupsWithReplicates,
    method
  )
  cor_intra <- data.table::as.data.table(repCor)
  cor_intra <- data.table::melt(cor_intra, measure.vars = colnames(cor_intra), variable.name = "Normalization", value.name = "Correlation")
  cor_intra$Normalization <- factor(cor_intra$Normalization, levels = sort(as.character(unique(cor_intra$Normalization))))

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  p <- ggplot2::ggplot(cor_intra, ggplot2::aes(x=get("Normalization"), y=get("Correlation"), fill=get("Normalization"))) +
    ggplot2::geom_boxplot() +
    ggplot2::stat_boxplot(geom ='errorbar', width = 0.4) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
    ggplot2::ylab("Pearson Correlation") +
    ggplot2::xlab("Normalization Method") +
    ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector)
  return(p)
}

#' Plot intragroup pooled coefficient of variation (PCV) of the normalized data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition Column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param diff Boolean indicating whether to visualize the reduction of intragroup variation (PCV) compared to "log" (TRUE) or a normal boxplot of intragroup variation (PCV) for each normalization method (FALSE).
#'
#' @return ggplot object (boxplot)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_intragroup_PCV(tuberculosis_TMT_se, ain = NULL,
#'                             condition = NULL, diff = FALSE)
#'
plot_intragroup_PCV <- function(se, ain = NULL, condition = NULL, diff = FALSE){
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }
  assays <- SummarizedExperiment::assays(se)
  exclude <- names(assays)[!names(assays) %in% ain]
  for(e in exclude){
    assays[e] <- NULL
  }

  assays["raw"] <- NULL

  # get condition
  condition <- get_condition_value(se, condition)

  # prepare data
  avgCVPerNormAndReplicates <- NormalyzerDE:::calculateReplicateCV(assays, data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]])
  PCV_intra <- data.table::as.data.table(avgCVPerNormAndReplicates)
  PCV_intra <- data.table::melt(PCV_intra, measure.vars = colnames(PCV_intra), variable.name = "Normalization", value.name = "PCV")
  PCV_intra$Normalization <- factor(PCV_intra$Normalization, levels = sort(as.character(unique(PCV_intra$Normalization))))

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  if(diff){
    # check if log in assays
    if(!"log2" %in% names(assays)){
      stop("Log2 data not in SummarizedExperiment! Difference not applicable.")
    }
    avgcvmempdiff <- NormalyzerDE:::calculatePercentageAvgDiffInMat(avgCVPerNormAndReplicates)
    PCV_diff <- data.table::data.table(Normalization = names(assays), PCV = avgcvmempdiff)
    PCV_diff$Normalization <- factor(PCV_diff$Normalization, levels=sort(as.character(unique(PCV_diff$Normalization))))
    p <- ggplot2::ggplot(PCV_diff, ggplot2::aes(x=get("Normalization"), y=get("PCV"),fill=get("Normalization"))) +
      ggplot2::geom_col() +
      ggplot2::ylab("") +
      ggplot2::xlab("Normalization Method") +
      ggplot2::geom_label(ggplot2::aes(label=round(get("PCV"), digits = 0)), fill="white", show.legend = FALSE, color="black") +
      ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  } else {
    p <- ggplot2::ggplot(PCV_intra, ggplot2::aes(x=get("Normalization"), y=get("PCV"), fill=get("Normalization"))) +
      ggplot2::geom_boxplot() +
      ggplot2::stat_boxplot(geom ='errorbar', width = 0.4) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
      ggplot2::ylab("PCV") +
      ggplot2::xlab("Normalization Method") +
      ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector)
  }
  return(p)
}

#' Plot intragroup pooled median absolute deviation (PMAD) of the normalized data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param diff Boolean indicating whether to visualize the reduction of intragroup variation (PMAD) compared to "log" (TRUE) or a normal boxplot of intragroup variation (PMAD) for each normalization method (FALSE).
#'
#' @return ggplot object (boxplot)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_intragroup_PMAD(tuberculosis_TMT_se, ain = NULL,
#'                             condition = NULL, diff = FALSE)
#'
plot_intragroup_PMAD <- function(se, ain=NULL, condition = NULL, diff = FALSE){
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }
  assays <- SummarizedExperiment::assays(se)
  exclude <- names(assays)[!names(assays) %in% ain]
  for(e in exclude){
    assays[e] <- NULL
  }

  assays["raw"] <- NULL

  # transform assays to matrix
  for(a in names(assays)){
    assays[[a]] <- as.matrix(assays[[a]])
  }

  # get condition
  condition <- get_condition_value(se, condition)

  # prepare data
  avgmadmem <- NormalyzerDE:::calculateAvgMadMem(assays, data.table::as.data.table(SummarizedExperiment::colData(se))[,condition])
  PMAD_intra <- data.table::as.data.table(avgmadmem)
  PMAD_intra <- data.table::melt(PMAD_intra, measure.vars = colnames(PMAD_intra), variable.name = "Normalization", value.name = "PMAD")
  PMAD_intra$Normalization <- factor(PMAD_intra$Normalization, levels = sort(as.character(unique(PMAD_intra$Normalization))))

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  if(diff){
    # check if log in assays
    if(!"log2" %in% names(assays)){
      stop("Log2 data not in SummarizedExperiment! Difference not applicable.")
    }
    avgmadmempdiff<- NormalyzerDE:::calculatePercentageAvgDiffInMat(avgmadmem)
    PMAD_diff <- data.table::data.table(Normalization = names(assays), PMAD = avgmadmempdiff)
    sort_methods <- sort(as.character(unique(PMAD_diff$Normalization)))
    PMAD_diff$Normalization <- factor(PMAD_diff$Normalization, levels=sort_methods)
    p <- ggplot2::ggplot(PMAD_diff, ggplot2::aes(x=get("Normalization"), y=get("PMAD"),fill=get("Normalization"))) +
      ggplot2::geom_col() +
      ggplot2::ylab("") +
      ggplot2::xlab("Normalization Method") +
      ggplot2::geom_label(ggplot2::aes(label=round(get("PMAD"), digits = 0)), fill="white", show.legend = FALSE, color="black") +
      ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  } else {
    p <- ggplot2::ggplot(PMAD_intra, ggplot2::aes(x=get("Normalization"), y=get("PMAD"), fill=get("Normalization"))) +
      ggplot2::geom_boxplot() +
      ggplot2::stat_boxplot(geom ='errorbar', width = 0.4) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
      ggplot2::ylab("PMAD") +
      ggplot2::xlab("Normalization Method") +
      ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector)
  }
  return(p)
}


#' Plot intragroup pooled estimate of variance (PEV) of the normalized data
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param diff Boolean indicating whether to visualize the reduction of intragroup variation (PEV) compared to "log" (TRUE) or a normal boxplot of intragroup variation (PEV) for each normalization method (FALSE).
#'
#' @return ggplot object (boxplot)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' plot_intragroup_PEV(tuberculosis_TMT_se, ain = NULL,
#'                             condition = NULL, diff = FALSE)
#'
plot_intragroup_PEV <- function(se, ain=NULL, condition = NULL, diff = FALSE){
  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }
  assays <- SummarizedExperiment::assays(se)
  exclude <- names(assays)[!names(assays) %in% ain]
  for(e in exclude){
    assays[e] <- NULL
  }
  assays["raw"] <- NULL

  # transform assays to matrix
  for(a in names(assays)){
    assays[[a]] <- as.matrix(assays[[a]])
  }

  # get condition
  condition <- get_condition_value(se, condition)

  # prepare data
  avgVarianceMat <- NormalyzerDE:::calculateAvgReplicateVariation(assays, data.table::as.data.table(SummarizedExperiment::colData(se))[[condition]])
  PEV_intra <- data.table::as.data.table(avgVarianceMat)
  PEV_intra <- data.table::melt(PEV_intra, measure.vars = colnames(PEV_intra), variable.name = "Normalization", value.name = "PEV")
  PEV_intra$Normalization <- factor(PEV_intra$Normalization, levels = sort(as.character(unique(PEV_intra$Normalization))))

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  if(diff){
    # check if log in assays
    if(!"log2" %in% names(assays)){
      stop("Log2 data not in SummarizedExperiment! Difference not applicable.")
    }
    avgvarmempdiff <- NormalyzerDE:::calculatePercentageAvgDiffInMat(avgVarianceMat)
    PEV_diff <- data.table::data.table(Normalization = names(assays), PEV = avgvarmempdiff)
    sort_methods <- sort(as.character(unique(PEV_diff$Normalization)))
    PEV_diff$Normalization <- factor(PEV_diff$Normalization, levels=sort_methods)
    p <- ggplot2::ggplot(PEV_diff, ggplot2::aes(x=get("Normalization"), y=get("PEV"),fill=get("Normalization"))) +
      ggplot2::geom_col() +
      ggplot2::ylab("") +
      ggplot2::xlab("Normalization Method") +
      ggplot2::geom_label(ggplot2::aes(label=round(get("PEV"), digits = 0)), fill="white", show.legend = FALSE, color="black") +
      ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  } else {
    p <- ggplot2::ggplot(PEV_intra, ggplot2::aes(x=get("Normalization"), y=get("PEV"), fill=get("Normalization"))) +
      ggplot2::geom_boxplot() +
      ggplot2::stat_boxplot(geom ='errorbar', width = 0.4) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
      ggplot2::ylab("PEV") +
      ggplot2::xlab("Normalization Method") +
      ggplot2::scale_fill_manual(name = "Normalization Method", values = col_vector)
  }
  return(p)
}
