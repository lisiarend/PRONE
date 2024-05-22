
#' Volcano plots of DE results
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param facet_norm Boolean indicating whether to facet by normalization method (TRUE) or not (FALSE)
#' @param facet_comparison Boolean indicating whether to facet by comparison (TRUE) or not (FALSE). Only valid if facet_norm = FALSE.
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' plot_volcano_DE(tuberculosis_TMT_de_res, ain = NULL,
#'                 comparisons = NULL, facet_norm = TRUE,
#'                 facet_comparison = FALSE)
#'
plot_volcano_DE <- function(de_res, ain = NULL, comparisons = NULL, facet_norm = TRUE, facet_comparison = FALSE){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]
  de_res <- de_res[de_res$Assay %in% ain,]
  # plot
  p <- list()
  if(facet_norm){
    # facet by normalization method
    for(comp in comparisons){
      dt <- de_res[de_res$Comparison == comp,]
      if("Significant Change" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
      } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
      } else {
        color_values <- c("No Change" = "grey")
      }
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("logFC"), y=-log10(get("P.Value")))) +
        ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(ggplot2::aes(col=get("Change"))) +
        ggplot2::labs(y = expression(-log[10]~"P-value"), x= "logFC") +
        ggplot2::scale_color_manual(values = color_values, name = "Change") +
        ggplot2::facet_wrap(~Assay, scales = "free")
      p[[comp]] <- tmp
    }
  } else if (facet_comparison){
    # facet by comparison
    for(method in ain){
      dt <- de_res[de_res$Assay == ain,]
      if("Significant Change" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
      } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
      } else {
        color_values <- c("No Change" = "grey")
      }
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("logFC"), y=-log10(get("P.Value")))) +
        ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(ggplot2::aes(col=get("Change"))) +
        ggplot2::labs(y = expression(-log[10]~"P-value"), x= "logFC") +
        ggplot2::scale_color_manual(values = color_values, name = "Change") +
        ggplot2::facet_wrap(~Comparison, scales = "free")
      p[[method]] <- tmp
    }
  } else {
    # no facet at all --> individual plots for each method and each comparison
    for(comp in comparisons){
      for(method in ain){
        dt <- de_res[de_res$Assay == ain,]
        dt <- dt[dt$Comparison == comp,]
        if("Significant Change" %in% dt$Change){
          color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
        } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
          color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
        } else {
          color_values <- c("No Change" = "grey")
        }
        tmp <- ggplot2::ggplot(dt, ggplot2::aes(x=get("logFC"), y=-log10(get("P.Value")))) +
          ggplot2::geom_vline(xintercept=0) +
          ggplot2::geom_point(ggplot2::aes(col=get("Change"))) +
          ggplot2::labs(y = expression(-log[10]~"P-value"), x= "logFC") +
          ggplot2::scale_color_manual(values = color_values, name = "Change")
        p[[paste0(comp, "_", method)]] <- tmp
      }
    }
  }
  return(p)
}

#' Heatmap of DE results
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set (including the normalized intensities)
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparison String of comparison (must be a valid comparison saved in de_res)
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#' @param pvalue_column column name of p-values in de_res
#'
#' @return list of ComplexHeatmaps for each method
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' data(tuberculosis_TMT_de_res)
#' plot_heatmap_DE(tuberculosis_TMT_se, tuberculosis_TMT_de_res, ain = NULL,
#'                 comparison = "PTB-HC",
#'                 condition = NULL, label_by = NULL,
#'                 pvalue_column = "adj.P.Val")
#'
plot_heatmap_DE <- function(se, de_res, ain, comparison, condition = NULL, label_by = NULL, pvalue_column = "adj.P.Val"){
  # check parameters
  if(length(comparison) != 1){
    stop("Only one comparison as input!")
  }
  ain <- check_input_assays(se, ain)
  tmp <- check_plot_DE_parameters(de_res, ain, comparison)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparison <- tmp[[3]]
  tmp <- get_label_value(se, label_by)
  show_sample_names <- tmp[[1]]
  label_by <- tmp[[2]]
  condition <- get_condition_value(se, condition)
  if(!pvalue_column %in% colnames(de_res)){
    stop("No column named ", pvalue_column, " in DE results. Consider changing the pvalue_column parameter!")
  }
  # TODO: check if comparison matches condition column
  # subset data
  de_res <- de_res[de_res$Comparison == as.character(comparison),]
  # plot
  plots <- list()
  for(method in ain){
    # extract data
    dt <- de_res[de_res$Assay == method,]
    # extract metadata
    cond_a <- strsplit(as.character(comparison), split="-")[[1]][1]
    cond_b <- strsplit(as.character(comparison), split="-")[[1]][2]
    md <- data.table::as.data.table(SummarizedExperiment::colData(se))
    md <- md[md[[condition]] %in% c(cond_a, cond_b),]
    # extract significant proteins
    sig <- dt[ dt$Change != "No Change",]
    if(nrow(sig) > 0){
      # extract intensities of significant proteins
      data <- data.table::as.data.table(SummarizedExperiment::assays(se)[[method]])
      data$IDs <- rownames(data)
      data <- data[data$IDs %in% sig$IDs,]
      data <- data[, c("IDs",md$Column)]
      rownames(data) <- data$ID
      data$IDs <- NULL
      if(show_sample_names){
        colnames(data) <- md[, label_by]
      }
      # color vector
      qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      col_vector <- rev(col_vector)

      # TODO: add protein groups or gene names (instead of ID numbers)
      # Heatmap annotation
      colors <- col_vector[seq_len(length(unique(md[,condition])))]
      names(colors) <- unique(md[, condition])
      col_annotation <- ComplexHeatmap::columnAnnotation(
        Condition = md[, condition],
        annotation_legend_param = list(Condition = list(title = condition)),
        col = list(Condition = colors)
      )
      row_annotation <- ComplexHeatmap::rowAnnotation(
        P.Value = sig[[pvalue_column]],
        annotation_legend_param = list(P.Value = list(title = pvalue_column)),
        col = list(P.Value = circlize::colorRamp2(c(min(sig[[pvalue_column]]), max(sig[[pvalue_column]])), c("green", "white")))
      )

      # Cluster
      cluster_cols <- FALSE
      cluster_rows <- FALSE
      tryCatch({
        tmp <- stats::hclust(stats::dist(as.matrix(data)))
        cluster_cols <- TRUE
        cluster_rows <- TRUE
      }, error = function(e){
        print(e)
      })

      p <- ComplexHeatmap::Heatmap(
        as.matrix(data),
        top_annotation = col_annotation,
        right_annotation = row_annotation,
        show_column_names = TRUE,
        show_row_names = TRUE,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_cols,
        name = paste0("Intensity (", method, ")")
      )
      plots[[method]] <- p

    } else {
      warning(paste0("No significant changes for method ", method, " and comparison ", comparison, ": nothing to plot."))
    }
  }
  return(plots)
}

#' Overview plots of DE results
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param plot_type String indicating whether to plot a single plot per comparison ("single"), facet by comparison ("facet_comp"), stack the number of DE per comparison ("stacked"), or stack the number of DE per comparison but facet by up- and down-regulated ("facet_regulation")
#'
#' @return list of ggplot objects or single object if plot_type = facet or stacked
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' plot_overview_DE_bar(tuberculosis_TMT_de_res, ain = NULL, comparisons = NULL,
#'                      plot_type = "facet_regulation")
#'
plot_overview_DE_bar <- function(de_res, ain = NULL, comparisons = NULL, plot_type = "single"){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]
  stopifnot(plot_type %in% c("single", "facet_comp", "stacked", "facet_regulation"))

  de_res <- de_res[de_res$Comparison %in% comparisons,]

  # get overview DE
  dt <- get_overview_DE(de_res)
  dt <- dt[dt$Assay %in% ain, ]
  melted_dt <- data.table::melt(dt, id.vars = c("Assay", "Comparison"), variable.name = "Change", value.name = "N")
  melted_dt$Assay <- factor(melted_dt$Assay, levels = sort(as.character(unique(melted_dt$Assay))))
  # plot
  p <- list()
  if(plot_type == "stacked"){
    # sum up- and down-regulated proteins
    plot_dt <- melted_dt %>% dplyr::group_by(Comparison, Assay) %>% dplyr::summarise(N = sum(N))
    p <-  ggplot2::ggplot(plot_dt, ggplot2::aes(x = get("N"), y = get("Assay"), fill = get("Comparison"), label = get("N"))) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_stack()) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Comparison") +
      ggplot2::labs(title = "Overview of DE results", x = "Number of proteins", y = "Assay") +
      ggplot2::theme_minimal()
  } else if(plot_type == "facet_regulation"){
    p <-  ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("N"), y = get("Assay"), fill = get("Comparison"), label = get("N"))) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_stack()) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Comparison") +
      ggplot2::facet_wrap(~Change, scales = "free_y") +
      ggplot2::labs(title = "Overview of DE results", x = "Number of proteins", y = "Assay") +
      ggplot2::theme_minimal()
  } else if(plot_type == "facet_comp"){
    if("Significant Change" %in% melted_dt$Change){
      color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
    } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
      color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
    } else {
      color_values <- c("No Change" = "grey")
    }
    p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("N"), y = get("Assay"), fill = get("Change"), label = get("N"))) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_stack()) +
      ggplot2::facet_wrap(~Comparison, scales = "free") +
      ggplot2::scale_fill_manual(values = color_values, name = "Change") +
      ggplot2::labs(title = "Overview of DE results", x = "Number of proteins", y = "Assay") +
      ggplot2::theme_minimal()
  } else {
    for(comp in comparisons){
      dt <- melted_dt[melted_dt$Comparison == comp,]
      if("Significant Change" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Significant Change" = "#D55E00")
      } else if ("Up Regulated" %in% dt$Change | "Down Regulated" %in% dt$Change){
        color_values <- c("No Change" = "grey", "Up Regulated" = "#D55E00", "Down Regulated" =  "#0072B2")
      } else {
        color_values <- c("No Change" = "grey")
      }
      tmp <- ggplot2::ggplot(dt, ggplot2::aes(x = get("N"), y = get("Assay"), fill = get("Change"), label = get("N"))) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_stack()) +
        ggplot2:: scale_fill_manual(values = color_values, name = "Change") +
        ggplot2::labs( y = "Normalization Method", x = "Number of DE Proteins")
      p[[comp]] <- tmp
    }
  }
  return(p)
}


#' Barplot of coverage of DE markers per normalization method in any comparison. (If you want to have a look at a specific comparison, just subset the de_res data table before plotting.)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param markers vector of the IDs of the markers to plot
#' @param id_column String specifying the column of the rowData of the SummarizedExperiment object which includes the IDs of the markers
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' data(tuberculosis_TMT_de_res)
#' plot_coverage_DE_markers(tuberculosis_TMT_se, tuberculosis_TMT_de_res,
#'                           ain = NULL, markers = c("Q7Z7F0", "Q13790"),
#'                           id_column = "Protein.IDs")
#'
plot_coverage_DE_markers <- function(se, de_res, ain, markers, id_column = "Protein.IDs"){
  # check parameters
  ain <- check_input_assays(se, ain)
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons = NULL)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]

  # check id_column
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  if(!id_column %in% colnames(rd)){
    # if id_column not in rowData
    stop(paste0(id_column, " not in rowData of the SummarizedExperiment object!"))
  }

  # prepare data
  de_res <- de_res[de_res$Change != "No Change",]
  de_res <- de_res[de_res$Assay %in% ain,]
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  found_markers <- stringr::str_detect(de_res[[id_column]], paste(markers, collapse = "|") )
  de_biomarkers <- de_res[found_markers,]
  de_biomarkers <- de_biomarkers[, c("Assay", id_column), with = FALSE]
  de_biomarkers <- unique(de_biomarkers)
  stats <- de_biomarkers %>% dplyr::group_by(Assay) %>% dplyr::summarise(DEPs = dplyr::n())
  if(nrow(de_biomarkers)==0){
    stats <- data.table::data.table(Assay = unique(de_res$Assay), DEPs = rep(0, nrow(stats)))
  }
  missing_methods <- unique(de_res$Assay)[! unique(de_res$Assay) %in% stats$Assay]
  if(length(missing_methods) > 0){
    tmp <- data.table::data.table(Assay = missing_methods, DEPs = rep(0, length(missing_methods)))
    stats <- rbind(stats, tmp)
  }
  stats$Coverage <- stats$DEPs / length(markers)
  stats$Label <- paste0(stats$DEPs, "/", length(markers))

  # color vector
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = get("Assay"), y = get("Coverage"), fill = get("Assay"))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "Assay", y = "Marker Coverage") +
    ggplot2::scale_fill_manual(values = col_vector, name = "Normalization Method") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
    ggplot2::geom_label(ggplot2::aes(label = get("Label")))
  return(p)
}

#' Overview heatmap plot of DE results
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' plot_overview_DE_tile(tuberculosis_TMT_de_res, ain = NULL,
#'                       comparisons = NULL)
#'
plot_overview_DE_tile <- function(de_res, ain = NULL, comparisons = NULL){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  # get overview DE
  dt <- de_res[de_res$Change != "No Change",]
  dt <- dt[dt$Assay %in% ain, ]
  dt <- dt[dt$Comparison %in% comparisons,]
  dt <- dt %>% dplyr::count(Comparison, Assay, sort = TRUE)

  # plot
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = get("Assay"), y = get("Comparison"), fill = get("n"))) +
    ggplot2::geom_tile() +
    ggplot2:: theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1)) +
    ggplot2::labs(x = "Normalization Method", y = "Comparison") +
    ggplot2::scale_fill_distiller(name = "DE Proteins", direction = 1)
  return(p)
}



