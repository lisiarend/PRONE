
#' Get performance metrics of DE results of spike-in data set.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#'
#' @return data table with multiple performance metrics of the DE results
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' stats <- get_spiked_stats_DE(spike_in_se, spike_in_de_res)
#'
get_spiked_stats_DE <- function(se, de_res){
  spike_column <- S4Vectors::metadata(se)$spike_column
  spike_val <- S4Vectors::metadata(se)$spike_value
  stats <- de_res %>% dplyr::group_by(Assay, Comparison) %>%
    dplyr::summarise(TP = sum(Change %in% c("Up Regulated", "Down Regulated", "Significant Change") & get(spike_column) == spike_val, na.rm=TRUE),
                     FP = sum(Change %in% c("Up Regulated", "Down Regulated", "Significant Change") & get(spike_column) != spike_val, na.rm=TRUE),
                     FN = sum(Change == "No Change" & get(spike_column) == spike_val, na.rm=TRUE),
                     TN = sum(Change == "No Change" & get(spike_column) != spike_val, na.rm=TRUE),
                     Sensitivity = TP / (TP + FN),
                     Specificity = TN /(TN + FP),
                     Precision = TP / (TP + FP),
                     FPR = FP / (FP + TN),
                     F1Score = 2 * (Precision * Sensitivity) / (Precision + Sensitivity),
                     Accuracy = (TP + TN)/(TP + TN + FP + FN)) %>%
    data.table::as.data.table()
  return(stats)
}


#' Heatmap of performance metrics for spike-in data sets
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#' @param metrics vector of Strings specifying the metrics (must be colnames of stats)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' stats <- get_spiked_stats_DE(spike_in_se, spike_in_de_res)
#' plot_stats_spiked_heatmap(stats, ain = NULL, comparisons = NULL,
#'                           metrics = c("F1Score", "Accuracy"))
#'
plot_stats_spiked_heatmap <- function(stats, ain = NULL, comparisons = NULL, metrics = c("Accuracy", "Precision", "F1Score")){
  # check if measures in stats
  if(is.null(metrics)){
    message("All available metrics will be used for plotting")
    valid_metrics <- colnames(stats)[!colnames(stats) %in% c("TP", "FN", "TN", "FP")]
  } else {
    metrics <- c(metrics)
    valid_metrics <- metrics[metrics %in% colnames(stats)]
    non_valid_metrics <- metrics[!metrics %in% colnames(stats)]
    if(length(valid_metrics)== 0){
      stop("No valid metrics! Check colnames of stats.")
    } else if(length(non_valid_metrics)>0){
      warning(paste0("Not all metrics are valid!", metrics, " are not visualized!"))
    }
  }
  metrics <- valid_metrics
  dt <- stats[, colnames(stats) %in% c("Assay", "Comparison", metrics)]

  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]

  # prepare data
  melted_dt <- data.table::melt(dt, variable.name = "Metric", value.name = "Value", id.vars = c("Assay", "Comparison"))
  p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("Assay"), y = get("Comparison"), fill = get("Value"))) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::labs(x = "Normalization Method", y = "Comparison") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 2, barheight = 5, title = "Value"))
  # plot
  if(length(metrics) == 1){
   p <- p + ggplot2::scale_fill_distiller(name = metrics[1],palette = "YlOrRd", direction = 1, limits = c(0,1))
  } else {
    p <- p + ggplot2::facet_wrap(~Metric, ncol = 2) + ggplot2::scale_fill_distiller(name = "Measure",palette = "YlOrRd", direction = 1, limits = c(0,1))
  }
  return(p)
}


#' Barplot of true and false positives for specific comparisons and normalization methods
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object (barplot)
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' stats <- get_spiked_stats_DE(spike_in_se, spike_in_de_res)
#' plot_TP_FP_spiked_bar(stats, ain = NULL, comparisons = NULL)
#'
plot_TP_FP_spiked_bar <- function(stats, ain = NULL, comparisons = NULL){
  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- stats[, c("Assay", "Comparison", "TP", "FP")]
  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]

  dt$TP <- as.integer(dt$TP)
  dt$FP <- as.integer((-1) * as.integer(dt$FP))

  melted_dt <- data.table::melt(dt, variable.name = "Class", value.name = "Value", measure.vars = c("TP", "FP"))

  # plot
  p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x=get("Assay"), y=get("Value"), fill = get("Class"))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(y="False Positives | True Positives", x = "Normalization Method") +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "Class") +
    ggplot2::facet_wrap(~Comparison, scales="free_y") +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::coord_flip()
  return(p)
}

#' Boxplot of true and false positives for specific comparisons and normalization methods
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object (barplot)
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' stats <- get_spiked_stats_DE(spike_in_se, spike_in_de_res)
#' plot_TP_FP_spiked_box(stats, ain = NULL, comparisons = NULL)
#'
plot_TP_FP_spiked_box <- function(stats, ain = NULL, comparisons = NULL){
  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- stats[, c("Assay", "Comparison", "TP", "FP")]
  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]

  dt$TP <- as.integer(dt$TP)
  dt$FP <- as.integer(dt$FP)

  melted_dt <- data.table::melt(dt, variable.name = "Class", value.name = "Value", measure.vars = c("TP", "FP"))

  # order methods
  tps <- melted_dt[melted_dt$Class == "TP",]
  meds <- tps %>% dplyr::group_by(Assay) %>% dplyr::summarise(Median = stats::median(Value, na.rm = TRUE)) %>% data.table::as.data.table()
  meds <- meds[order(meds$Median),]
  melted_dt_1 <- melted_dt
  melted_dt_1$Assay <- factor(melted_dt_1$Assay, levels = meds$Assay)

  # plot
  p <- ggplot2::ggplot(melted_dt_1, ggplot2::aes(x = get("Assay"), y = get("Value"), fill = get("Class"))) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "Class") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle =90, vjust =0.5)) +
    ggplot2::labs(x = "Normalization Method", y="Number of Proteins")
  return(p)
}

#' Scatterplot of true positives and false positives (median with errorbars as Q1, and Q3) for all comparisons
#'
#' @param stats data table with multiple metrics of the DE results (resulting of get_spiked_stats_DE)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' stats <- get_spiked_stats_DE(spike_in_se, spike_in_de_res)
#' plot_TP_FP_spiked_scatter(stats, ain = NULL, comparisons = NULL)
#'
plot_TP_FP_spiked_scatter <- function(stats, ain = NULL, comparisons = NULL){
  # check stats, ain, comparisons
  tmp <- check_stats_spiked_DE_parameters(stats, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  dt <- stats[, c("Assay", "Comparison", "TP", "FP")]
  dt <- dt[dt$Assay %in% ain,]
  dt <- dt[dt$Comparison %in% comparisons,]
  dt$TP <- as.integer(dt$TP)
  dt$FP <- as.integer(dt$FP)
  summarized_stats <- dt %>% dplyr::group_by(Assay) %>% dplyr::summarize(Median_TP = stats::median(TP, na.rm = TRUE), Median_FP = stats::median(FP, na.rm = TRUE), Q1_TP = stats::quantile(TP, c(0.25), na.rm = TRUE), Q3_TP = stats::quantile(TP,c(0.75), na.rm = TRUE),Q1_FP = stats::quantile(FP, c(0.25), na.rm = TRUE), Q3_FP = stats::quantile(FP,c(0.75), na.rm = TRUE))

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  p <- ggplot2::ggplot(summarized_stats, ggplot2::aes(x = Median_TP, y = Median_FP, color = Assay)) +
    ggplot2::geom_point(size = 6, shape = 15) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Q1_FP, ymax = Q3_FP), linewidth = 1) +
    ggplot2::scale_y_continuous(trans = "log10", labels = scales::scientific) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Q1_TP, xmax = Q3_TP), linewidth = 1) +
    ggplot2::scale_x_continuous(trans = "log10", labels = scales::scientific) +
    ggplot2::scale_color_manual(name = "Normalization", values = col_vector) +
    ggplot2::labs(x = "log10(TPs)", y = "log10(FPs)")
  return(p)
}

#' Boxplot of log fold changes of spike-in and background proteins for specific normalization methods and comparisons. The ground truth (calculated based on the concentrations of the spike-ins) is shown as a horizontal line.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' plot_fold_changes_spiked(spike_in_se, spike_in_de_res,
#'                          condition = "Condition", ain = NULL,
#'                          comparisons = NULL)
#'
plot_fold_changes_spiked <- function(se, de_res, condition, ain = NULL, comparisons = NULL){
  # check input
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  de_res <- de_res[de_res$Assay %in% ain,]
  de_res <- de_res[de_res$Comparison %in% comparisons,]

  # check if spike_column in de_res
  spike_column <- S4Vectors::metadata(se)$spike_column
  if(!spike_column %in% colnames(de_res)){
    stop(paste0(spike_column, " not in de_res. Please perform run_DE again!"))
  }

  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))

  # check condition
  condition <- get_condition_value(se, condition)

  # retrieve concentration
  spike_concentration <- S4Vectors::metadata(se)$spike_concentration
  concentrations <- unique(coldata[, c(condition, spike_concentration)])

  comps <- data.table::data.table(Comparison = unique(de_res$Comparison))
  comps$SampleA <- sapply(strsplit(as.character(comps$Comparison),"-"), "[", 1)
  comps$SampleB <- sapply(strsplit(as.character(comps$Comparison),"-"), "[", 2)
  comps <- merge(comps, concentrations, by.x = "SampleA", by.y = condition)
  comps <- comps %>% dplyr::rename("ConcentrationA" = Concentration)
  comps <- merge(comps, concentrations, by.x = "SampleB", by.y = condition)
  comps <- comps %>% dplyr::rename("ConcentrationB" = Concentration)
  comps$Truth <- log2(comps$ConcentrationA / comps$ConcentrationB)

  spike_value <- S4Vectors::metadata(se)$spike_value
  types <- unique(de_res[[spike_column]])
  levels <- c(spike_value, types[!types %in% spike_value])
  de_res[[spike_column]] <- factor(de_res[[spike_column]], levels = levels)

  de_res <- merge(de_res, comps[, c("Comparison", "Truth")], by="Comparison")

  p <- ggplot2::ggplot(de_res, ggplot2::aes(x=get("Assay"), y=get("logFC"), fill=get(spike_column))) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~Comparison, ncol=1, scales = "free_y") +
    ggplot2::labs(fill = spike_column, x = "Normalization Method", y = "LogFC") +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2")) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour="#0072B2", alpha=0.8)  +
    ggplot2::geom_hline(ggplot2::aes(yintercept = get("Truth")), linetype = "dashed", color="#D55E00", alpha=0.8)
  return(p)
}


#' Boxplot of p-values of spike-in and background proteins for specific normalization methods and comparisons. The ground truth (calculated based on the concentrations of the spike-ins) is shown as a horizontal line.
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' plot_pvalues_spiked(spike_in_se, spike_in_de_res, ain = NULL,
#'                     comparisons = NULL)
#'
plot_pvalues_spiked <- function(se, de_res, ain = NULL, comparisons = NULL){
  # check input
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  de_res <- de_res[de_res$Assay %in% ain,]
  de_res <- de_res[de_res$Comparison %in% comparisons,]

  # check if spike_column in de_res
  spike_column <- S4Vectors::metadata(se)$spike_column
  if(!spike_column %in% colnames(de_res)){
    stop(paste0(spike_column, " not in de_res. Please perform run_DE again!"))
  }

  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))

  spike_value <- S4Vectors::metadata(se)$spike_value
  types <- unique(de_res[[spike_column]])
  levels <- c(spike_value, types[!types %in% spike_value])
  de_res[[spike_column]] <- factor(de_res[[spike_column]], levels = levels)


  p <- ggplot2::ggplot(de_res, ggplot2::aes(x=get("Assay"), y=log10(get("P.Value")), fill=get(spike_column))) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~Comparison, ncol=1, scales = "free_y") +
    ggplot2::labs(fill = spike_column, x = "Normalization Method", y = expression(-log[10]~"P-value")) +
    ggplot2::scale_fill_manual(values = c("#D55E00", "#0072B2"))
  return(p)
}


#' Line plot of number of true and false positives when applying different logFC thresholds
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#' @param nrow number of rows for facet wrap
#' @param alpha threshold for adjusted p-values
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' plot_logFC_thresholds_spiked(spike_in_se, spike_in_de_res,
#'                              condition = "Condition", ain = NULL,
#'                              comparisons = NULL, nrow = 2, alpha = 0.05)
#'
plot_logFC_thresholds_spiked <- function(se, de_res, condition, ain = NULL, comparisons = NULL, nrow = 2, alpha = 0.05){
  # check input
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  de_res <- de_res[de_res$Assay %in% ain,]
  de_res <- de_res[de_res$Comparison %in% comparisons,]

  # check condition
  condition <- get_condition_value(se, condition)

  # retrieve concentrations to calculate ground truth logFC per comparison
  coldata <- data.table::as.data.table(SummarizedExperiment::colData(se))
  concentration <- S4Vectors::metadata(se)$spike_concentration
  condition <- S4Vectors::metadata(se)$condition
  concentrations <- unique(coldata[, c(condition, concentration), with=FALSE])
  comps <- data.table::data.table(Comparison = unique(de_res$Comparison))
  comps$SampleA <- sapply(strsplit(as.character(comps$Comparison),"-"), "[", 1)
  comps$SampleB <- sapply(strsplit(as.character(comps$Comparison),"-"), "[", 2)
  comps <- merge(comps, concentrations, by.x = "SampleA", by.y = condition)
  comps <- comps %>% dplyr::rename("ConcentrationA" = Concentration)
  comps <- merge(comps, concentrations, by.x = "SampleB", by.y = condition)
  comps <- comps %>% dplyr::rename("ConcentrationB" = Concentration)
  comps$Truth <- log2(comps$ConcentrationA / comps$ConcentrationB)

  step <- ifelse(round(max(comps$Truth)) < 0, -0.5, 0.5)
  logFCs <- seq(0, round(max(comps$Truth)), step)
  alpha <- 0.05

  stats <- NULL
  for(thr in logFCs){
    de_res$Change <- ifelse(de_res$adj.P.Val <= alpha & de_res$logFC >= thr , "Significant Change", "No Change")
    stats_chunk <- get_spiked_stats_DE(se, de_res)
    stats_chunk$logFCthr <- thr
    if(is.null(stats)){
      stats <- stats_chunk
    } else {
      stats <- rbind(stats, stats_chunk)
    }
  }

  # add ground truth
  stats <- merge(stats, comps[, c("Comparison", "Truth")])

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual",]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  # plot
  p_tp <- ggplot2::ggplot(stats, ggplot2::aes(x = get("logFCthr"), y= get("TP"), color = get("Assay"))) +
    ggplot2::geom_line(lwd = 1) +
    ggplot2::facet_wrap(~Comparison, scales = "free_y", nrow = nrow) +
    ggplot2::labs(x = "Fold Change Threshold", y="Number of True Positives") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = get("Truth")), linetype="dotted") +
    ggplot2::scale_x_continuous(breaks = seq(0, round(max(stats$Truth)))) +
    ggplot2::scale_color_manual(name = "Normalization Method", values = col_vector)

  p_fp <- ggplot2::ggplot(stats, ggplot2::aes(x = get("logFCthr"), y= get("FP"), color = get("Assay"))) +
    ggplot2::geom_line(lwd =1) +
    ggplot2::facet_wrap(~Comparison, scales = "free_y", nrow = nrow) +
    ggplot2::labs(x = "Fold Change Threshold", y="Number of False Positives") +
    ggplot2::scale_x_continuous(breaks = seq(0, round(max(stats$Truth)))) +
    ggplot2::scale_color_manual(name = "Normalization Method", values = col_vector)
  return(list("TP" = p_tp, "FP" = p_fp))
}


#' Plot ROC curve and barplot of AUC values for each method for a specific comparion or for all comparisons
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in stats)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in stats)
#'
#' @return list of ggplot objects
#' @export
#'
#' @examples
#' data(spike_in_se)
#' data(spike_in_de_res)
#' plot_ROC_AUC_spiked(spike_in_se, spike_in_de_res)
#'
plot_ROC_AUC_spiked <- function(se, de_res, ain = NULL, comparisons = NULL){
  # check input
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  stats <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]

  de_res <- de_res[de_res$Assay %in% ain,]
  de_res <- de_res[de_res$Comparison %in% comparisons,]

  # extract adjusted p-values from limma
  p.adj <- data.table::data.table("Protein.IDs" = de_res$Protein.IDs, "P.Adj" = de_res$adj.P.Val, "Assay" = de_res$Assay, "Comparison" = de_res$Comparison)
  # extract ground truth
  spike <- S4Vectors::metadata(se)$spike_column
  spike_val <- S4Vectors::metadata(se)$spike_value
  rowdata <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  truth <- data.table::data.table(Protein.IDs = rowdata$Protein.IDs, spike = rowdata[[spike]])
  colnames(truth) <- c("Protein.IDs", spike)
  truth$Truth <- ifelse(truth[[spike]] == spike_val, 1, 0)

  dt <- merge(p.adj, truth, by ="Protein.IDs")
  dt <- dt[order(dt$P.Adj, decreasing = TRUE),]

  cuts <- 0
  round <- 3
  ymin <- 0
  ymax <- 1

  dt$P.Adj <- 1 - as.numeric(dt$P.Adj)

  # colors
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- rev(col_vector)

  roc <- ggplot2::ggplot(dt, ggplot2::aes(d=get("Truth"), m=get("P.Adj"), color = Assay)) +
    plotROC::geom_roc(n.cuts = cuts, labelround = round) +
    ggplot2::ylim(ymin, ymax) +
    ggplot2::xlab ("FPR") +
    ggplot2::ylab("TPR") +
    ggplot2::scale_color_manual(name = "Normalization", values =col_vector) +
    ggplot2::facet_wrap(~Comparison)

  auc_dt <- plotROC::calc_auc(roc)
  auc_bars <- ggplot2::ggplot(auc_dt, ggplot2::aes(x=get("Assay"), y=get("AUC"), fill=get("Assay"))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~Comparison) +
    ggplot2::scale_fill_manual(name = "Normalization", values = col_vector) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5)) +
    ggplot2::labs(x = "Normalization")
  auc_box <- ggplot2::ggplot(auc_dt, ggplot2::aes(x = Assay, y=AUC, fill=Assay)) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(name = "Normalization", values = col_vector) +
    ggplot2::labs(x = "Normalization")

  return(list("ROC" = roc, "AUC_bars" = auc_bars, "AUC_box" = auc_box,"AUC_dt" = auc_dt))
}



