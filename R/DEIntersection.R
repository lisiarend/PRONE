#' Upset plots of DE results of the different normalization methods
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param min_degree Minimal degree of an intersection for it to be included
#' @param plot_type String indicating whether to plot a single plot per comparison ("single") or stack the number of DE per comparison ("stacked)
#'
#' @return list of plots and intersection tables (split by comparison if plot_type == "single)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' plot_upset_DE(tuberculosis_TMT_de_res,
#'               ain = c("IRS_on_RobNorm", "IRS_on_Median"),
#'               comparisons = NULL, min_degree = 2,
#'               plot_type = "stacked")
#'
plot_upset_DE <- function(de_res, ain = NULL, comparisons = NULL, min_degree = 2, plot_type = "single") {
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]
  stopifnot(plot_type %in% c("single", "stacked"))
  de_res <- de_res[de_res$Assay %in% ain, ]
  de_res <- de_res[de_res$Comparison %in% comparisons,]
  p <- list()

  # prepare data for upset
  dt <- de_res[de_res$Change %in% c("Up Regulated", "Down Regulated", "Significant Change"), ]
  dt <- dt[, c("Assay", "Comparison", "Protein.IDs"), with = FALSE]
  dt <- data.table::dcast(dt, Protein.IDs + Comparison ~ Assay, fun.aggregate = length) %>% as.data.frame()

  # prepare table of intersections
  available_ains <- colnames(dt[,!colnames(dt) %in% c("Protein.IDs", "Comparison")])
  nr_methods <- data.frame("Protein.IDs" = dt$Protein.IDs, "Comparison" = dt$Comparison, "Nr" = rowSums(dt[,available_ains]))
  t <- purrr::map2_df(dt, names(dt), ~  replace(.x, .x==TRUE, .y))
  t[t == 0] <- NA
  t <- as.data.frame(t)
  res <- t %>% tidyr::unite(., available_ains, col = "Assays", na.rm=TRUE, sep = ",")
  res <- merge(res, nr_methods, by = c("Protein.IDs", "Comparison"))
  res <- res[order(-res$Nr),]
  res <- res[, c("Protein.IDs", "Nr", "Assays", "Comparison")]
  colnames(res) <- c("Protein.IDs", "Number of Intersected Assays", "Assays", "Comparison")

  if(plot_type == "stacked"){
    upset <- ComplexUpset::upset(dt, intersect = available_ains, min_degree = min_degree,
                                 base_annotations = list("Intersection Size" = ComplexUpset::intersection_size(counts = TRUE,
                                                                                                               bar_number_threshold = 1,
                                                                                                               mapping = ggplot2::aes(fill = Comparison)) + ggplot2::scale_fill_brewer(palette = "Set2", name = "Comparison")))
    p <- list("upset" = upset, "table" = res)
  } else {

    for (comp in comparisons) {
      comp_dt <- dt[dt$Comparison == comp, ]
      comp_res <- res[res$Comparison == comp,]
      if (nrow(dt) == 0) {
        warning(paste0("No significant changes for comparison ", comp, ": nothing to plot."))
      } else {
        # plot
        upset <- ComplexUpset::upset(
          comp_dt,
          available_ains,
          name = "",
          set_sizes = ComplexUpset::upset_set_size(position = "right") + ggplot2::ylab("Set Size"),
          sort_sets = FALSE,
          keep_empty_groups = FALSE,
          sort_intersections = "descending",
          min_degree = min_degree,
          base_annotations = list("Intersection Size" =
                                    ComplexUpset::intersection_size(text = list(size = 3))),
          themes = ComplexUpset::upset_default_themes(text = ggplot2::element_text(size = 12))
        )

        p[[comp]] <- list("upset" = upset, "table" = comp_res)
      }
    }
  }
  return(p)
}



#' Jaccard similarity heatmap of DE proteins of the different normalization methods
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param plot_type String indicating whether to plot a single plot per comparison ("single"), facet by comparison ("facet_comp"), or include all comparisons in a single plot ("all")
#'
#' @return ggplot object (list of objects if plot_type == "single)
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' plot_jaccard_heatmap(tuberculosis_TMT_de_res, ain = NULL,
#'                     comparisons = NULL, plot_type = "all")
#'
plot_jaccard_heatmap <- function(de_res, ain = NULL, comparisons = NULL, plot_type = "single"){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]
  stopifnot(plot_type %in% c("single", "facet_comp", "all"))
  de_res <- de_res[de_res$Assay %in% ain, ]
  p <- list()

  # prepare data for upset
  dt <- de_res[de_res$Change != "No Change",]
  dt <- dt[dt$Comparison %in% comparisons,]
  if(plot_type == "all"){
    # take all comparisons
    dt_cast <- data.table::dcast(dt, Assay ~ Protein.IDs, fun.aggregate = length)
    assays <- dt_cast$Assay
    intersections <- as.data.frame(lapply(dt_cast[, -1], function(x) as.integer(x > 0)))
    intersections_m <- as.matrix(intersections)
    intersections_m <- t(intersections_m)
    intersections_dt <- as.data.frame(intersections_m)
    colnames(intersections_dt) <- assays

    jaccard <- function(x, y){
      M.11 <- sum(x == 1 & y == 1)
      M.10 <- sum(x == 1 & y == 0)
      M.01 <- sum(x == 0 & y == 1)
      return(M.11 / (M.11 + M.10 + M.01))
    }

    m <- matrix(data = NA, nrow = ncol(intersections_dt), ncol = ncol(intersections_dt))
    for(i in seq_len(ncol(intersections_dt))){
      for(j in seq_len(ncol(intersections_dt))){
        col1 <- colnames(intersections_dt)[i]
        col2 <- colnames(intersections_dt)[j]
        if(col1 == col2){
          m[i, j] <- 1
        } else if(i > j){
          m[i, j] <- jaccard(intersections_dt[, col1], intersections_dt[, col2])
        }
      }
    }
    colnames(m) <- colnames(intersections_dt)
    rownames(m) <- colnames(intersections_dt)
    melted_m <- reshape2::melt(m, measure.vars = colnames(m), na.rm = TRUE)

    p <- ggplot2::ggplot(melted_m, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2))) +
      ggplot2::scale_fill_gradient(low = "white", high = "#4daabd") +
      ggplot2::labs(x = "Normalization Method", y = "Normalization Method", fill = "Jaccard Similarity") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  } else{
    # facet by comparison -> need to include all comparisons
    complete_m <- NULL
    for(comp in comparisons){
      dt_comp <- dt[dt$Comparison == comp,]
      # take all comparisons
      dt_cast <- data.table::dcast(dt_comp, Assay ~ Protein.IDs, fun.aggregate = length)
      assays <- dt_cast$Assay
      intersections <- as.data.frame(lapply(dt_cast[, -1], function(x) as.integer(x > 0)))
      intersections_m <- as.matrix(intersections)
      intersections_m <- t(intersections_m)
      intersections_dt <- as.data.frame(intersections_m)
      colnames(intersections_dt) <- assays

      jaccard <- function(x, y){
        M.11 <- sum(x == 1 & y == 1)
        M.10 <- sum(x == 1 & y == 0)
        M.01 <- sum(x == 0 & y == 1)
        return(M.11 / (M.11 + M.10 + M.01))
      }

      m <- matrix(data = NA, nrow = ncol(intersections_dt), ncol = ncol(intersections_dt))
      for(i in seq_len(ncol(intersections_dt))){
        for(j in seq_len(ncol(intersections_dt))){
          col1 <- colnames(intersections_dt)[i]
          col2 <- colnames(intersections_dt)[j]
          if(col1 == col2){
            m[i, j] <- 1
          } else if(i > j){
            m[i, j] <- jaccard(intersections_dt[, col1], intersections_dt[, col2])
          }
        }
      }
      colnames(m) <- colnames(intersections_dt)
      rownames(m) <- colnames(intersections_dt)
      melted_m <- reshape2::melt(m, measure.vars = colnames(m), na.rm = TRUE)
      melted_m$Comparison <- comp
      if(nrow(melted_m) > 1){
        if(is.null(complete_m)){
          complete_m <- melted_m
        } else {
          complete_m <- rbind(complete_m, melted_m)
        }
      } else {
        warning(paste0("Comparison ", comp, " not included in visualization because only 1 normalization approach found DEPs!"))
      }
    }
    if(plot_type == "facet_comp"){
      p <- ggplot2::ggplot(complete_m, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::facet_wrap(~ Comparison, scales = "free") +
        ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2))) +
        ggplot2::scale_fill_gradient(low = "white", high = "#4daabd") +
        ggplot2::labs(x = "Normalization Method", y = "Normalization Method", fill = "Jaccard Similarity") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
    } else {
      p <- list()
      for(comp in comparisons){
        tmp <- complete_m[complete_m$Comparison == comp,]
        if(nrow(tmp) > 1){
          p_comp <- ggplot2::ggplot(tmp, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2))) +
            ggplot2::scale_fill_gradient(low = "white", high = "#4daabd") +
            ggplot2::labs(x = "Normalization Method", y = "Normalization Method", fill = "Jaccard Similarity") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
          p[[comp]] <- p_comp
        }
      }
    }
  }
  return(p)
}


#' Extract consensus DE candidates
#'
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param norm_thr Threshold for the number of normalization methods that must agree on a DE candidate
#' @param per_comparison Logical indicating if the consensus should be calculated per comparison
#'
#' @return data table with consensus DE candidates
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' extract_consensus_DE_candidates(tuberculosis_TMT_de_res, ain = NULL,
#'                   comparisons = NULL, norm_thr = 0.8, per_comparison = TRUE)
#'
extract_consensus_DE_candidates <- function(de_res, ain = NULL, comparisons = NULL, norm_thr = 0.8, per_comparison = FALSE){
  # extract only significant changes
  dt <- de_res[de_res$Change != "No Change",]
  if(!is.null(ain)){
    dt <- dt[dt$Assay %in% ain,]
  } else {
    ain <- unique(dt$Assay)
  }
  # check comparisons
  if(!is.null(comparisons)){
    dt <- dt[dt$Comparison %in% comparisons,]
  }
  if(per_comparison){
    dt <- dt[, c("Protein.IDs", "Assay", "Comparison")]
    count_dt <- dt %>% dplyr::group_by(Protein.IDs, Comparison) %>% dplyr::summarise("Intersections" = dplyr::n())
    count_dt$Percentage <- count_dt$Intersections / length(ain)
    consensus_de <- count_dt[count_dt$Percentage >= norm_thr,]
  } else {
    dt <- dt[, c("Protein.IDs", "Assay")]
    count_dt <- dt %>% dplyr::group_by(Protein.IDs) %>% dplyr::summarise("Intersections" = n())
    count_dt$Percentage <- count_dt$Intersections / length(ain)
    consensus_de <- count_dt[count_dt$Percentage >= norm_thr,]
  }
  consensus_de$Percentage <- consensus_de$Percentage * 100
  return(consensus_de)
}

