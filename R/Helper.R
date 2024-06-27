

#' Helper function to read example data
#'
#' @param path NULL to get all example data set files, otherwise specify the file name
#'
#' @return If path=NULL a character vector with the file names, otherwise the path to the specific file
#' @export
#'
#' @examples
#' readPRONE_example()
readPRONE_example <- function(path = NULL){
  if(is.null(path)){
    dir(system.file("extdata", package = "PRONE"))
  } else {
    system.file("extdata", path, package = "PRONE", mustWork = TRUE)
  }
}

#' Helper function to check the condition value
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#'
#' @return String of column for condition
#'
#' @keywords internal
get_condition_value <- function(se, condition){
  if(is.null(condition)){
    # check if condition in metadata of se
    if(is.null(S4Vectors::metadata(se)$condition)){
      # print error
      stop("No condition provided!")
    } else {
      condition <- S4Vectors::metadata(se)$condition
      message("Condition of SummarizedExperiment used!")
    }
  } else {
    # check if condition in colData
    if(! condition %in% colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))){
      stop(paste0("No column named ", condition, " in SummarizedExperiment!"))
    }
  }
  return(condition)
}


#' Helper function to get correct value for coloration of plots (color_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param color_by String specifying the column to color the samples (If NULL, the condition column of the SummarizedExperiment object is used. If "No", no color bar added.)
#'
#' @return String of column to color or NULL if no color should be applied
#'
#' @keywords internal
get_color_value <- function(se, color_by){
  if(is.null(color_by)){
    # check if condition in metadata of se
    if(is.null(S4Vectors::metadata(se)$condition)){
      color_by <- NULL
      message("No condition provided. Hence, no color bar added.")
    } else {
      color_by <- S4Vectors::metadata(se)$condition
      message("Condition of SummarizedExperiment used!")
    }
  } else if(length(color_by)==1){
    if(color_by == "No"){
      color_by <- NULL
      message("No color bar added.")
    }
  }
  if(!is.null(color_by)){
    #check if color_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(sum(!color_by %in% cols) > 0){
      stop("The color_by value(s) not a valid column in the SummarizedExperiment object!")
    }
  }
  return(color_by)
}



#' Helper function to get correct value for sample labeling of plots (label_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param label_by String specifying the column to label the samples (If NULL, the labels column of the SummarizedExperiment object is used. If "No", no labeling of samples done.)
#'
#' @return String of column to label or NULL if no label should be applied
#'
#' @keywords internal
get_label_value <- function(se, label_by){
  show_sample_names <- TRUE
  if(is.null(label_by)){
    # check if condition in metadata of se
    if(is.null(S4Vectors::metadata(se)$label)){
      label_by <- "Column"
      show_sample_names <- FALSE
      message("No label provided. Hence, no labeling of samples.")
    } else {
      label_by <- S4Vectors::metadata(se)$label
      message("Label of SummarizedExperiment used!")
    }
  } else if (label_by == "No") {
    label_by <- "Column"
    show_sample_names <- FALSE
    message("No labeling of samples.")
  } else {
    # check if label_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!label_by %in% cols){
      stop("Label_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(list(show_sample_names, label_by))
}

#' Helper function to get correct value for shaping of plots (shape_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param shape_by String specifying the column to shape the samples (If NULL or "No", no shaping is done.)
#'
#' @return String of column to shape or NULL if no shaping should be done
#'
#' @keywords internal
get_shape_value <- function(se, shape_by){
  if(is.null(shape_by)){
      message("No shaping done.")
  } else if(shape_by == "No"){
    shape_by <- NULL
    message("No shaping done.")
  } else {
    # check if shape_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!shape_by %in% cols){
      stop("Shape_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(shape_by)
}


#' Helper function to get correct value for faceting of plots (facet_by parameter)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param facet_by String specifying the column to facet the samples (If NULL or "No", no faceting is done.)
#'
#' @return String of column to facet or NULL if no faceting should be done
#'
#' @keywords internal
get_facet_value <- function(se, facet_by){
  if(is.null(facet_by)){
    message("No faceting done.")
  } else if(facet_by == "No"){
    facet_by <- NULL
    message("No faceting done.")
  } else {
    # check if facet_by a valid column of the SummarizedExperiment object
    cols <- colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))
    if(!facet_by %in% cols){
      stop("Facet_by value not a valid column in the SummarizedExperiment object!")
    }
  }
  return(facet_by)
}

#' Helper function to transform an expression data frame to a data table
#'
#' @param expr_data Expression data frame containing the expression data
#' @param column_names Column names of the expression data
#' @param row_names Row names of the expression data
#'
#' @return Data table containing the expression data
#'
#' @keywords internal
expressToDT <- function(expr_data, column_names, row_names) {
  data_df <- data.table::as.data.table(t(data.frame(expr_data)))
  colnames(data_df) <- column_names
  rownames(data_df) <- row_names
  return(data_df)
}


#' Helper function to transform a tibble to a data table
#'
#' @param expr_data Tibble data frame containing the expression data
#' @param column_names Column names of the expression data
#' @param row_names Row names of the expression data
#'
#' @return Data table containing the expression data
#'
#' @keywords internal
tibToDF <- function(expr_data, column_names, row_names) {
  data_df <- data.table::as.data.table(expr_data)
  colnames(data_df) <- column_names
  rownames(data_df) <- row_names
  return(data_df)
}


#' Helper function to check whether all given assays are in SummarizedExperiment object
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#'
#' @return NULL if no methods in SummarizedExperiment object, else all available methods ready for visualization
#'
#' @keywords internal
check_input_assays <- function(se, ain) {
  if (is.null(ain)) {
    assays_se <- names(SummarizedExperiment::assays(se))
    message("All assays of the SummarizedExperiment will be used.")
    return(assays_se)
  } else {
    assays_se <- names(SummarizedExperiment::assays(se))
    not_existing_assays <- ain[!ain %in% assays_se]
    existing_assays <- ain[ain %in% assays_se]
    if (length(not_existing_assays) > 0) {
      # some methods of ain are not in SummarizedExperiment object
      # check if there are still some methods of ain in the SummarizedExperiment object
      if (length(existing_assays) > 0) {
        warning(
          paste0(
            paste0(not_existing_assays, collapse = ", "),
            " not in SummarizedExperiment object. Check with names(assays(se)) which assays can be used for visualization!"
          )
        )
        message(paste0(
          paste0(existing_assays, collapse = ", "),
          " used for visualization."
        ))
        return(existing_assays)
      } else {
        stop(
          paste0(
            paste0(not_existing_assays, collapse = ", "),
            " not in SummarizedExperiment object. Check with names(assays(se)) which assays can be used for visualization!"
          )
        )
        return(NULL)
      }
    } else {
      return(ain)
    }
  }
}


#' Check parameters for DE analysis
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param comparisons Vector of comparisons that are performed in the DE analysis (from specify_comparisons method)
#' @param DE_method String specifying which DE method should be applied (limma, ROTS)
#' @param covariate String specifying which column to include as covariate into limma
#' @param logFC Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)
#' @param logFC_up Upper log2 fold change threshold (dividing into up regulated)
#' @param logFC_down Lower log2 fold change threshold (dividing into down regulated)
#' @param p_adj Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)
#' @param alpha Threshold for adjusted p-values or p-values
#' @param p_adj_method String specifying the method for adjusted p-values
#' @param B Number of bootstrapping for ROTS
#' @param K Number of top-ranked features for reproducibility optimization
#'
#' @return list of checked assays and condition column name
#'
#' @keywords internal
check_DE_parameters <- function(se, ain = NULL, condition = NULL, comparisons = NULL, DE_method = "limma", covariate = NULL, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, p_adj_method = "BH", alpha = 0.05, B = 100, K = 500){
  # check input parameters
  stopifnot(DE_method %in% c("limma", "ROTS"))
  stopifnot(p_adj_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")) # TODO: check which values available
  stopifnot(methods::is(logFC_up, "numeric"))
  stopifnot(methods::is(logFC_down, "numeric"))
  stopifnot(methods::is(B, "numeric"))
  stopifnot(methods::is(K, "numeric"))
  stopifnot(methods::is(alpha, "numeric"))
  stopifnot(methods::is(p_adj, "logical"))
  stopifnot(methods::is(logFC, "logical"))

  # check covariate
  if(!is.null(covariate)){
    if(!covariate %in% colnames(data.table::as.data.table(SummarizedExperiment::colData(se)))){
      stop(paste0("No column named ", covariate, " in SummarizedExperiment!"))
    }
  }

  # check comparisons
  # TODO

  # get condition
  condition <- get_condition_value(se, condition)

  # check input
  ain <- check_input_assays(se, ain)
  if(is.null(ain)){
    return(NULL)
  }
  return(list("ain" = ain, "condition" = condition))
}

#' Helper function to check the parameters for plotting the DE results
#'
#' @param de_res data table resulting of run_DE
#' @param ain String of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#'
#' @return list of valid inputs for plotting functions
#'
#' @keywords internal
check_plot_DE_parameters <- function(de_res, ain, comparisons){
  # check if logFC, P.Value, Change, Comparison, and Assay in de_res
  stopifnot("logFC" %in% colnames(de_res))
  stopifnot("Change" %in% colnames(de_res))
  stopifnot("P.Value" %in% colnames(de_res))
  stopifnot("Comparison" %in% colnames(de_res))
  stopifnot("Assay" %in% colnames(de_res))

  # check comparisons
  if(is.null(comparisons)){
    comparisons <- unique(de_res$Comparison)
    message("All comparisons of de_res will be visualized.")
  }
  valid_comparisons <- comparisons[comparisons %in% unique(de_res$Comparison)]
  not_valid_comparisons <- comparisons[!comparisons %in% unique(de_res$Comparison)]
  if(length(not_valid_comparisons) > 0){
    if(length(not_valid_comparisons) == length(c(comparisons))){
      # no valid comparisons
      stop("No valid comparison! Please change the comparisons parameter. Check with unique(de_res$Comparisons) which comparisons you can visualize.")
    } else {
      # some invalid comparisons --> notification
      warning(paste0(paste0(not_valid_comparisons, collapse = ", "), ": not valid comparisons. Only valid comparisons will be visualized."))
    }
  }
  comparisons <- valid_comparisons

  # check ain
  if(is.null(ain)){
    ain <- unique(de_res$Assay)
    message("All normalization methods of de_res will be visualized.")
  }
  valid_ain <- ain[ain %in% unique(de_res$Assay)]
  not_valid_ain <- ain[!ain %in% unique(de_res$Assay)]
  if(length(not_valid_ain) > 0){
    if(length(not_valid_ain) == length(c(ain))){
      # no valid input
      stop("No valid normalization methods! Please change the ain parameter. Check with unique(de_res$Assays) which normalization methods can be visualized.")
    } else {
      # some invalid inputs --> notification
      warning(paste0(paste0(not_valid_ain, collapse = ", "), ": not valid normalization methods. Only valid normalization methods will be visualized."))
    }
  }
  ain <- valid_ain
  return(list(de_res, ain, comparisons))
}

#' Helper function to check the parameters for plotting the DE stats of spike-in data sets
#'
#' @param stats data table resulting of get_spiked_stats_DE
#' @param ain String of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#'
#' @return list of valid inputs for plotting functions
#'
#' @keywords internal
check_stats_spiked_DE_parameters <- function(stats, ain, comparisons){
  # check Comparison, and Assay in stats
  stopifnot("Assay" %in% colnames(stats))
  stopifnot("Comparison" %in% colnames(stats))

  # check comparisons
  if(is.null(comparisons)){
    comparisons <- unique(stats$Comparison)
    message("All comparisons of stats will be visualized.")
  }
  valid_comparisons <- comparisons[comparisons %in% unique(stats$Comparison)]
  not_valid_comparisons <- comparisons[!comparisons %in% unique(stats$Comparison)]
  if(length(not_valid_comparisons) > 0){
    if(length(not_valid_comparisons) == length(c(comparisons))){
      # no valid comparisons
      stop("No valid comparison! Please change the comparisons parameter. Check with unique(stats$Comparisons) which comparisons you can visualize.")
    } else {
      # some invalid comparisons --> notification
      warning(paste0(paste0(not_valid_comparisons, collapse = ", "), ": not valid comparisons. Only valid comparisons will be visualized."))
    }
  }
  comparisons <- valid_comparisons

  # check ain
  if(is.null(ain)){
    ain <- unique(stats$Assay)
    message("All normalization methods of de_res will be visualized.")
  }
  valid_ain <- ain[ain %in% unique(stats$Assay)]
  not_valid_ain <- ain[!ain %in% unique(stats$Assay)]
  if(length(not_valid_ain) > 0){
    if(length(not_valid_ain) == length(c(ain))){
      # no valid input
      stop("No valid normalization methods! Please change the ain parameter. Check with unique(de_res$Assays) which normalization methods can be visualized.")
    } else {
      # some invalid inputs --> notification
      warning(paste0(paste0(not_valid_ain, collapse = ", "), ": not valid normalization methods. Only valid normalization methods will be visualized."))
    }
  }
  ain <- valid_ain
  return(list(stats, ain, comparisons))
}

