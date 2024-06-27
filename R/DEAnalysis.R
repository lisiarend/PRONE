
#' Create vector of comparisons for DE analysis (either by single condition (sep = NULL) or by combined condition)
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param condition Column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param sep Separator that separates both groups in the condition vector (NULL if condition composed only of single group)
#' @param control String of control samples (how the control condition is named) (NULL if no control sample)
#'
#' @return Vector of comparisons for DE analysis
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' comparisons <- specify_comparisons(tuberculosis_TMT_se, condition = NULL,
#'                             sep = NULL, control = NULL)
#'
specify_comparisons <- function(se, condition = NULL, sep = NULL, control = NULL){
  # get condition
  condition <- get_condition_value(se, condition)

  # get condition vector
  condition_vector <- unique(SummarizedExperiment::colData(se)[[condition]])

  comparisons <- c()
  for(index_a in seq_len(length(condition_vector)-1)){
    for(index_b in seq(index_a+1, length(condition_vector))){
      sample_a <- condition_vector[index_a]
      sample_b <- condition_vector[index_b]

      if(is.null(sep)){
        # only one condition to compare
        comparisons <- c(comparisons, paste0(sample_a, "-", sample_b))
      } else {
        # check if a or b are controls
        if(!is.null(control)){
          if((sample_a == control) || (sample_b == control)){
            next
          }
        }

        # check of how many individual groups the condition is composed
        splits_a <- stringr::str_split(sample_a, sep)[[1]]
        splits_b <- stringr::str_split(sample_b, sep)[[1]]
        must_match <- length(splits_a) - 1 # one group should stay variable

        # check how many groups are the same --> n - 1 groups need to be the same to add the comparison
        match_count <- 0
        for(i in seq_len(length(splits_a))){
          split_a <- splits_a[i]
          split_b <- splits_b[i]
          if(split_a == split_b){
            match_count <- match_count + 1
          }
        }
        if(match_count == must_match){
          comparisons <- c(comparisons, paste0(sample_b, "-", sample_a))
        }
      }
    }
  }

  # add control comparisons
  if(!is.null(control)){
    for(index_a in seq(length(condition_vector))){
      sample_a <- condition_vector[index_a]
      if(sample_a != control){
        comparisons <- c(comparisons, paste0(sample_a, "-", control))
      }
    }
  }

  return(as.factor(comparisons))
}

#' Fitting a linear model using limma
#'
#' @param data Data table of intensities (rows = proteins, cols = samples)
#' @param condition_vector Vector of experimental design specifying the condition(s) to compare
#' @param comparisons Vector of comparisons that are performed in the DE analysis (from specify_comparisons method)
#' @param covariate String specifying which column to include as covariate into limma
#' @param trend logical, should an intensity-dependent trend be allowed for the prior variance? If FALSE then the prior variance is constant. Alternatively, trend can be a row-wise numeric vector, which will be used as the covariate for the prior variance.
#' @param robust logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?
#'
#' @return eBayes object
#'
perform_limma <- function(data, condition_vector, comparisons, covariate = NULL, trend = TRUE, robust = TRUE){
  # Create design matrix
  groupsM <- as.factor(condition_vector)
  if(is.null(covariate)){
    designM <- stats::model.matrix(~0+groupsM)
    colnames(designM) <- levels(groupsM)
  } else {
    designM <- stats::model.matrix(~0+groupsM+covariate)
    colnames(designM)[seq(length(levels(groupsM)))] <- levels(groupsM)
  }

  # Fit linear moel
  fit <- limma::lmFit(data, designM)

  # Create contrasts
  contr <- limma::makeContrasts(contrasts = comparisons, levels = colnames(stats::coef(fit)))

  # Contrast fit and ebayes
  fit2 <- limma::contrasts.fit(fit, contr)
  ebfit <- limma::eBayes(fit2, trend=trend, robust = robust)
  return(ebfit)
}


#' Extract the DE results from eBayes fit of perform_limma function.
#'
#' @param fit eBayes object resulting from perform_limma method
#' @param comparisons Vector of comparisons that are performed in the DE analysis (from specify_comparisons method)
#' @param logFC Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)
#' @param logFC_up Upper log2 fold change threshold (dividing into up regulated)
#' @param logFC_down Lower log2 fold change threshold (dividing into down regulated)
#' @param p_adj Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)
#' @param alpha Threshold for adjusted p-values or p-values
#'
#' @return Data table with limma DE results
#'
extract_limma_DE <- function(fit, comparisons, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05){
  de_res <- NULL
  for (i in seq_len(length(comparisons))){
    # get results for each comparison
    top.table <- limma::topTable(fit, sort.by = "P", number=Inf, coef=c(i), adjust.method = "BH")
    gene_reg <- data.table::setDT(top.table, keep.rownames = "ID") # save row names as column
    if(logFC){
      # logFC
      if(p_adj){
        # p.adjust
        gene_reg$Change <- ifelse(gene_reg$logFC >= logFC_up & gene_reg$adj.P.Val <= alpha, "Up Regulated", ifelse(gene_reg$logFC <= logFC_down & gene_reg$adj.P.Val <= alpha, "Down Regulated", "No Change"))
      } else {
        # p.value
        gene_reg$Change <- ifelse(gene_reg$logFC >= logFC_up & gene_reg$P.Value <= alpha, "Up Regulated", ifelse(gene_reg$logFC <= logFC_down & gene_reg$P.Value <= alpha, "Down Regulated", "No Change"))
      }
    } else {
      # no logFC
      if(p_adj){
        # p.adjust
        gene_reg$Change <- ifelse(gene_reg$adj.P.Val <= alpha, "Significant Change", "No Change")
      } else {
        # p.value
        gene_reg$Change <- ifelse(gene_reg$P.Value <= alpha, "Significant Change", "No Change")
      }
    }
    gene_reg$Change <- as.factor(gene_reg$Change)
    gene_reg$Comparison <- rep(comparisons[i], nrow(gene_reg))
    if(is.null(de_res)){
      de_res <- gene_reg
    } else {
      de_res <- rbind(de_res, gene_reg)
    }
  }
  return(de_res)
}


#' Performing ROTS
#'
#' @param data Data table of intensities (rows = proteins, cols = samples)
#' @param condition Vector of experimental design specifying the condition(s) to compare
#' @param comparisons Vector of comparisons that are performed in the DE analysis (from specify_comparisons method)
#' @param condition_name String of name of condition in colData
#' @param coldata colData of the SummarizedExperiment
#' @param logFC Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)
#' @param logFC_up Upper log2 fold change threshold (dividing into up regulated)
#' @param logFC_down Lower log2 fold change threshold (dividing into down regulated)
#' @param p_adj Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)
#' @param alpha Threshold for adjusted p-values or p-values
#' @param B Number of bootstrapping for ROTS
#' @param K Number of top-ranked features for reproducibility optimization
#'
#' @return Data table with DE results
#'
perform_ROTS <- function(data, condition, comparisons, condition_name, coldata, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05, B =100, K = 500){
  de_res <- NULL
  for (comp in comparisons){
    # extract data of comparison
    sampleA <- strsplit(comp, "-")[[1]][1]
    sampleB <- strsplit(comp, "-")[[1]][2]
    coldata_chunk <- coldata[coldata[[condition_name]] %in% c(sampleA, sampleB),]
    coldata_chunk <- coldata_chunk %>% dplyr::arrange(factor(get(condition_name), levels = c(sampleA, sampleB)))
    dt_chunk <- data[, coldata_chunk$Column]
    # specify group vector
    group <- factor(coldata_chunk[[condition_name]], levels = c(sampleA, sampleB))
    group <- as.numeric(group)
    # remove protein groups with less than 2 valid values per group
    for (i in unique(group)){
      to_remove <- rowSums(is.na(dt_chunk[, which(group == i)])) >= 2
      dt_chunk <- dt_chunk[!to_remove,]
    }
    # perform ROTS
    rots <- ROTS::ROTS(data = dt_chunk, groups = group, B = B, K = K, seed = 1234)
    # extract ROTS results
    gene_reg <- data.table::data.table(ID = rownames(rots$data), logFC = rots$logfc, adj.P.Val = rots$FDR, P.Value = rots$pvalue)
    gene_reg$Comparison <- comp

    if(logFC){
      # logFC
      if(p_adj){
        # p.adjust
        gene_reg$Change <- ifelse(gene_reg$logFC >= logFC_up & gene_reg$adj.P.Val <= alpha, "Up Regulated", ifelse(gene_reg$logFC <= logFC_down & gene_reg$adj.P.Val <= alpha, "Down Regulated", "No Change"))
      } else {
        # p.value
        gene_reg$Change <- ifelse(gene_reg$logFC >= logFC_up & gene_reg$P.Value <= alpha, "Up Regulated", ifelse(gene_reg$logFC <= logFC_down & gene_reg$P.Value <= alpha, "Down Regulated", "No Change"))
      }
    } else {
      # no logFC
      if(p_adj){
        # p.adjust
        gene_reg$Change <- ifelse(gene_reg$adj.P.Val <= alpha, "Significant Change", "No Change")
      } else {
        # p.value
        gene_reg$Change <- ifelse(gene_reg$P.Value <= alpha, "Significant Change", "No Change")
      }
    }
    gene_reg$Change <- as.factor(gene_reg$Change)
    #gene_reg$Comparison <- rep(comparisons[i], nrow(gene_reg))
    if(is.null(de_res)){
      de_res <- gene_reg
    } else {
      de_res <- rbind(de_res, gene_reg)
    }
  }
  return(de_res)
}

#' Run DE analysis on a single normalized data set
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param method String specifying which assay should be used as input
#' @param comparisons Vector of comparisons that are performed in the DE analysis (from specify_comparisons method)
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param DE_method String specifying which DE method should be applied (limma, ROTS)
#' @param covariate String specifying which column to include as covariate into limma
#' @param logFC Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)
#' @param logFC_up Upper log2 fold change threshold (dividing into up regulated)
#' @param logFC_down Lower log2 fold change threshold (dividing into down regulated)
#' @param p_adj Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)
#' @param alpha Threshold for adjusted p-values or p-values
#' @param B Number of bootstrapping for ROTS
#' @param K Number of top-ranked features for reproducibility optimization
#' @param trend logical, should an intensity-dependent trend be allowed for the prior variance? If FALSE then the prior variance is constant. Alternatively, trend can be a row-wise numeric vector, which will be used as the covariate for the prior variance.
#' @param robust logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?
#'
#' @return Data table of DE results
#'
run_DE_single <- function(se, method, comparisons, condition = NULL, DE_method = "limma", covariate = NULL, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05, B = 100, K = 500, trend = TRUE, robust = TRUE){
  # check parameters (if users are calling this function instead of run_DE)
  params <- check_DE_parameters(se, ain = method, condition = condition, comparisons = comparisons, DE_method = DE_method, covariate = covariate, logFC = logFC, logFC_up = logFC_up, logFC_down = logFC_down, p_adj = p_adj, alpha = alpha, B = B, K = K)
  method <- params[["ain"]]
  condition <- params[["condition"]]

  # get covariate
  if(!is.null(covariate)){
    covariate <- SummarizedExperiment::colData(se)[[covariate]]
  }

  # prepare data
  dt <- as.data.frame(SummarizedExperiment::assays(se)[[method]])
  rownames(dt) <- data.table::as.data.table(SummarizedExperiment::rowData(se))$Protein.IDs
  dt <- dt[rowSums(is.na(dt)) != ncol(dt), ] # remove complete NAs
  condition_vector <- SummarizedExperiment::colData(se)[[condition]]

  # run DE
  if(DE_method == "limma"){
    # run limma
    fit <- perform_limma(data = dt, condition_vector = condition_vector, comparisons = comparisons, covariate = covariate, trend = trend, robust = robust)
    de_chunk <- extract_limma_DE(fit = fit, comparisons = comparisons, logFC = logFC, logFC_up = logFC_up, logFC_down = logFC_down, p_adj = p_adj, alpha = alpha)
  } else if (DE_method == "ROTS"){
    # run ROTS
    de_chunk <- perform_ROTS(data = dt, condition = condition, comparisons = comparisons, condition_name = condition, coldata = data.table::as.data.table(SummarizedExperiment::colData(se)), logFC = logFC, logFC_up = logFC_up, logFC_down = logFC_down, p_adj = p_adj, alpha = alpha, B = B, K = K)
  }

  # add other information
  de_chunk$Assay <- method
  de_chunk$Protein.IDs <- de_chunk$ID
  de_chunk$ID <- NULL
  de_chunk <- de_chunk[, c("Protein.IDs", "logFC", "P.Value", "adj.P.Val", "Change", "Comparison", "Assay")]

  # add missing proteins (that could not be calculated)
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  for(comp in comparisons){
    # take comparison
    de_res_comp <- de_chunk[de_chunk$Comparison == comp,]
    if(sum(rd$Protein.IDs %in% de_res_comp$Protein.IDs) != nrow(rd)){
      missing <- rd$Protein.IDs[! rd$Protein.IDs %in% de_res_comp$Protein.IDs]
      if(length(missing) > 0){
        de <- data.table::data.table("Protein.IDs" = missing, "logFC" = rep(NA, length(missing)), "P.Value" = rep(NA, length(missing)), "adj.P.Val" = rep(NA, length(missing)), "Change" = rep("No Change", length(missing)), "Comparison" = rep(comp, length(missing)), "Assay" = rep(method, length(missing)))
        de_chunk <- rbind(de_chunk, de)
      }
    }
  }
  de_chunk$Change[is.na(de_chunk$Change)] <- "No Change"

  if(is.null(S4Vectors::metadata(se)$spike_column)){
    rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))[,c("Protein.IDs", "Gene.Names", "IDs")]
  } else {
    # add spike info
    spike_column <- S4Vectors::metadata(se)$spike_column
    cols <- c("Protein.IDs", "Gene.Names", "IDs", spike_column)
    rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))[, cols]
  }

  de_chunk <- merge(de_chunk, rd, by="Protein.IDs")
  de_chunk$Comparison <- factor(de_chunk$Comparison, levels = comparisons)
  return(de_chunk)
}

#' Run DE analysis of a selection of normalized data sets
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param comparisons Vector of comparisons that are performed in the DE analysis (from specify_comparisons method)
#' @param ain Vector of strings which assay should be used as input (default NULL).
#'            If NULL then all normalization of the se object are plotted next to each other.
#' @param condition column name of condition (if NULL, condition saved in SummarizedExperiment will be taken)
#' @param DE_method String specifying which DE method should be applied (limma, ROTS)
#' @param covariate String specifying which column to include as covariate into limma
#' @param logFC Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)
#' @param logFC_up Upper log2 fold change threshold (dividing into up regulated)
#' @param logFC_down Lower log2 fold change threshold (dividing into down regulated)
#' @param p_adj Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)
#' @param alpha Threshold for adjusted p-values or p-values
#' @param B Number of bootstrapping for ROTS
#' @param K Number of top-ranked features for reproducibility optimization
#' @param trend logical, should an intensity-dependent trend be allowed for the prior variance? If FALSE then the prior variance is constant. Alternatively, trend can be a row-wise numeric vector, which will be used as the covariate for the prior variance.
#' @param robust logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?
#'
#' @return Data table of DE results of selected normalized data sets
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' comparisons <- specify_comparisons(tuberculosis_TMT_se, condition = NULL,
#'                             sep = NULL, control = NULL)
#' de_res <- run_DE(tuberculosis_TMT_se, comparisons,
#'               ain = NULL, condition = NULL, DE_method = "limma",
#'               logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE,
#'               alpha = 0.05, B = 100, K = 500, trend = TRUE, robust = TRUE)
#'
run_DE <- function(se, comparisons, ain = NULL, condition = NULL, DE_method = "limma", covariate = NULL, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05, B = 100, K = 500, trend = TRUE, robust = TRUE){
  # check parameters
  params <- check_DE_parameters(se, ain = ain, condition = condition, comparisons = comparisons, DE_method = DE_method, covariate = covariate, logFC = logFC, logFC_up = logFC_up, logFC_down = logFC_down, p_adj = p_adj, alpha = alpha, B = B, K = K)
  ain <- params[["ain"]]
  condition <- params[["condition"]]

  if("raw" %in% ain){
    message("DE Analysis will not be performed on raw data.")
    ain <- ain[ain != "raw"]
  }

  # run DE
  de_res <- NULL
  for(method in c(ain)){
    de_chunk <- run_DE_single(se, method = method, condition = condition, comparisons = comparisons, DE_method = DE_method, covariate = covariate, logFC = logFC, logFC_up = logFC_up, logFC_down = logFC_down, p_adj = p_adj, alpha = alpha, B = B, K = K, trend = trend, robust = robust)
    # add to overall results
    if(is.null(de_res)){
      de_res <- de_chunk
    } else {
      de_res <- rbind(de_res, de_chunk)
    }
    message(paste0(method, ": DE analysis completed."))
  }
  return(de_res)
}


#' Apply other thresholds to DE results
#'
#' @param de_res data table resulting of run_DE
#' @param logFC Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)
#' @param logFC_up Upper log2 fold change threshold (dividing into up regulated)
#' @param logFC_down Lower log2 fold change threshold (dividing into down regulated)
#' @param p_adj Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)
#' @param alpha Threshold for adjusted p-values or p-values
#'
#' @return data table updating the Change column with the newly applied thresholds
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' de_res <- apply_thresholds(tuberculosis_TMT_de_res, logFC = FALSE,
#'                            p_adj = TRUE, alpha = 0.01)
#'
apply_thresholds <- function(de_res, logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05){
  stopifnot("Change" %in% colnames(de_res))

  if(logFC){
    # logFC
    if(p_adj){
      # p.adjust
      de_res$Change <- ifelse(de_res$logFC >= logFC_up & de_res$adj.P.Val <= alpha, "Up Regulated", ifelse(de_res$logFC <= logFC_down & de_res$adj.P.Val <= alpha, "Down Regulated", "No Change"))
    } else {
      # p.value
      de_res$Change <- ifelse(de_res$logFC >= logFC_up & de_res$P.Value <= alpha, "Up Regulated", ifelse(de_res$logFC <= logFC_down & de_res$P.Value <= alpha, "Down Regulated", "No Change"))
    }
  } else {
    # no logFC
    if(p_adj){
      # p.adjust
      de_res$Change <- ifelse(de_res$adj.P.Val <= alpha, "Significant Change", "No Change")
    } else {
      # p.value
      de_res$Change <- ifelse(de_res$P.Value <= alpha, "Significant Change", "No Change")
    }
  }
  de_res$Change <- as.factor(de_res$Change)
  return(de_res)
}


#' Get overview table of DE results
#'
#' @param de_res data table resulting of run_DE
#'
#' @return data table of numbers of DE proteins per comparison and per normalization method
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_de_res)
#' get_overview_DE(tuberculosis_TMT_de_res)
#'
get_overview_DE <- function(de_res){
  if("Significant Change" %in% de_res$Change){
    dt <- de_res %>% dplyr::group_by(Assay, Comparison) %>% dplyr::summarize('Significant Change' = sum(Change == "Significant Change", na.rm=TRUE)) %>% data.table::as.data.table()
  } else if ("Up Regulated" %in% de_res$Change | "Down Regulated" %in% de_res$Change){
    dt <- de_res %>% dplyr::group_by(Assay, Comparison) %>% dplyr::summarize('Up Regulated' = sum(Change == "Up Regulated", na.rm=TRUE),
                                                                                  'Down Regulated' = sum(Change == "Down Regulated", na.rm=TRUE)) %>% data.table::as.data.table()
  } else {
    stop("No DE proteins found for any comparison and normalization method!")
  }
  dt$Assay <- factor(dt$Assay, levels = unique(de_res$Assay))
  return(dt)
}
