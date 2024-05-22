

#' Intersect top N enrichment terms per normalization method
#'
#' @param se SummarizedExperiment containing all necessary information of the proteomics data set
#' @param de_res data table resulting of run_DE
#' @param ain Vector of strings of normalization methods to visualize (must be valid normalization methods saved in de_res)
#' @param comparisons Vector of comparisons (must be valid comparisons saved in de_res)
#' @param id_column String specifying the column of the rowData of the SummarizedExperiment object which includes the gene names
#' @param organism Organism name (gprofiler parameter)
#' @param per_comparison Boolean specifying whether the enrichment analysis should be performed per comparison (TRUE) or on all given comparisons together (FALSE)
#' @param sources Vector of data sources to use (gprofiler parameter)
#' @param top Number of enrichment terms to extract for each normalization method
#'
#' @return list of ggplot objects or single ggplot object
#' @export
#'
#' @examples
#' data(tuberculosis_TMT_se)
#' data(tuberculosis_TMT_de_res)
#' plot_intersection_enrichment(tuberculosis_TMT_se, tuberculosis_TMT_de_res,
#'                 ain = c("IRS_on_RobNorm", "IRS_on_Median"),
#'                 comparisons = NULL, id_column = "Gene.Names",
#'                 organism = "hsapiens", per_comparison = TRUE,
#'                 sources = c("GO:BP", "GO:MF", "GO:CC"), top = 10)
#'
plot_intersection_enrichment <- function(se, de_res, ain = NULL, comparisons = NULL, id_column = "Gene.Names", organism = "hsapiens", per_comparison = TRUE, sources = c("GO:BP", "GO:MF", "GO:CC"), top = 10){
  # check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, comparisons)
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparisons <- tmp[[3]]
  # check top10
  stopifnot(is.numeric(top))

  # check id_column
  rd <- data.table::as.data.table(SummarizedExperiment::rowData(se))
  if(!id_column %in% colnames(rd)){
    # if id_column not in rowData
    stop(paste0(id_column, " not in rowData of the SummarizedExperiment object!"))
  }

  de_res <- de_res[de_res$Change != "No Change",]
  de_res <- de_res[de_res$Assay %in% ain,]
  de_res <- de_res[de_res$Comparison %in% comparisons,]

  if(per_comparison){
    plots <- list()
    for(comp in comparisons){
      dt <- de_res[de_res$Comparison == comp,]
      # extract queries
      queries <- lapply(ain, function(method){
        dt <- dt[dt$Assay == method,]
        query <- dt[[id_column]]
        query <- unique(query[query != ""])
        query
      })
      names(queries) <- ain

      # gprofiler
      gres <- gprofiler2::gost(queries, organism = organism, sources = sources)
      gres <- gres$result

      # extract top
      top_go <- NULL
      for(method in ain){
        dt <- gres[gres$query == method,]
        dt <- dt %>% dplyr::arrange(p_value) %>% dplyr::group_by(source) %>% dplyr::slice(seq_len(top))
        dt <- dt[, c("query", "term_name", "source")]
        if(is.null(top_go)){
          top_go <- dt
        } else {
          top_go <- rbind(top_go, dt)
        }
      }
      top_go$present <- TRUE
      dt <- data.table::dcast(top_go, ...~query, value.var = "present")
      dt <- dt[order(dt$source),]
      terms <- paste0(dt$term_name, " (", dt$source, ")")
      source <- dt$source
      dt$term_name <- NULL
      dt$source <- NULL
      dt <- as.data.frame(dt)
      rownames(dt) <- terms
      dt[dt == TRUE] <- "Yes"
      dt[is.na(dt)] <- "No"
      dt$Term <- rownames(dt)
      melted_dt <- data.table::melt(data.table::as.data.table(dt), value.name = "Present", variable.name = "Assay", measure.vars = colnames(dt)[colnames(dt) != "Term"])
      # plot
      p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("Assay"), y = get("Term"), fill = get("Present"))) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_manual(name = "Significant", values = c("No" = "grey80", "Yes" = "#D55E00")) +
        ggplot2::labs(x = "Normalization Method", y = "Terms") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
      plots[[comp]] <- p
    }
    return(plots)
  } else {
    # extract queries
    queries <- lapply(ain, function(method){
      dt <- de_res[de_res$Assay == method,]
      query <- dt[[id_column]]
      # split by ";"
      query <- unlist(strsplit(query, ";"))
      query <- unique(query[query != ""])
      query
    })
    names(queries) <- ain

    # gprofiler
    gres <- gprofiler2::gost(queries, organism = organism, sources = sources)
    gres <- gres$result

    # extract top
    top_go <- NULL
    for(method in ain){
      dt <- gres[gres$query == method,]
      dt <- dt %>% dplyr::arrange(p_value) %>% dplyr::group_by(source) %>% dplyr::slice(seq_len(top))
      dt <- dt[, c("query", "term_name", "source")]
      if(is.null(top_go)){
        top_go <- dt
      } else {
        top_go <- rbind(top_go, dt)
      }
    }
    top_go$present <- TRUE
    dt <- data.table::dcast(data.table::as.data.table(top_go), ...~query, value.var = "present")
    dt <- dt[order(dt$source),]
    terms <- dt$term_name
    source <- dt$source
    dt$term_name <- NULL
    dt$source <- NULL
    dt <- as.data.frame(dt)
    rownames(dt) <- terms
    dt[dt == TRUE] <- "Yes"
    dt[is.na(dt)] <- "No"
    dt$Term <- rownames(dt)
    melted_dt <- data.table::melt(data.table::as.data.table(dt), value.name = "Present", variable.name = "Assay", measure.vars = colnames(dt)[colnames(dt) != "Term"])

    # plot
    p <- ggplot2::ggplot(melted_dt, ggplot2::aes(x = get("Assay"), y = get("Term"), fill = get("Present"))) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_manual(name = "Significant", values = c("grey80", "#D55E00")) +
      ggplot2::labs(x = "Normalization Method", y = "Terms") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
    return(p)
  }
}
