#'
#' @rdname tuberculosis_TMT_se
#' @name tuberculosis_TMT_se
#' @title Example SummarizedExperiment of a real-world proteomics data set
#' @description
#' A \code{SummarizedExperiment} containing the raw and log2-scaled data of 301 proteins measured in 20 samples
#' @source Biadglegne et al. Mycobacterium tuberculosis Affects Protein and Lipid Content of Circulating Exosomes in Infected Patients Depending on Tuberculosis Disease State. Biomedicines 10.4 (Mar. 2022), p. 783. doi: 10.3390/ biomedicines10040783.
#' @usage data(tuberculosis_TMT_se)
"tuberculosis_TMT_se"

#'
#' @rdname tuberculosis_TMT_de_res
#' @name tuberculosis_TMT_de_res
#' @title Example data.table of DE results of a real-world proteomics data set
#' @description
#' A \code{data.table} containing the DE results of the tuberculosis_TMT_se data set (limma, logFC > 1, logFC < -1, p.adj < 0.05)
#' @source Biadglegne et al. Mycobacterium tuberculosis Affects Protein and Lipid Content of Circulating Exosomes in Infected Patients Depending on Tuberculosis Disease State. Biomedicines 10.4 (Mar. 2022), p. 783. doi: 10.3390/ biomedicines10040783.
#' @usage data(tuberculosis_TMT_de_res)
"tuberculosis_TMT_de_res"

#'
#' @rdname spike_in_se
#' @name spike_in_se
#' @title Example SummarizedExperiment of a spike-in proteomics data set
#' @description
#' A \code{SummarizedExperiment} containing the raw and log2-scaled data of 301 proteins measured in 20 samples. Due to size restriction, we only included the relevant columns of the original proteinGroups.txt of MaxQuant. 
#' @source Jürgen Cox, Marco Y. Hein, Christian A. Luber, Igor Paron, Nagarjuna Nagaraj, and Matthias Mann.Accurate Proteome-wide Label-free Quantification by Delayed Normalization and Maximal Peptide Ratio Extraction, Termed MaxLFQ. Molecular & Cellular Proteomics 13.9 (Sept. 2014), pp. 2513–2526. <https://doi.org/10.1074/mcp.M113.031591>.
#' @usage data(spike_in_se)
"spike_in_se"

#'
#' @rdname spike_in_de_res
#' @name spike_in_de_res
#' @title Example data.table of DE results of a spike-in proteomics data set
#' @description
#' A \code{data.table} containing the DE results of the spike_in_se data set (limma, logFC > 1, logFC < -1, p.adj < 0.05)
#' @source Jürgen Cox, Marco Y. Hein, Christian A. Luber, Igor Paron, Nagarjuna Nagaraj, and Matthias Mann.Accurate Proteome-wide Label-free Quantification by Delayed Normalization and Maximal Peptide Ratio Extraction, Termed MaxLFQ. Molecular & Cellular Proteomics 13.9 (Sept. 2014), pp. 2513–2526. <https://doi.org/10.1074/mcp.M113.031591>.
#' @usage data(spike_in_de_res)
"spike_in_de_res"




