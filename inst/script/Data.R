## This file contains code for recreating the data files in the extdata folder. Details about the sources of the individual data files are provided in R/data.R and vignettes.

# Tuberculosis TMT data 
# @source Biadglegne et al. Mycobacterium tuberculosis Affects Protein and Lipid Content of Circulating Exosomes in Infected Patients Depending on Tuberculosis Disease State. Biomedicines 10.4 (Mar. 2022), p. 783. doi: 10.3390/ biomedicines10040783.

data_path <- readPRONE_example("tuberculosis_TMT_proteinGroups.txt")
md_path <- readPRONE_example("tuberculosis_TMT_metadata.txt")

data <- read.csv(data_path, sep ="\t")
md <- read.csv(md_path, sep = "\t")
md$Group <- plyr::mapvalues(md$Group, from=c("Common.reference","Healthy.control", "Pulmonary.tuberculosis", "Tuberculous.lymphadenitis", "Treated.pulmonary.tuberculosis"), to = c("ref", "HC", "PTB", "TBL", "Rx"))
md$Column <- stringr::str_replace_all(md$Column, " ", ".")

ref_samples <- md[md$Group == "ref",]$Column

tuberculosis_TMT_se <- load_data(data, md, protein_column = "Protein.IDs", gene_column = "Gene.names", ref_samples = ref_samples, batch_column = "Pool", condition_column = "Group", label_column = "Label")

# Spike-in Data 
# @source Jürgen Cox, Marco Y. Hein, Christian A. Luber, Igor Paron, Nagarjuna Nagaraj, and Matthias Mann.Accurate Proteome-wide Label-free Quantification by Delayed Normalization and Maximal Peptide Ratio Extraction, Termed MaxLFQ. Molecular & Cellular Proteomics 13.9 (Sept. 2014), pp. 2513–2526. <https://doi.org/10.1074/mcp.M113.031591>.

data_path <- readPRONE_example("spike_in_proteinGroups.txt")
md_path <- readPRONE_example("spike_in_metadata.txt")

data <- read.csv(data_path, sep = "\t")
md <- read.csv(md_path, sep = "\t")

# Check if some protein groups are mixed
mixed <- grepl("Homo sapiens.*Escherichia|Escherichia.*Homo sapiens", data$Fasta.headers)
data <- data[!mixed,]

data$Spiked <- rep("HUMAN", nrow(data))
data$Spiked[grepl("ECOLI", data$Fasta.headers)] <- "ECOLI"

spike_in_se <- load_spike_data(data, md, spike_column = "Spiked", spike_value = "ECOLI", spike_concentration = "Concentration",protein_column = "Protein.IDs", gene_column = "Gene.names", ref_samples = NULL, batch_column = NULL, condition_column = "Condition", label_column = "Label")

# Normalized Data
tuberculosis_TMT_se <- normalize_se(tuberculosis_TMT_se, methods = c("IRS_on_RobNorm", "IRS_on_Median"))
spike_in_se <- normalize_se(spike_in_se, methods = c("RobNorm", "Median", "LoessF", "VSN"))

# Differential Expression Results
se <- remove_reference_samples(tuberculosis_TMT_se)
comparisons <- specify_comparisons(se, condition = "Group", sep = "_", control = NULL)
tuberculosis_TMT_de_res <- run_DE(se = se,
                 comparisons = comparisons,
                 ain = NULL,
                 condition = NULL,
                 DE_method = "limma",
                 covariate = NULL,
                 logFC = TRUE,
                 logFC_up = 1,
                 logFC_down = -1,
                 p_adj = TRUE,
                 alpha = 0.05,
                 B = 100,
                 K = 500)

comparisons <- specify_comparisons(spike_in_se, condition = "Condition", sep = "_", control = NULL)
spike_in_de_res <- run_DE(se = spike_in_se,
                                  comparisons = comparisons,
                                  ain = NULL,
                                  condition = NULL,
                                  DE_method = "limma",
                                  covariate = NULL,
                                  logFC = TRUE,
                                  logFC_up = 1,
                                  logFC_down = -1,
                                  p_adj = TRUE,
                                  alpha = 0.05,
                                  B = 100,
                                  K = 500)


usethis::use_data(tuberculosis_TMT_se, overwrite = TRUE, compress = "xz")
usethis::use_data(tuberculosis_TMT_de_res, overwrite = TRUE, compress = "xz")
usethis::use_data(spike_in_se, overwrite = TRUE, compress = "xz")
usethis::use_data(spike_in_de_res, overwrite = TRUE, compress = "xz")


