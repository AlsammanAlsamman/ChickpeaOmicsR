#' Get the data for a specific stress study
#' @param stress The stress study to get the data for
#' @return A list with the counts and metadata for the study
#' @export
#' @examples
#' get_stress_data("Salinity")
get_study_data <- function(study = "all") {
  #if vector provided print error
  if (is.vector(study)) {
    stop("Please provide a single study name")
  }
  #study="Salinity"
  # create a list to store the data
  data <- list()
  # Load the data
  loadCaExpressionCounts()
  loadCaExpressionMetadata()
  if (study == "all") {
    data$counts <- CaExpressionCounts
    data$metadata <- CaExpressionMetadata
  } else {
    # samples of this study
    studysamples<-rownames(CaExpressionMetadata[CaExpressionMetadata$study == study, ])
    data$counts <- CaExpressionCounts[, colnames(CaExpressionCounts) %in% studysamples]
    data$metadata <- CaExpressionMetadata[CaExpressionMetadata$study == study, ]
  }
  return(data)
}

#' Get the data for a specific set of samples
#' @param samples The samples to get the data for
#' @return A list with the counts and metadata for the samples
#' @export
#' @examples
#' get_samples_data(c("SRR10331638", "SRR10331639", "SRR10331640", "SRR10331641", "SRR10331642", "SRR10331643", "SRR10331644", "SRR10331645")
get_samples_data <- function(samples = c()) {
  # create a list to store the data
  data <- list()
  # Load the data
  loadCaExpressionCounts()
  loadCaExpressionMetadata()
  if (length(samples) == 0) {
    data$counts <- CaExpressionCounts
    data$metadata <- CaExpressionMetadata
  } else {
    # subset the data
    data$counts <- CaExpressionCounts[, colnames(CaExpressionCounts) %in% samples]
    data$metadata <- CaExpressionMetadata[rownames(CaExpressionMetadata) %in% samples, ]
  }
  return(data)
}

#' Get Metadata for a specific samples, group, experiment, or study
#' @param samples The samples to get the metadata for
#' @param group The group to get the metadata for
#' @param experiment The experiment to get the metadata for
#' @param study The study to get the metadata for
#' @return Metadata for the samples, group, experiment, or study
#' @export
#' @examples
#' get_metadata(samples = c("SRR10331638", "SRR10331639", "SRR10331640", "SRR10331641", "SRR10331642", "SRR10331643", "SRR10331644", "SRR10331645"))

get_metadata <- function(samples = NULL, group = NULL, experiment = NULL, study = NULL) {
  # check that only one argument is provided
  if (sum(is.null(c(samples, group, experiment, study))) > 1) {
    stop("Please provide only one argument")
  }
  # Load the data
  loadCaExpressionMetadata()
  metadata <- CaExpressionMetadata
  # Filter the metadata based on the input
  if (!is.null(samples)) {
    metadata <- metadata[rownames(metadata) %in% samples, ]
  }
  if (!is.null(group)) {
    metadata <- metadata[metadata$group == group, ]
  }
  if (!is.null(experiment)) {
    metadata <- metadata[metadata$Experiment == experiment, ]
  }
  if (!is.null(study)) {
    metadata <- metadata[metadata$study == study, ]
  }
  return(metadata)
}

#' Get GWAS data for specific genes
#' @param geneids The genes to get the GWAS data for
#' @param pvalue The p-value threshold to use
#' @param phenotype The phenotype to get the GWAS data for
#' @export
#' @return GWAS data for the genes
get_gwas_data <- function(geneids, pvalue = 0.05, phenotype = NULL) {
  # Load the data
  loadCaGWAS()
  # Filter the data based on the input
  gwas_data <- CaGWAS[CaGWAS$geneid %in% geneids & CaGWAS$pvalue < pvalue, ]
  if (!is.null(phenotype)) {
    gwas_data <- gwas_data[gwas_data$phenotype == phenotype, ]
  }
  return(gwas_data)
}