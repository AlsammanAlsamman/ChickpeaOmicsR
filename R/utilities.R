#' Match gene expression data with metadata
#' @param dataList A list containing metadata and gene expression data (counts) or individual metadata and counts as separate arguments
#' @param counts Optional argument. Gene expression data (counts), required if metadata is provided separately from dataList.
#' @param metadata Optional argument. Metadata, required if counts is provided separately from dataList.
#' @export
#' @return A list containing the matched counts and metadata
align_gene_expression_metadata<-function(dataList = NULL, metadata = NULL, counts = NULL){
  # Create a list to store the gene expression data and metadata
  gene_meta <- list()
  # If dataList is provided, use it
  if (!is.null(dataList)) {
    if (length(dataList) == 0) {
      stop("Please provide a list of gene expression data")
    }
    counts <- dataList$counts
    metadata <- dataList$metadata
  }
  
  # If metadata and counts are provided separately
  if (!is.null(metadata) && !is.null(counts)) {
    if (length(metadata) == 0) {
      stop("Please provide metadata")
    }
    if (length(counts) == 0) {
      stop("Please provide counts")
    }
  }
  
  # Ensure that both metadata and counts are available for matching
  if (is.null(metadata) || is.null(counts)) {
    stop("Both metadata and counts must be provided.")
  }
  
  # Remove samples not in metadata
  counts <- counts[, colnames(counts) %in% rownames(metadata)]
  
  # Remove samples not in data
  metadata <- metadata[rownames(metadata) %in% colnames(counts),]
  
  gene_meta$counts <- counts
  gene_meta$metadata <- metadata
  
  return(gene_meta)
}

#' Get report for data
#' @param data The data to get the report for
#' @return A report for the data
#' @export
#' @examples
#' get_report(data)
#' get_report(get_study_data("Salinity"))
get_data_report <- function(data) {
  # create a list to store the report
  report <- list()
  # get the number of samples
  report$samples <- ncol(data$counts)
  # get the number of genes
  report$genes <- nrow(data$counts)
  # get the number of studies
  report$studies <- length(unique(data$metadata$study))
  # get the number of experiments
  report$experiments <- length(unique(data$metadata$experiment))
  # print the report
  print("These data have been generated from the ChickpeaOmicsR package.")
  print(paste("The data has", report$samples, "samples and", report$genes, "genes."))
  print(paste("The data is from", report$studies, "studies", "and", report$experiments, "experiments."))
  return(report)
}