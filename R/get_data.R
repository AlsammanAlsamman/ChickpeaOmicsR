#' @title Loading required data
#' @description ChickpeaOmicsR is a package that provides gene expression data for chickpea under different stress conditions.
#' @export
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ChickpeaOmicsR!")
  print("This package provides gene expression data for chickpea under different stress conditions.")
  print("Several Data MayBe Loaded to the Environment.")
}

#' Load the gene expression count data
#' @export
loadCaExpressionCounts <- function() {
 # check if the data is already loaded
  if (!exists("CaExpressionCounts")) {
    # create global variable of the data
    data(CaExpressionCounts)
    CaExpressionCounts<<-CaExpressionCounts
  }
}

#' Load CaExpressionMetadata
#' @export
#' @examples
#' loadCaExpressionMetadata()
loadCaExpressionMetadata <- function() {
  # check if the data is already loaded
  if (!exists("CaExpressionMetadata")) {
    # create global variable of the data
    data(CaExpressionMetadata)
    CaExpressionMetadata<<-CaExpressionMetadata
  }
}

#' Load CaGenomeAnnotation
#' @export
#' @examples
#' loadCaGenomeAnnotation()
loadCaGenomeAnnotation <- function() {
  # check if the data is already loaded
  if (!exists("CaGenomeAnnotation")) {
    # create global variable of the data
    data(CaGenomeAnnotation)
    CaGenomeAnnotation<<-CaGenomeAnnotation
  }
}

#data("CaProteinEnrichment")

#' Load CaProteinEnrichment
#' @export
#' @examples
#' loadCaProteinEnrichment()
loadCaProteinEnrichment <- function() {
  # check if the data is already loaded
  if (!exists("CaProteinEnrichment")) {
    # create global variable of the data
    data(CaProteinEnrichment)
    CaProteinEnrichment<<-CaProteinEnrichment
  }
}


#' Get the data for a specific stress study
#' @param stress The stress study to get the data for
#' @return A list with the counts and metadata for the study
#' @export
#' @examples
#' get_stress_data("Salinity")
get_study_data <- function(study = "all") {
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
    samples <- unique(CaExpressionMetadata[CaExpressionMetadata$study == study, ]$sampleid)
    data$counts <- CaExpressionCounts[, samples]
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
    data$counts <- CaExpressionCounts[, samples]
    data$metadata <- CaExpressionMetadata[CaExpressionMetadata$sampleid %in% samples, ]
  }
  return(data)
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
  # get the dimensions of the data
  report$dimensions <- dim(data$counts)
  # get the number of samples
  report$samples <- ncol(data$counts)
  # get the number of genes
  report$genes <- nrow(data$counts)
  # get the number of studies
  report$studies <- length(unique(data$metadata$study))
  # get the number of experiments
  report$experiments <- length(unique(data$metadata$Experiment))
  # print the report
  print("These data have been generated from the ChickpeaOmicsR package.")
  print(paste("The data has", report$samples, "samples and", report$genes, "genes."))
  print(paste("The data is from", report$studies, "studies", "and", report$experiments, "experiments."))
  return(report)
}

#' Get gene expression counts for a specific geneid
#' @param genenames The gene or genes to get the counts for
#' @return Expression Count matrix for the gene or genes
#' @export
#' @examples
#' get_gene_express_bygeneid(c("GRX3", "HIS2A","HSFA1B","HSFA2","HSFA3"))
get_gene_express_bygeneid <- function(geneids, gene_express_count=NULL){
  # Load the data
  geneExpressionCounts <- NULL
  if (is.null(gene_express_count)) {
    print("Gene expression count data is not provided. Using default data.")
    loadCaExpressionCounts()
    geneExpressionCounts <- CaExpressionCounts
  } else {
    geneExpressionCounts <- gene_express_count
  }
  # subset the data
  gene_express <- geneExpressionCounts[geneids, ]
  return(gene_express)
}

#' Get gene expression counts for a specific samples
#' @param samples The samples to get the counts for
#' @return Expression Count matrix for the samples
#' @export
#' @examples
#' get_gene_express_bysamples(c("SRR10331638", "SRR10331639", "SRR10331640", "SRR10331641", "SRR10331642", "SRR10331643", "SRR10331644", "SRR10331645"))
get_gene_express_bysamples<- function(samples){
  # Load the data
  loadCaExpressionCounts()
  # subset the data
  gene_express <- CaExpressionCounts[, samples]
  return(gene_express)
}


#' Get gene expression data using GO term id
#' @param go_id The GO term id to get the data for
#' @param gene_express_count The gene expression count
#' @return Gene expression data
#' @export
#' @examples
#' get_gene_express_bygoid("GO:0071704","GO:1901564","GO:0005622","GO:0005634","GO:0043226","GO:0043227","GO:0043229","GO:0043231","GO:0110165")
get_gene_express_bygoid <- function(go_id, gene_express_count=NULL){
  # Load the data
  loadCaProteinEnrichment()
  geneExpressionCounts <- NULL
  if (is.null(gene_express_count)) {
    print("Gene expression count data is not provided. Using default data.")
    loadCaExpressionCounts()
    geneExpressionCounts <- CaExpressionCounts
  } else {
    geneExpressionCounts <- gene_express_count
  }  
  # subset the data
  geneids <- CaProteinEnrichment[CaProteinEnrichment$termid %in% go_id, ]$geneid
  # remove NA values
  geneids <- geneids[!is.na(geneids)]
  geneids <- unique(geneids)
  gene_express <- geneExpressionCounts[geneids, ]
  return(gene_express)
}

#' Get gene expression data using GO term keywords
#' @param keywords The keywords to search for
#' @param gene_expression_count The gene expression count
#' @return Gene expression data
#' @export
#' @examples
#' get_gene_express_gotermkeywords(c("targeting"))
#' get_gene_express_gotermkeywords(c("targeting"), CaExpressionCounts)
get_gene_express_gokeywords<- function(keywords, gene_expression_count=NULL){
  #keywords<-"F-box"
  goids<-get_goids_bykewords(keywords)
  GeneExpressionCounts<-NULL
  if (is.null(gene_expression_count)) {
    print("Gene expression count data is not provided. Using default data.")
    loadCaExpressionCounts()
    GeneExpressionCounts<-CaExpressionCounts
  } else {
    GeneExpressionCounts<-gene_expression_count
  }
  geneids<-get_geneids_bygoids(goids)
  gene_express<-GeneExpressionCounts[geneids,]
  return(gene_express)
}


#' Get gene ids by GO ids
#' @param goids The GO ids to search for
#' @return Gene ids
#' @export
#' @examples
#' get_geneids_bygoids(c("GO:0071704","GO:1901564","GO:0005622","GO:0005634","GO:0043226","GO:0043227","GO:0043229","GO:0043231","GO:0110165"))
get_geneids_bygoids<- function(goids){
  loadCaProteinEnrichment()
  geneids<-CaProteinEnrichment[CaProteinEnrichment$termid %in% goids,]$geneid
  geneids<-unique(geneids)
  # remove NA values
  geneids<-geneids[!is.na(geneids)]
  return(geneids)
}


#' Get GO ids by keywords
#' @param keywords The keywords to search for
#' @return GO ids
#' @export
#' @import stringr
#' @examples
#' get_goids_bykewords(c("targeting"))
get_goids_bykewords <- function(keywords){
  loadCaProteinEnrichment()
  terms<-CaProteinEnrichment$term
  terms<-tolower(terms)
  keywords<-tolower(keywords)
  goids<-CaProteinEnrichment[grepl(paste(keywords, collapse="|"), terms),]$termid
  goids<-unique(goids)
  return(goids)
}


#' Get gene annotation data for a specific geneid
#' @return Gene annotation data
#' @param geneids The gene or genes to get the annotation for
#' @export
#' @examples
#' get_gene_annotation(c("GRX3", "HIS2A","HSFA1B","HSFA2","HSFA3"))
#' get_gene_annotation()
get_gene_annotation_bygeneid <- function(geneids = c()){
  # Load the data
  loadCaGenomeAnnotation()
  if (length(geneids) == 0) {
    return(CaGenomeAnnotation)
  } else {
    # subset the data
    gene_annotation <- CaGenomeAnnotation[CaGenomeAnnotation$geneid %in% geneids, ]
    return(gene_annotation)
  }
}

