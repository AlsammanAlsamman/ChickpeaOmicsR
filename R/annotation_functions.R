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

#' Get gene enrichment data for a specific geneid
#' @param geneids The gene or genes to get the enrichment for
#' @return Gene enrichment data
#' @export
get_gene_enrichment_bygeneid<-function(geneids){
  loadCaProteinEnrichment()
  gene_enrichment<-CaProteinEnrichment[CaProteinEnrichment$geneid %in% geneids,]
  return(gene_enrichment)
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


#' Get gene ids by GO term keywords
#' @param keywords The keywords to search for
#' @return Gene ids
#' @export
#' @examples
#' get_geneids_bygotermkeywords(c("targeting"))
get_geneids_bygotermkeywords<-function(keywords){
  goids<-get_goids_bykewords(keywords)
  geneids<-get_geneids_bygoids(goids)
  return(geneids)
}

#' Get gene ids by list of Go ids
#' @param goids The GO ids to search for
#' @return Gene ids
#' @export
get_geneids_bygoids<- function(goids){
  loadCaProteinEnrichment()
  #create a list to store the gene ids
  geneids<-list()
  for (goid in goids){
    geneids[[goid]]<-CaProteinEnrichment[CaProteinEnrichment$termid==goid,]$geneid
  }
  return(geneids)
}
