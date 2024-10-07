
#' Get gene expression counts for a specific geneid
#' @param genenames The gene or genes to get the counts for
#' @return Expression Count matrix for the gene or genes
#' @export
#' @examples
#' get_gene_express_bygeneid(c("GRX3", "HIS2A","HSFA1B","HSFA2","HSFA3"))
get_gene_express_bygeneid <- function(geneids, gene_express_count=NULL, study=NULL, experiment=NULL){
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
  gene_express <- geneExpressionCounts[rownames(geneExpressionCounts) %in% geneids, ]
  #if study or experiment is provided
  if (!is.null(study)) {
   # load metadata
    loadCaExpressionMetadata()
    # get samples of this study
    studysamples<-rownames(CaExpressionMetadata[CaExpressionMetadata$study == study, ])
    gene_express <- gene_express[, colnames(gene_express) %in% studysamples]
  }
  if (!is.null(experiment)) {
    # load metadata
    loadCaExpressionMetadata()
    # get samples of this experiment
    experimentsamples<-rownames(CaExpressionMetadata[CaExpressionMetadata$experiment == experiment, ])
    gene_express <- gene_express[, colnames(gene_express) %in% experimentsamples]
  }
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
  gene_express <- CaExpressionCounts[, colnames(CaExpressionCounts) %in% samples]
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
  gene_express <- geneExpressionCounts[rownames(geneExpressionCounts) %in% geneids, ]
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
  gene_express<-GeneExpressionCounts[rownames(GeneExpressionCounts) %in% geneids,]
  return(gene_express)
}

#' Filter gene expression data based minimum expression level
#' @param data gene expression matrix
#' @return filtered gene expression matrix
#' @export
#' @examples
#' data("CaExpressionCounts")
filter_gene_expression<-function(data, min_expression=10){
  # replace NA with 0
  data[is.na(data)]<-0
  #remove genes with 0 expression
  data<-data[rowSums(data)>min_expression,]
  return(data)
}

#' Plot gene expression heatmap
#' @param geneExp Gene expression data
#' @param plot_height The height of the plot
#' @param plot_width The width of the plot
#' @param file The file to save the plot
#' @import ggplot2
#' @export
#' @return The plot
plot_gene_expression_heatmap<-function(geneExp,plot_height=10,plot_width=10,file="gene_expression_heatmap.pdf"){
  p<-ggplot(geneExp, aes(x=geneid, y=study, fill=logfc)) + 
    geom_tile() +
    # scale blue yellow red
    scale_fill_gradient(low = "blue", high = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Gene Expression Heatmap", x="Gene", y="Study")
  ggsave(file, plot=p, height=plot_height, width=plot_width)
  print(paste("Heatmap saved to ", file))
  p
}
############## Dealing with Gene expression analysis

#' Get gene expression significance data
#' @return Gene logfc data
#' @param genesSample The gene or genes to get the logfc for
#' @param pvalue The pvalue threshold
#' @param fdr The fdr threshold
#' @param study The study to get the logfc for
#' @param experiment The experiment to get the logfc for
#' @export
get_gene_logfc<-function(genesSample, pvalue=1, fdr=1, study="all", experiment="all"){
  # load gene expression data if not loaded
  loadCaExpressionSignificance()
  # select by geneid
  geneExp<-CaExpressionSignificance[CaExpressionSignificance$geneid %in% genesSample,]
  # if pvalue is not all
  if (pvalue != 1) {
    geneExp<-geneExp[geneExp$pvalue<=pvalue,]
    print("some genes may have been removed due to pvalue")
  }
  # if fdr is not 1
  if (fdr != 1) {
    geneExp<-geneExp[geneExp$fdr<=fdr,]
    print("some genes may have been removed due to fdr")
  }
  # if study is not all
  if (study != "all") {
    geneExp<-geneExp[geneExp$study==study,]
  }
  # if experiment is not all
  if (experiment != "all") {
    geneExp<-geneExp[geneExp$experiment==experiment,]
  }
  # if the gene have several logfc values, get the mean
  print("some genes may have several logfc values, the mean is taken")
  geneExp<-aggregate(logfc~geneid, geneExp, mean)
  return(geneExp)
}

#' Get the mean logfc for gene list across categories
#' @param geneids The gene or genes to get the logfc for
#' @param pvalue The pvalue threshold
#' @param fdr The fdr threshold
#' @param category The category to get the logfc for
#' @export
get_gene_logfc_bycategory<-function(geneids, pvalue=1, fdr=1, category="study"){
  # load gene expression data if not loaded
  loadCaExpressionSignificance()
  # select by geneid
  geneExp<-CaExpressionSignificance[CaExpressionSignificance$geneid %in% geneids,]
  # if pvalue is not all
  if (pvalue != 1) {
    geneExp<-geneExp[geneExp$pvalue<=pvalue,]
    print("some genes may have been removed due to pvalue")
  }
  # if fdr is not 1
  if (fdr != 1) {
    geneExp<-geneExp[geneExp$fdr<=fdr,]
    print("some genes may have been removed due to fdr")
  }
  # if category is study
  if (category == "study") {
    geneExp<-aggregate(logfc~geneid+study, geneExp, mean)
  }
  # if category is experiment
  if (category == "experiment") {
    geneExp<-aggregate(logfc~geneid+experiment, geneExp, mean)
  }
  return(geneExp)
}
