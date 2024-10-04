rm(list=ls())
library(ChickpeaOmicsR)
## # Load data
#data("CaExpressionCounts")
#data("CaExpressionMetadata")
#data("CaGenomeAnnotation")
#data("CaProteinEnrichment")









salinity<-get_study_data("Salinity")
counts<-salinity$counts
metadata<-salinity$metadata
# get only genes in "F-box" family
genes<-get_gene_express_gokeywords(c("F-box"),counts)
# remove genes with low expression
genes<-genes[rowSums(genes)>100,]

# reaplace NA with 0
genes[is.na(genes)]<-0
# remove genes with low expression
genes<-genes[rowSums(genes)>100,]
genes
library(pheatmap)
library(limma)
# color blue yellow red


heatmap(genes, scale="row", col=heat.colors(256))
# color blue yellow red
pheatmap(genes, scale="row", col=colorRampPalette(c("blue", "yellow", "red"))(256))


get_gene_annotation(c("GRX3", "HIS2A","HSFA1B","HSFA2","HSFA3"))
data<-get_gene_express_bygoid(c("GO:0030246"))
nrow(data)
data
get_goids_bykewords(c("targeting"))
get_gene_express_gotermkeywords(c("targeting"))


get_geneids_bygoids(c("CL:30078"))
