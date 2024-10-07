# Function to map expression values to 5 negative and 10 positive scales
map_expression_to_color <- function(expr, min_expr, max_expr) {
  # Define 5 color scale for negative values (blue to white)
  negative_colors <- colorRampPalette(c('blue', 'white'))(10)
  # Define 10 color scale for positive values (white to red)
  positive_colors <- colorRampPalette(c('yellow', 'red'))(10)
  
  # Map the expression values to discrete scales
  if (expr < 0) {
    # Scale negative values between min_expr and 0 (5 bins)
    bin <- cut(expr, breaks = seq(min_expr, 0, length.out = 6), include.lowest = TRUE, labels = FALSE)
    return(negative_colors[bin])
  } else {
    # Scale positive values between 0 and max_expr (10 bins)
    bin <- cut(expr, breaks = seq(0, max_expr, length.out = 11), include.lowest = TRUE, labels = FALSE)
    return(positive_colors[bin])
  }
}

#' Plot a protein interaction network with gene expression data
#' @param network_table The protein interaction network table
#' @param geneExp The gene expression data
#' @param targetCol The column in geneExp to use for coloring nodes
#' @param title The title of the plot
#' @param file The file path to save the plot (optional)
#' @param plot_width The width of the plot (default: 10)
#' @param plot_height The height of the plot (default: 10)
#' @param minnodes The minimum number of connections a node must have to be included in the plot (default: 3)
#' @return A plot of the protein interaction network with gene expression data
#' @export

plot_protein_interaction_network_with_gene_expression <- function(network_table, geneExp, targetCol="score", title="Protein Interaction Network", file=NULL, plot_width = 10, plot_height=10, minnodes=3) {
  # Select only the required columns for the edges
  edges <- network_table[, c("p1_geneid", "p2_geneid", targetCol)]
  #just to be sure
  colnames(geneExp)<-c("geneid","expression")
  # Create a graph object from the edges
  g <- graph_from_data_frame(d = edges, directed = FALSE)
  # print a note that some nodes were removed if minnodes > 0
  if (minnodes > 0) {
    #remove nodes with less than 3 connections
    g<-delete_vertices(g, V(g)[degree(g) < minnodes])
    print(paste("Note: Removed nodes with less than", minnodes, "connections"))
  }
  # Find the min and max of the expression values for scaling
  min_expr <- min(geneExp$expression, na.rm = TRUE)
  max_expr <- max(geneExp$expression, na.rm = TRUE)
  
  # Check if all gene names are present in the geneExp table
  unmatched_genes <- setdiff(V(g)$name, geneExp$geneid)
  if (length(unmatched_genes) > 0) {
    print(paste("Warning: The following genes are missing expression data:", paste(unmatched_genes, collapse = ", ")))
  }
  
  # Map gene expression to node colors using the new mapping function
  V(g)$color <- sapply(V(g)$name, function(gene) {
    expr <- geneExp$expression[geneExp$geneid == gene]
    if (length(expr) > 0 && !is.na(expr)) {
      # Call the custom function to map expression to color
      map_expression_to_color(expr, min_expr, max_expr)
    } else {
      "grey"  # If no expression data is found, color the node grey
    }
  })
  
  # Normalize the 'score' column to set as edge lengths (inverse, as higher scores = longer distances)
  edge_lengths <- 1 / (E(g)$score + 0.01)  # Avoid division by zero by adding a small constant
  
  # Assign edge colors based on the 'score' or 'coexpression' values (scaling them similarly)
  E(g)$color <- apply(colorRamp(c('blue', 'yellow', 'red'))((E(g)$score - min(E(g)$score)) / (max(E(g)$score) - min(E(g)$score))), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
  
  # Use a layout that takes edge weights (lengths) into account
  layout <- layout_with_fr(g, weights = edge_lengths)
  
  # Plot the network with node and edge colors
  plot(g, layout = layout, vertex.color = V(g)$color, edge.color = E(g)$color, 
       vertex.size = 10, edge.width = 2)
  
  # Add a title to the plot
  title(main = title)
  
  # Add note that edge lengths represent scores
  mtext(paste("Edge lengths represent", targetCol), side = 1, line = 3, cex = 0.8, col="red")
  
  if (!is.null(file)) {
    # Save the plot as a file
    dev.copy(pdf, file, width = plot_width, height = plot_height)
    dev.off()
  }
}

#' Generate a heatmap for a given GO term keyword
#' @param keyword The keyword to search for
#' @param study The study to filter the gene expression data
#' @param experiment The experiment to filter the gene expression data
#' @param outfolder The folder to save the heatmap and data files
#' @param plotwidth The width of the plot (default: 10)
#' @param plotheight The height of the plot (default: 10)
#' @import pheatmap
#' @export
#' @examples
#' generate_heatmap_by_go_term_keyword(c("Flowering"),study="Salinity", experiment=c("PRJNA232700","PRJNA288473"), outfolder="heatmaps", plotwidth=10, plotheight=20)
generate_heatmap_by_go_term_keyword<-function(keyword, study=NULL, experiment=NULL, outfolder="", plotwidth=10, plotheight=10){
  # get the gene ids related to the keywords
  genes<-get_geneids_bygotermkeywords(keyword)
  # get the gene annotation data
  gene_annotation<-get_gene_annotation_bygeneid(unlist(genes))
  # get the gene enrichment data
  gene_enrichment<-get_gene_enrichment_bygeneid(unlist(genes))
  # get the gene expression data
  genesExpress<-get_gene_express_bygeneid(unlist(genes), study=study, experiment=experiment)
  # filter the gene expression data
  genesExpress<-filter_gene_expression(genesExpress, min_expression=10)
  # get the metadata for the gene expression data
  metadata<-get_metadata(colnames(genesExpress))
  # select only study and experiment columns
  metadata<-metadata[,c("group","study","experiment")]
  # plot the heatmap
  pheatmap(as.matrix(genesExpress), scale="row",annotation_col=metadata,
           #color scale blue yellow red
           col=colorRampPalette(c("blue", "yellow", "red"))(256))
  if (outfolder!="") {
    # if folder does not exist create it
    if (!file.exists(outfolder)) {
      dir.create(outfolder)
    }
    # get the last on device
    dev.copy(pdf, paste(outfolder, paste0(keyword, "_heatmap.pdf"), sep="/"), width=plotwidth, height=plotheight)
    dev.off()
    # save annotation data
    write.csv(gene_annotation, paste(outfolder, paste0(keyword, "_annotation",".csv"), sep="/"))
    # save enrichment data
    write.csv(gene_enrichment, paste(outfolder, paste0(keyword, "_enrichment",".csv"), sep="/"))
  }
}

#' Get gene expression data using GO term keywords
#' @param keywords The keywords to search for
#' @param pvlaue The pvalue threshold
#' @param phenotype The phenotype to filter the data
#' @return gene Pheno matrix from GWAS
#' @import reshape2
#' @export
gene_pheno_GWASmatrix_byGOKeyword<-function(keyword, pvalue=0.05,phenotype=NULL){
  genes<-get_geneids_bygotermkeywords(keyword)
  genes<-unlist(genes)
  GWASdata<-get_gwas_data(genes,pvalue = pvalue, phenotype = phenotype)
  gene_pheno_matrix<-reshape2::dcast(GWASdata, geneid~phenotype, value.var="pvalue")
  return(gene_pheno_matrix)
}


#' Plot gene-phenotype heatmap from GWAS data
#' @param gene_pheno_matrix The gene-phenotype matrix
#' @param outfolder The folder to save the plot
#' @param plotwidth The width of the plot (default: 10)
#' @param plotheight The height of the plot (default: 10)
#' @param keyword The keyword to search for
#' @param pvalue The pvalue threshold
#' @param phenotype The phenotype to filter the data
#' @import ggplot2
#' @import reshape2
#' @export
plot_gene_pheno_GWASmatrix_byGOKeyword<-function(gene_pheno_matrix=NULL, outfolder=NULL, plotwidth=10, plotheight=10, keyword=NULL, pvalue=0.05,phenotype=NULL)
{
  if (is.null(gene_pheno_matrix)) {
    gene_pheno_matrix<-gene_pheno_GWASmatrix_byGOKeyword(keyword, pvalue=pvalue,phenotype=phenotype)
  }
  gene_pheno_matrix <- gene_pheno_GWASmatrix_byGOKeyword("Flowering", pvalue = 0.05)
  gene.phenoCount<-reshape2::melt(gene_pheno_matrix)
  colnames(gene.phenoCount)<-c("gene","phenotype","value")
  p<-ggplot(gene.phenoCount, aes(x=gene, y=phenotype, size=value, color=as.factor(phenotype))) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Gene-Phenotype Frequency Heatmap", x="Gene", y="Phenotype")
  
  if (!is.null(outfolder)) {
    if (!file.exists(outfolder)) {
      dir.create(outfolder)
    }
    ggplot2::ggsave(p, file=paste(outfolder, "gene_pheno_heatmap.pdf", sep="/"), width=plotwidth, height=plotheight)
  }
}


