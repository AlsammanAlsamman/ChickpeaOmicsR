#' Load the gene expression count data
#' @export
loadCaProteinInteractionNetwork <- function() {
  # check if the data is already loaded
  if (!exists("CaProteinInteractionNetwork")) {
    # load the data
    data(CaProteinInteractionNetwork)
    CaProteinInteractionNetwork<<-CaProteinInteractionNetwork
  }
}

#' Get the protein interaction network table for specific genes
#' @param genes a character vector of gene names
#' @return a data frame with protein interaction network information
#' @examples
#' get_protein_interaction_table(c("LOC101501102","LOC101507179"))
#' @export
get_protein_interaction_table <- function(genes, only_interacting=TRUE) {
  # load the data
  loadCaProteinInteractionNetwork()
  if (only_interacting) {
    # get the protein interaction network for the genes where p1, and p2 are in the list of genes
    network<-CaProteinInteractionNetwork[CaProteinInteractionNetwork$p1_geneid%in%genes & CaProteinInteractionNetwork$p2_geneid%in%genes,]
  } else {
    # get the protein interaction network for the genes where p1, or p2 is in the list of genes
    network<-CaProteinInteractionNetwork[CaProteinInteractionNetwork$p1_geneid%in%genes | CaProteinInteractionNetwork$p2_geneid%in%genes,]
  }
  return(network)
}

#' Plot a protein interaction network without gene expression colors
#' @param network_table a data frame with protein interaction network information
#' @param targetCol the column name for the target score
#' @param title the title for the plot
#' @import igraph
#' @export
plot_protein_interaction_network <- function(network_table, targetCol="score", title="Protein Interaction Network", file=NULL,plot_width = 10, plot_height=10) {
  # select only the required columns for the edges
  edges <- network_table[, c("p1_geneid", "p2_geneid", targetCol)]
  # Create a graph object from the edges
  g <- graph_from_data_frame(d = edges, directed = FALSE)
  # Normalize the 'score' column to set as edge lengths (inverse, as higher scores = longer distances)
  edge_lengths <- 1 / (E(g)$score + 0.01)  # Avoid division by zero by adding a small constant
  # use a layout that takes edge weights (lengths) into account
  layout <- layout_with_fr(g, weights = edge_lengths)
  # plot the network without node colors
  plot(g, layout = layout, edge.width = 2)
  # add a title to the plot
  title(main = title)
  # add note that edge lengths represent scores
  mtext(paste("Edge lengths represent", targetCol), side = 1, line = 3, cex = 0.8,col="red")
  if (!is.null(file)) {
    # save the plot as a file
    dev.copy(pdf, file, width = plot_width, height = plot_height)
    dev.off()
  }
}