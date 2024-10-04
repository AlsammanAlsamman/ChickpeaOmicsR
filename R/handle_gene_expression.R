#' Filter Gene Expression Data
#' Filter gene expression data based minimum expression level
#' @param data gene expression matrix
#' @param minexp minimum expression level
#' @return filtered gene expression matrix
#' @examples
#' data("CaExpressionCounts")
filter_gene_expression<-function(data){
  data<-genes
  # convert NA to 0
  data[is.na(data)]<-0
  # convert list to linear
  datalinear<-as.vector(data)
  datalinear<-unique(datalinear)
  #get the mean expression using boxplot
  meanexp<-boxplot(datalinear)$stats[3]
  maxexp<-boxplot(datalinear)$stats[5]
  minexp<-boxplot(datalinear)$stats[1]
  # values less the minimum expression level convert to 0
  data[data<minexp]<-0
  # values greater than the maximum expression level convert to 0
  data[data>maxexp]<-0
  # plot the boxplot
  boxplot(data)
  
  return(data)
}
