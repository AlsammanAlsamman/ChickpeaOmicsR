#' Gene Expression Counts
#' Gene expression counts for 28891 gene and 197 samples
#' @format matrix with 28891 rows and 197 columns
#' @source https://www.ncbi.nlm.nih.gov
#' @examples
#' data(CaExpressionCounts)
"CaExpressionCounts"

#' Gene Expression Metadata
#' Metadata for 197 samples across 14 different studies and 6 different stress conditions
#' @format data.frame with 197 rows and 6 columns
#' @source https://www.ncbi.nlm.nih.gov
#' @examples
#' data(CaExpressionMetadata)
"CaExpressionMetadata"

#' Gene annotation data
#' Annotation data for 24681 genes
#' @format data.frame with 24681 rows and 6 columns
#' @examples
#' data(CaGenomeAnnotation)
"CaGenomeAnnotation"

#' Protein enrichment data
#' Protein enrichment data for 24254 genes
#' @format data.frame with 24254 rows and 5 columns
#' @examples
#' data(CaProteinEnrichment)
"CaProteinEnrichment"

#' Protein sequence information
#' Protein sequence information for 24681 genes
#' @format data.frame with 24681 rows and 5 columns
#' @examples
#' data(CaProteinSequenceInfo)
#' head(CaProteinSequenceInfo)
"CaProteinSequenceInfo"

#' Protein interaction network
#' Protein interaction network for 24681 genes
#' @format data.frame 9 columns
#' @examples
#' data(CaProteinInteractionNetwork)
"CaProteinInteractionNetwork"

#' Significance of gene expression
#' Significance of gene expression for 28891 genes
#' @format data.frame "geneid", "pvalue", "FDR", "logfc", "exp.direction", "study", "experiment"
#' @examples
#' data(CaExpressionSignificance)
"CaExpressionSignificance"

#' GWAS data
#' GWAS data for 1052 genes
"CaGWAS"