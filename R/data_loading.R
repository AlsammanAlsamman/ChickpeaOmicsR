#' @title Loading required data
#' @description ChickpeaOmicsR is a package that provides gene expression data for chickpea under different stress conditions.
#' @export
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ChickpeaOmicsR!")
  print("This package provides gene expression data for chickpea under different stress conditions.")
  print("Several Data MayBe Loaded to the Environment.")
  print("Gene Expression counts are normalized")
  print("Alsamman M. Alsamman, aalsamman100[at and not .]gmail[now .]com")
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
loadCaGenomeAnnotation <- function() {
  # check if the data is already loaded
  if (!exists("CaGenomeAnnotation")) {
    # create global variable of the data
    data(CaGenomeAnnotation)
    CaGenomeAnnotation<<-CaGenomeAnnotation
  }
}

#' Load CaProteinEnrichment
#' @export
loadCaProteinEnrichment <- function() {
  # check if the data is already loaded
  if (!exists("CaProteinEnrichment")) {
    # create global variable of the data
    data(CaProteinEnrichment)
    CaProteinEnrichment<<-CaProteinEnrichment
  }
}

#' Load CaExpressionSignificance
#' @export
loadCaExpressionSignificance <- function() {
  # check if the data is already loaded
  if (!exists("CaExpressionSignificance")) {
    # create global variable of the data
    data(CaExpressionSignificance)
    CaExpressionSignificance<<-CaExpressionSignificance
  }
}

#' Load CaGWAS
#' @export
loadCaGWAS <- function() {
  # check if the data is already loaded
  if (!exists("CaGWAS")) {
    # create global variable of the data
    data(CaGWAS)
    CaGWAS<<-CaGWAS
  }
}