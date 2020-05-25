# cobiclust R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRAE, Universite Paris-Saclay, 75005, Paris, France
####################################################################################
#' Creation of the cobiclustering class.
#' @export cobiclustering
#' @keywords internal


cobiclustering <- function(data = matrix(nrow = 3, ncol = 3, NA), K = 2, G = 2, classification = list(length=2), strategy = list(), parameters = list(), info = list())
{

  object <- list(
    data = "matrix",
    K = "numeric",
    G = "numeric",
    classification = "list",
    strategy = "list",
    parameters = "list",
    info = "list"
  )
  #UseMethod("new",object)
  ## Set the name for the class
  class(object) <- "cobiclustering"
  return(object)
}

#' Is an object of class cobiclustering ?
#'
#' @param object an object of class cobiclustering.
#' @export is.cobiclustering
#' @keywords internal
is.cobiclustering <- function(object) {
  return(class(object) == "cobiclustering")
}

#' Summary of an object of class Cobiclust
#'
#' @param object an object of class cobiclustering.
#' @param ... ignored
#' @keywords internal
#' @method summary cobiclustering
#' @export
summary.cobiclustering <- function (object, ...){
  cat("-----------------------------------------------------------\n")
  cat("Number of biclusters: K = ", object$K, "* G = ", object$G,"\n")
  cat("Row proportions:", round(object$parameters$pi, 3) ,"\n")
  cat("Column proportions:", round(object$parameters$rho, 3) ,"\n")
  cat("-----------------------------------------------------------\n")
  cat("Matrix of alpha_kg interactions terms: \n")
  print(matrix(nrow = object$K, ncol = object$G, round(object$parameters$alpha, 3)))
  cat("-----------------------------------------------------------\n")
  cat("Frequencies of MAP classification \n")
  cat("Rows")
  print(table(object$classification$rowclass))
  cat("Columns")
  print(table(object$classification$colclass))
  cat(" \n")
  NextMethod("summary",object)
  cat("-----------------------------------------------------------\n")
}
