#' @title Compute the Determinant of a Matrix
#' @description \code{Det} computes the determinant of a square matrix.
#' If the matrix is complex, the determinant is computed as the product of the eigenvalues; if the matrix
#' is real, \code{Det} calls the base function \code{det} for maximum efficiency.
#' @param M a square matrix, real or complex.
#' @return The determinant of M.
#' @author Albert Dorador
#' @export
#' @examples
#' A <- matrix(c(1, 2, 2+3i, 5), ncol = 2) #complex matrix
#' B <- matrix(1:4, ncol = 2) #real matrix
#'
#' Det(A)
#' Det(B)
#'
Det <- function(M){
  if (nrow(M) != ncol(M))
    stop("Supplied matrix is not square")
  if(is.complex(M) == TRUE){
    prod(eigen(M, only.values = TRUE)$values)
  }else{
    det(M)
  }
}
