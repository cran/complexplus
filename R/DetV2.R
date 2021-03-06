#' @title Compute the Determinant of a Matrix
#' @description \code{Det} computes the determinant of a square matrix.
#' This function first checks whether the matrix is full rank or not; if not,
#' the value 0 is returned. This avoids relatively frequent numerical errors
#' that produce a non-zero determinant when in fact it is zero.
#' Only if the matrix is full rank does the algorithm proceed to compute the determinant.
#' If the matrix is complex, the determinant is computed as the product of the eigenvalues; if the matrix
#' is real, \code{Det} calls the base function \code{det} for maximum efficiency.
#' @param M a square matrix, real or complex.
#' @return The determinant of M.
#' @author Albert Dorador
#' @export
#' @importFrom Matrix rankMatrix
#' @examples
#' A <- matrix(c(1, 2, 2+3i, 5), ncol = 2) #complex matrix
#' B <- matrix(1:4, ncol = 2) #real matrix
#' S <- matrix(c(3, 4+3i, 0, 0), ncol = 2) #Singular matrix
#'
#' Det(A)
#' Det(B)
#' Det(S)
#'
Det <- function(M){
  nrw = nrow(M)
  if (nrw != ncol(M))
    stop("Supplied matrix is not square")
  if (rankMatrix(M) != nrw){
    return(0)
  }
  if(is.complex(M)){
    prod(eigen(M, only.values = TRUE)$values)
  }else{
    det(M)
  }
}
