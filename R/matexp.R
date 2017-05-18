#' @title Matrix Exponential
#' @description \code{matexp} computes the exponential of a square matrix A i.e. \eqn{exp(A)}.
#' @details This function adapts function \code{expm} from package \pkg{expm}
#' to be able to handle complex matrices, by decomposing the original
#' matrix into real and purely imaginary matrices and creating a real
#' block matrix that function \code{expm} can successfully process.
#' If the original matrix is real, \code{matexp} calls \code{expm} directly for maximum efficiency.
#' @param A a square matrix, real or complex.
#' @param ... arguments passed to or from other methods.
#' @return The matrix exponential of A. Method used may be chosen from the options available in \code{expm}.
#' @author Uffe Høgsbro Thygesen
#' @export
#' @import expm
#' @seealso \code{\link[expm]{expm}}
#' @examples
#' A <- matrix(c(1, 2, 2+3i, 5), ncol = 2)  # complex matrix
#' B <- matrix(1:4, ncol = 2)  # real matrix
#'
#' matexp(A)
#' matexp(A, "Ward77")  # uses Ward77's method in function expm
#' matexp(B)
#' matexp(B, "Taylor")  # uses Taylor's method in function expm
#'
matexp <- function(A, ...)
{
  d <- dim(A)
  n <- d[1]
  p <- d[2]
  if (n != p)
    stop("Supplied matrix is not square")
  if (is.complex(A)){
    # Extract real and imaginary parts
    Ar <- Re(A)
    Ai <- Im(A)

    # Construct real but extended matrix
    E <- rbind( cbind(Ar,-Ai), cbind(Ai,Ar))

    # Compute exponential of that
    eE <- expm(E, ...)

    # Construct expm(A) by extracting real and imaginary parts from blocks
    eA <- eE[1:n,1:n] + 1i*eE[(1:n)+n,1:n]

    eA <- matrix(Imzap(eA), ncol = n) #to guarantee best-looking solutions

    return(eA)
  }else{
    expm(A, ...)
  }
}
