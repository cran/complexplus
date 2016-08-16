#' @title Matrix Logarithm
#' @description \code{matlog} computes the (principal) matrix logarithm of a square matrix.
#' @details This function adapts function \code{logm} from package \pkg{expm}
#' to be able to handle complex matrices, by decomposing the original matrix
#' into real and purely imaginary matrices and creating a real
#' block matrix that function \code{logm} can successfully process.
#' If the original matrix is real, \code{matlog} calls \code{logm} directly for maximum efficiency.
#' Hence, for real matrices, \code{matlog} can compute the matrix logarithm in the same instances as \code{logm};
#' for complex matrices, \code{matlog} can compute the matrix logarithm as long as all real
#' eigenvalues are positive: zero eigenvalues imply singularity (and therefore the log
#' does not exist) and negative eigenvalues can be problematic as it may be hard and
#' numerically unstable to calculate Jordan blocks. See references below.
#' @references For more on the matrix logarithm, visit
#' \url{https://en.wikipedia.org/wiki/Logarithm_of_a_matrix}
#' @param A a square matrix, real or complex.
#' @param ... arguments passed to or from other methods.
#' @return The matrix logarithm of A. Method used may be chosen from the options available in \code{logm}.
#' @author Albert Dorador
#' @export
#' @import expm
#' @seealso \code{\link[expm]{logm}}
#' @examples
#' A <- matrix(c(1, 2, 2+3i, 5), ncol = 2)  # complex matrix
#' B <- matrix(c(2, 0, 3, 4), ncol = 2)  # real matrix with positive eigenvalues
#'
#' matlog(A)
#' matlog(A, "Eigen")  # uses Eigen method in function logm
#' matlog(B)
#' matlog(B, "Eigen")  # uses Eigen method in function logm
#'
matlog <- function(A, ...)
{
  d <- dim(A)
  n <- d[1]
  p <- d[2]
  if (n != p)
    stop("Supplied matrix is not square")
  if (is.complex(A) == TRUE){
    counter <- integer(1) #integer takes half the memory space of numeric (double type)
    e.val <- eigen(A, only.values = TRUE)$values
    for (i in 1:n){
      if(class(Imzap(e.val[i])) != "complex"){ #i.e. if eigenvalue is real
        if(Imzap(e.val[i]) > 0){ #0 would make the matrix singular, and negative is likely to cause trouble
          counter <- counter + 1
        }
      }else{
        counter <- counter + 1 #if complex eigenvalue, always OK individually
      }
    }
    if (counter == n){ #i.e. only if taking the log is safe
      # Extract real and imaginary parts
      Ar <- Re(A)
      Ai <- Im(A)

      # Construct real but extended matrix
      L <- rbind( cbind(Ar,-Ai), cbind(Ai,Ar))

      # Compute exponential of that
      lL <- logm(L, ...)

      # Construct logm(A) by extracting real and imaginary parts from blocks
      lA <- lL[1:n,1:n] + 1i*lL[(1:n)+n,1:n]

      lA <- matrix(Imzap(lA), ncol = n) # to guarantee best-looking solutions

      return(lA)
    }else{
      stop("The matrix logarithm may not exist for this matrix.")
    }
  }else{
    logm(A, ...)
  }
}
