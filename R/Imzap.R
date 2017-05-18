#' @title Rounding of Null Imaginary Part of a Complex Number
#' @description imaginary parts with values very close to 0 are 'zapped', i.e. treated as 0.
#' Therefore, the number becomes real and changes its class from complex to numeric.
#' @param x a scalar or vector, real or complex.
#' @param tol a tolerance, \eqn{10^{-6}}{10^-6} by default. Prevents possible numerical problems.
#' Can be set to 0 if desired.
#' @author Albert Dorador
#' @export
#' @examples
#' x1 <- 1:100
#' x2 <- c(1:98,2+3i,0-5i)
#' x3 <- c(1:98,2+0i,7-0i)
#' x4 <- complex(real = rnorm(100), imaginary = rnorm(100))
#'
#' Imzap(x1)  # inocuous with real vectors
#' Imzap(x2)  # 1 single element is enough to turn the vector into complex
#' Imzap(x3)  # removes extra 0i and changes class from from complex to numeric
#' Imzap(x4)  # inocuous with complex vectors with non-null complex part
#'
Imzap <- function(x, tol = 1e-6) {
  if (all(abs(Im(z <- zapsmall(x))) <= tol)) as.double(x) else x
}
