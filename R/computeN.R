#' @title Compute N
#' 
#' @description Calculates N given gamma.
#' 
#' @details Computation is performed following Theorem 15 of Rubinstein, Benjamin IP, and Francesco Aldà. "Pain-free random differential privacy with sensitivity sampling." International Conference on Machine Learning. PMLR, 2017.
#'
#' @param gamma \code{numeric} Value of gamma to compute N
#'
#' @return A number to be used as N on re-sampling methods.
#' @export
computeN <- function(gamma) {
  rho <- min(gamma, 0.5)
  m <- (1)/(2 * (gamma - rho)^2) * log10(1/rho)
  return(m)
}