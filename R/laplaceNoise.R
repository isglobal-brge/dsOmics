#' @title Random laplace generator
#' 
#' @description Extracted from https://rpubs.com/mengxu/gaussian-and-laplace-noise
#'
#' @param n \code{integer} (default \code{1}) Number of samples to generate
#' @param m \code{integer} (default \code{0}) Mean of Laplace distribution
#' @param s \code{integer} (default \code{1}) Standard deviation of Laplace distribution
#'
#' @return Vector of numbers
#' @export

rlaplace <- function (n = 1, m = 0, s = 1){
  if (any(s <= 0)) 
    stop("s must be positive")
  q <- runif(n)
  ifelse(q < 0.5, s * log(2 * q) + m, -s * log(2 * (1 - q)) + m)
}

#' @title Random laplace distribution generator
#' 
#' @description Extracted from https://rpubs.com/mengxu/gaussian-and-laplace-noise
#' 
#' @param y \code{numeric} Samples
#' @param m \code{integer} (default \code{0}) Mean of Laplace distribution
#' @param s \code{integer} (default \code{1}) Standard deviation of Laplace distribution
#' @param log \code{bool} (default \code{FALSE}) Apply \code{exp} transformation to the distribution
#'
#' @return Vector of numbers
#' @export

dlaplace <- function (y, m = 0, s = 1, log = FALSE) {
  if (any(s <= 0)) 
    stop("s must be positive")
  density <- -abs(y - m)/s - log(2 * s)
  if (!log) 
    density <- exp(density)
  density
}

#' @title Random laplace generator
#' 
#' @description Extracted from https://rpubs.com/mengxu/gaussian-and-laplace-noise
#' 
#' @param mu \code{integer} Mean of Laplace distribution
#' @param b \code{integer} Standard deviation of Laplace distribution
#' @param n.noise \code{integer} Number of samples to generate
#'
#' @return Vector of numbers
#' @export

Laplace_noise_generator <- function(mu, b, n.noise){
  # generate n noise points
  noise.points <- rlaplace(n = n.noise, m = mu, s = b)
  noise.density <- dlaplace(noise.points, m = mu, s = b)
  
  return(noise = noise.points)
}