# Laplace noise generator from https://rpubs.com/mengxu/gaussian-and-laplace-noise

#' Title
#'
#' @param n 
#' @param m 
#' @param s 
#'
#' @return
#' @export
#'
#' @examples
rlaplace <- function (n = 1, m = 0, s = 1){
  if (any(s <= 0)) 
    stop("s must be positive")
  q <- runif(n)
  ifelse(q < 0.5, s * log(2 * q) + m, -s * log(2 * (1 - q)) + m)
}

#' Title
#'
#' @param y 
#' @param m 
#' @param s 
#' @param log 
#'
#' @return
#' @export
#'
#' @examples
dlaplace <- function (y, m = 0, s = 1, log = FALSE) {
  if (any(s <= 0)) 
    stop("s must be positive")
  density <- -abs(y - m)/s - log(2 * s)
  if (!log) 
    density <- exp(density)
  density
}

#' Title
#'
#' @param mu 
#' @param b 
#' @param n.noise 
#'
#' @return
#' @export
#'
#' @examples
Laplace_noise_generator <- function(mu, b, n.noise){
  # generate n noise points
  noise.points <- rlaplace(n = n.noise, m = mu, s = b)
  noise.density <- dlaplace(noise.points, m = mu, s = b)
  
  return(noise = noise.points)
}