% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/laplaceNoise.R
\name{Laplace_noise_generator}
\alias{Laplace_noise_generator}
\title{Random laplace generator}
\usage{
Laplace_noise_generator(mu, b, n.noise)
}
\arguments{
\item{mu}{\code{integer} Mean of Laplace distribution}

\item{b}{\code{integer} Standard deviation of Laplace distribution}

\item{n.noise}{\code{integer} Number of samples to generate}
}
\value{
Vector of numbers
}
\description{
Extracted from https://rpubs.com/mengxu/gaussian-and-laplace-noise
}
