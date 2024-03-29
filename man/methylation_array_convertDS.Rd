% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylation_array_convertDS.R
\name{methylation_array_convertDS}
\alias{methylation_array_convertDS}
\title{Convert a Methylation ExpressionSet from EPIC to 450k Array}
\usage{
methylation_array_convertDS(eSet, objective_array)
}
\arguments{
\item{eSet}{A methylation ExpressionSet object containing either EPIC or 450k array data.}

\item{objective_array}{A character string specifying the target array type. Can be either "450k" or "epic".}
}
\value{
An ExpressionSet object containing the subset of probes based on the target array type.
}
\description{
This function subsets a methylation ExpressionSet from an EPIC array (850k probes) to a 450k array or vice versa.
}
\examples{
# Load example data as 'eSet'
# Convert the eSet from EPIC array to 450k array
eSet_450k <- methylation_array_convertDS(eSet, objective_array = "450k")
# Convert the eSet from 450k array to EPIC array
eSet_epic <- methylation_array_convertDS(eSet, objective_array = "epic")
}
