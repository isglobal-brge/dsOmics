% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limmaDS.R
\name{limmaDS}
\alias{limmaDS}
\title{Differential expression analysis using limma on the server-side}
\usage{
limmaDS(
  Set,
  variable_names,
  covariable_names,
  type,
  contrasts,
  levels,
  coef,
  sva,
  annotCols,
  method,
  robust,
  normalization,
  voomQualityWeights,
  big,
  sort.by
)
}
\arguments{
\item{Set}{either a \code{ExpressionSet} or a \code{RangedSummarizedExperiment}}

\item{variable_names}{grouping variable used to perform differential expression analysis}

\item{covariable_names}{name of variables used in the adjusted models}

\item{type}{...}

\item{contrasts}{...}

\item{levels}{...}

\item{coef}{...}

\item{sva}{should differential expression analysis be adjusted by SVA?}

\item{annotCols}{variables from the annotation data used in the output}

\item{method}{String indicating the method used in the regression (e.g. lmFit function of limma: "ls" or 
"robust". (Default: "ls")}

\item{robust}{Logical value indicating whether robust method is applied in the eBayes function of limma. Default is FALSE.}

\item{normalization}{String indicating the normalize method used when using voom for RNAseq data (see normalized.method argument in limma::vomm)}

\item{voomQualityWeights}{Logical value indicating whether limma::voomWithQualityWeights should be used instead of
limma::voom.}

\item{big}{Logical value indicating whether SmartSVA should be used instead of SVA 
(TRUE recommended for methylation or when having large number of samples).}
}
\value{
a matrix with genes ordered by p-value
}
\description{
Performs differential expression analysis using LIMMA
}
\author{
Gonzalez, JR.
}
