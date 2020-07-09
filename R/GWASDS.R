#' @title Genome-wide association analysis (GWAS)
#' 
#' @description Performs GWAS using GENESIS
#' @param genoData a \code{GenotypeData} object which is a container for storing genotype data
#' from a GWAS toghether with the metadata associated with the subjects (i.e. phenotypes and/or covariates)
#' and SNPs
#' @param outcome A character string specifying the name of the outcome variable in \code{genoData}
#' @param covars A vector of character strings specifying the names of the fixed effect covariates 
#' in \code{genoData}; an intercept term is automatically included. If NULL (default) the only fixed effect 
#' covariate is the intercept term
#' @param family A description of the generalized linear model used. The defatul is "binomial" that is defined
#' for case/control studies. Quantitative traits can be analyzed by using "gaussian". Other values are accepted.
#' @param snpBlock an integer specifying the number of SNPs in an iteration block. See \code{GenotypeIterator} 
#' function in GWASTools package.
#'  
#' @param ... other arguments of \code{fitNullModel} function in GENESIS package
#' @return a matrix with SNPs ordered by p-values
#' 
#' @author Gonzalez, JR.
#'
#' @export 
#' 
GWASDS <- function(genoData, outcome, covars=NULL, family="binomial", snpBlock, ...){
  covars <- unlist(strsplit(covars, split=","))
  nullmod <- GENESIS::fitNullModel(genoData, outcome = outcome, 
                                   covars = covars, 
                                   family = family)
  genoIterator <- GWASTools::GenotypeBlockIterator(genoData, snpBlock=snpBlock)
  assoc <- GENESIS::assocTestSingle(genoIterator, null.model = nullmod)
  assoc$rs<-GWASTools::getVariable(genoData, "snp.rs.id")[assoc$variant.id]
  ans <- assoc %>% as_tibble() %>%
    arrange(Score.pval) %>% select(variant.id, rs, everything())
  
  return(ans)
}
