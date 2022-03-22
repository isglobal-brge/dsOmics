#' @title Extract Geno iterator and Pheno information from GenotypeData GDS containers
#'
#' @param genoData \code{GenotypeData / vector of GenotypeData} container for storing genotype data
#' from a GWAS toghether with the metadata associated with the subjects (i.e. phenotypes and/or covariates)
#' @param type \code{character} What to extract from the genoData, the GenoBlockIterators [geno] 
#' or the phenotypes [pheno]
#' @param snpBlock \code{numeric} Block size for dividing the genotype data, it equals to the 
#' number of SNPs used on each iteration, depending on the servers RAM it may perform faster using lower or greater 
#' block sizes, do some testing to assess it.
#'
#' @return \code{PhenoInfoTable} or \code{GenotypeBlockIterator}
#' @export
#'

fastGWAS_S <- function(genoData, type, snpBlock){
  if(type == "geno"){
    geno <- lapply(genoData, function(x){
      GWASTools::GenotypeBlockIterator(x, snpBlock=snpBlock)
    })
    return(geno)
  } else if (type == "pheno") {
    phenotypes <- lapply(genoData, function(x){
      vv <- GWASTools::getScanVariableNames(x)
      GWASTools::getScanVariable(x, vv)
    })
    if(length(phenotypes) > 1 & !all(sapply(phenotypes[2:length(phenotypes)], FUN = identical, phenotypes[[1]]))){
      stop('Provided genoData objects have different phenotypes associated, make sure
         when creating them using "ds.GenotypeData" they match.')
    } else {
      result <- data.frame(phenotypes[[1]])
      class(result) <- c(class(result), "PhenoInfoTable")
      return(result)
    }
  } else {
    stop('Wrong type, options are ["geno", "pheno"]')
  }
}

#' @title Remove individuals from the phenotype that have NAs on any of the objective variable or covars
#'
#' @param pheno \code{PhenoInfoTable} to which remove individuals
#' @param objective \code{character} Objective variable
#' @param covars \code{character vector} Covariables
#'
#' @return \code{PhenoInfoTable} with individuals trimmed
#' @export
#'

fastGWAS_PHENO_removeNAindiv <- function(pheno, objective, covars){
  if(!inherits(pheno, "PhenoInfoTable")){
    stop('[pheno] argument shoud be of class "PhenoInfoTable", generated using dsOmics::fastGWAS_S')
  }
  return(pheno[complete.cases(pheno[, c(objective, covars)]), c(objective, covars, "scanID")])
}

#' @title Get mean of the genotype by individual
#'
#' @param geno \code{GenotypeBlockIterator} Genotype information
#' @param pheno \code{PhenoInfoTable} Phenotype information
#'
#' @return \code{named numeric vector} with the means by individual. The names are the individual IDs
#' @export
#'

fastGWAS_S_means <- function(geno, pheno){
  if(!all(unlist(lapply(geno, function(x){inherits(x, "GenotypeBlockIterator")})))){
    stop('[geno] argument shoud be of class "GenotypeBlockIterator", generated using dsOmics::fastGWAS_S')
  }
  if(!inherits(pheno, "PhenoInfoTable")){
    stop('[pheno] argument shoud be of class "PhenoInfoTable", generated using dsOmics::fastGWAS_S')
  }
  ids <- pheno$scanID
  sample.index <- which(getScanID(geno[[1]]) %in% ids)
  results <- Reduce(cbind, lapply(geno, function(x){
    GWASTools::resetIterator(x)
    sums <- rowSums(getGenotypeSelection(x, scan=sample.index, order="selection",transpose=TRUE), na.rm = T)
    while(GWASTools::iterateFilter(x)){
      sums <- cbind(sums, rowSums(getGenotypeSelection(x, scan=sample.index, order="selection",transpose=TRUE), na.rm = T))
    }
    return(sums)
  }))
  
  n_snps <- Reduce("+", lapply(geno, function(x){
    tail(x@snpFilter[[length(x@snpFilter)]], 1)
  }))
  # diff P
  #############################################################
  # CAPTURE THE diffP SETTINGS
  nfilter.diffP.epsilon <- getOption("default.nfilter.diffP.epsilon")
  #############################################################
  if(!is.null(nfilter.diffP.epsilon)){
    # l1-sensitivity = 2 (geno encoding 0,1,2; when performing rowSums max difference when extracting
    # one individual is 2)
    laplace_noise <- Laplace_noise_generator(m = 0, 
                            b = 2/nfilter.diffP.epsilon, 
                            n.noise = length(results))
    results <- results + laplace_noise
  }
  output <- rowSums(results) / n_snps
  return(output)
}

#' @title Obtain fitted values given coefficients, family and values table
#'
#' @param x \code{PhenoInfoTable} Pheno data to get the fitted values from
#' @param covars \code{character vector} Covariables
#' @param mod_names \code{character vector} Names of the mod_values
#' @param output_family \code{character} A description of the generalized linear model used. "binomial" is defined
#' for case/control studies. Quantitative traits can be analyzed by using "gaussian"
#' @param mod_values \code{numeric vector} Values of model 0
#'
#' @return \code{numeric vector} of fitted values
#' @export
#'

fastGWAS_getFitted.values <- function(x, covars, mod_names, output_family, mod_values){
  if(!inherits(x, "PhenoInfoTable")){
    stop('[x] argument shoud be of class "PhenoInfoTable", generated using dsOmics::fastGWAS_S')
  }
  # Check with variables are to be recoded into dummies
  covars_to_recode <- covars[!(covars %in% mod_names)]
  
  # Recode dummies
  x_dummies <- dummies::dummy.data.frame(x, names = covars_to_recode, sep = "")
  
  if(output_family == "gaussian"){
    fitted.values <- .fittedValues_linear(mod_values, mod_names, x_dummies)
  } else if (output_family == "binomial") {
    fitted.values <- .fittedValues_exponential(mod_values, mod_names, x_dummies)
  } else {
    stop()
  }
  return(unlist(fitted.values))
}

#' @title Fit values for linear model
#' @description Internal function

.fittedValues_linear <- function(mod_values, mod_values_names, x){
  fitted.values <- mod_values[1] + Reduce("+", lapply(1:length(mod_values[-1]), function(i){
    mod_values[-1][i]*x[mod_values_names[i]]}))
}

#' @title Fit values for binomial model
#' @description Internal function

.fittedValues_exponential <- function(mod_values, mod_values_names, x){
  fitted.values <- (exp(mod_values[1] + 
                          Reduce("+",lapply(1:length(mod_values[-1]), function(i){
                            mod_values[-1][i]*x[mod_values_names[i]]
                          }))))/(1+exp(mod_values[1] + Reduce("+",lapply(1:length(mod_values[-1]), function(i){
                            mod_values[-1][i]*x[mod_values_names[i]]
                          }))))
}

#' @title Obtain residual values given fitted values, family and objective variable
#'
#' @param x \code{PhenoInfoTable} Pheno data to get the fitted values from
#' @param fitted.values \code{numeric vector} Fitted values
#' @param objective_variable \code{character} Objective variable
#' @param output_family \code{character} A description of the generalized linear model used. "binomial" is defined
#' for case/control studies. Quantitative traits can be analyzed by using "gaussian"
#'
#' @return \code{numeric vectors} of residual values
#' @export
#'

fastGWAS_getResiduals <- function(x, fitted.values, objective_variable, output_family){
  if(!inherits(x, "PhenoInfoTable")){
    stop('[x] argument shoud be of class "PhenoInfoTable", generated using dsOmics::fastGWAS_S')
  }
  # TODO fer que funcioni tambe per gaussian, ara nomes esta per binomial
  # Working residuals
  if(output_family == "gaussian"){
    return(unlist(x[objective_variable] - fitted.values))
  } else if (output_family == "binomial") {
    return(unlist((x[objective_variable]-fitted.values) / (fitted.values*(1-fitted.values))))
  } else {
    stop('Invalid family argument')
  }
  
}

#' @title Get colsums needed for fast GWAS
#'
#' @param table1 \code{numeric vector} To be squared or multiplied to genoData, depending on type
#' @param type \code{character} ["square_vect"] to square and sum table1 or ["GWASsums"] to get all the colsums 
#' needed for fast GWAS
#' @param do.par \code{bool} Whether to use parallelization on the servers, to do so the servers 
#' have to have the package \code{doParallel} installed and run on a POSIX OS (Mac, Linux, Unix, BSD); Windows 
#' is not supported. This parallelization computes in parallel each \code{genoData} object, therefore it is only useful 
#' when the genoData is divided by chromosome. 
#' @param n.cores \code{numeric} Numbers of cores to use when \code{do.par} is \code{TRUE}. If 
#' \code{NULL} the number of cores used will be the maximum available minus one.
#' @param means \code{numeric vector} of geno means by individual
#' @param pheno \code{PhenoInfoTable} Pheno data
#' @param geno \code{GenotypeBlockIterator} Genotype information
#'
#' @return \code{numeric vectors}
#' @export
#'

fastGWAS_ColSums <- function(table1, type, do.par = FALSE, n.cores = NULL, means = NULL, pheno = NULL, geno = NULL){
  # TODO check que el genoinfotable sigui del tipus que toque!!! sino el colsums altanto
  # TODO ficar potser algun filtre de minim delements (rows) a fer colsum
  
  # TODO check que el type sigui o "std" o "square" o "crossprod"
  if (type == "square_vect") {
    sum(table1^2)
  } else if (type == "GWASsums") {
    ids <- pheno$scanID
    sample.index <- which(getScanID(geno[[1]]) %in% ids)
    
    # if(do.par){
    #   if(is.null(n.cores)){n.cores <- parallel::detectCores() - 1}
    #   results <- Reduce(c, parallel::mclapply(geno, function(x){
    #     GWASTools::resetIterator(x)
    #     genoData_temp <- t(replace.NA(getGenotypeSelection(x, scan=sample.index, order="selection"), means, F))
    #     crossprod <- colSums(table1 * genoData_temp)
    #     sums <- colSums(genoData_temp)
    #     squared_sums <- colSums(genoData_temp ^ 2)
    #     
    #     while(GWASTools::iterateFilter(x)){
    #       genoData_temp <- t(replace.NA(getGenotypeSelection(x, scan=sample.index, order="selection"), means, F))
    #       crossprod <- c(crossprod, colSums(table1 * genoData_temp))
    #       sums <- c(sums, colSums(genoData_temp))
    #       squared_sums <- c(squared_sums, colSums(genoData_temp ^ 2))
    #     }
    #     
    #     return(list(crossprod = crossprod, sums = sums, squared_sums = squared_sums))
    #   }, mc.cores = n.cores))
    # } else {
    
    # maxes <- Reduce(cbind,lapply(geno, function(x){
    #   geno_temp <- GWASTools::getVariable(x, "genotype")
    #   return(max(geno_temp))
    # }))
    # browser()
      results <- Reduce(c, lapply(geno, function(x){
        GWASTools::resetIterator(x)
        genoData_temp <- t(replace.NA(getGenotypeSelection(x, scan=sample.index, order="selection"), means, F))
        crossprod <- colSums(table1 * genoData_temp)
        sums <- colSums(genoData_temp)
        squared_sums <- colSums(genoData_temp ^ 2)
        geno_maxs <- matrixStats::colMaxs(genoData_temp)
        
        while(GWASTools::iterateFilter(x)){
          genoData_temp <- t(replace.NA(getGenotypeSelection(x, scan=sample.index, order="selection"), means, F))
          crossprod <- c(crossprod, colSums(table1 * genoData_temp))
          sums <- c(sums, colSums(genoData_temp))
          squared_sums <- c(squared_sums, colSums(genoData_temp ^ 2))
          geno_maxs <- c(geno_maxs, matrixStats::colMaxs(genoData_temp))
        }
        
        return(list(crossprod = crossprod, sums = sums, squared_sums = squared_sums, geno_maxs = geno_maxs))
      }))
    # }
    results_combined <- tapply(unlist(results, use.names = FALSE), rep(names(results), lengths(results)), FUN = c)
    # diff P
    #############################################################
    # CAPTURE THE diffP SETTINGS
    nfilter.diffP.epsilon <- getOption("default.nfilter.diffP.epsilon")
    #############################################################
    if(!is.null(nfilter.diffP.epsilon)){
      # l1-sensitivity = 2 (geno encoding 0,1,2; when performing colSums max difference when extracting
      # one individual is 2)
        # Sums
      results_combined$geno_maxs[which(results_combined$geno_maxs == 0)] <- .Machine$double.xmin
      laplace_noise <- Laplace_noise_generator(m = 0,
                                               b = results_combined$geno_maxs/nfilter.diffP.epsilon,
                                               n.noise = length(results_combined$sums))
      results_combined$sums <- results_combined$sums + laplace_noise
      # l1-sensitivity = 2^2
        # Squares sums
      results_combined$squared_sums <- results_combined$squared_sums + laplace_noise^2
      # l1-sensitivity = max(table1) * 2
      l2.sens <- max(table1) * results_combined$geno_maxs
      laplace_noise <- Laplace_noise_generator(m = 0,
                                               b = l2.sens/nfilter.diffP.epsilon,
                                               n.noise = length(results_combined$sums))
      results_combined$crossprod <- results_combined$crossprod + laplace_noise
    }
    return(results_combined)
  } else {
    stop()
  }
  
}

##### Auxiliary function to fast replace NAs on matrix copied from the MCRestimate package to avoid the dependency
#' @export

replace.NA <- function(x, replacement , byRow = TRUE){
  if (byRow==TRUE){
    norows=nrow(x)
    lengthy=length(replacement)
    if (norows==lengthy) {
      for (i in 1:norows) {
        x[i,is.na(x[i,])]=replacement[i]
      }
    } else {
      print("Error: length of replacement vector does not match number of rows")
    }
  }
  else {
    nocol=ncol(x)
    lengthy=length(replacement)
    if (nocol==lengthy) {
      for (i in 1:nocol) {
        x[is.na(x[,i]),i]=replacement[i]
      }
    }
    else {
      print("Error: length of replacement vector does not match number of columns")
    }
  }
  return(x)
}