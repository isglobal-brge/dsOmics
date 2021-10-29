#' Title
#'
#' @param x 
#' @param replacement 
#' @param byRow 
#'
#' @return
#' @export
#'
#' @examples
replace.NA <- function(x, replacement , byRow = TRUE){ # Extracted from the MCRestimate package to avoid the dependency
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

#' Title
#'
#' @param x 
#' @param y 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
fastGWAS_getFitted.values <- function(x, mod_names, output_family, ...){
  # TODO ara aixo nomes funcione per una variable objectiu i una covariable,
  # ha de poder funcionar amb N covariables i sense covariables!
  # TODO tambe ha de funcionar per gaussian i binomial1!!!!!
  mod_values <- as.numeric(unlist(list(...)))
  if(output_family == "gaussian"){
    fitted.values <- .fittedValues_linear(mod_values, mod_names, x)
  } else if (output_family == "binomial") {
    fitted.values <- .fittedValues_exponential(mod_values, mod_names, x)
  } else {
    stop()
  }
  return(unlist(fitted.values))
}

#' Title
#'
#' @param mod_values 
#' @param mod_values_names 
#' @param x 
#' @param intercept 
#'
#' @return
#' @export
#'
#' @examples
.fittedValues_linear <- function(mod_values, mod_values_names, x){
  fitted.values <- mod_values[1] + Reduce("+", lapply(1:length(mod_values[-1]), function(i){
    mod_values[-1][i]*x[mod_values_names[i]]}))
}

#' Title
#'
#' @param mod_values 
#' @param mod_values_names 
#' @param x 
#' @param intercept 
#'
#' @return
#' @export
#'
#' @examples
.fittedValues_exponential <- function(mod_values, mod_values_names, x){
  fitted.values <- (exp(mod_values[1] + 
                          Reduce("+",lapply(1:length(mod_values[-1]), function(i){
                            mod_values[-1][i]*x[mod_values_names[i]]
                          }))))/(1+exp(mod_values[1] + Reduce("+",lapply(1:length(mod_values[-1]), function(i){
                            mod_values[-1][i]*x[mod_values_names[i]]
                          }))))
}

#' Title
#'
#' @param x 
#' @param y 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
fastGWAS_getResiduals <- function(x, fitted.values, objective_variable, output_family){
  # TODO fer que funcioni tambe per gaussian, ara nomes esta per binomial
  # Working residuals
  if(output_family == "gaussian"){
    return(unlist(x[objective_variable] - fitted.values))
  } else if (output_family == "binomial") {
    return(unlist((x[objective_variable]-fitted.values) / (fitted.values*(1-fitted.values))))
  } else {
    stop()
  }
  
}

#' Title
#'
#' @param genoData 
#'
#' @return
#' @export
#'
#' @examples
fastGWAS_S <- function(genoData, type, snpBlock = NULL){
  if(type == "geno"){
    # Get genotype data
    geno <- lapply(genoData, function(x){
      GWASTools::GenotypeBlockIterator(x, snpBlock=snpBlock)
    })
    return(geno)
    # geno <- lapply(genoData, function(x){
    #   ids_x <- GWASTools::getVariable(x, "snp.rs.id")
    #   geno_x <- t(data.frame(GWASTools::getGenotype(x)))
    #   sample_id_x <- GWASTools::getVariable(x, "sample.id")
    #   colnames(geno_x) <- ids_x
    #   rownames(geno_x) <- sample_id_x
    #   return(geno_x)
    # })
    # result <- Reduce(cbind, geno)
    # class(result) <- c(class(result), "GenoInfoTable")
    # return(result)
  } else if (type == "pheno") {
    phenotypes <- lapply(genoData, function(x){
      vv <- GWASTools::getScanVariableNames(x)
      GWASTools::getScanVariable(x, vv)
    })
    if(!all(sapply(phenotypes[2:length(phenotypes)], FUN = identical, phenotypes[[1]]))){
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

#' Title
#'
#' @return
#' @export
#'
#' @examples
fastGWAS_ColSums <- function(table1, table2, type, means = NULL, pheno = NULL, geno = NULL){
  # TODO check que el genoinfotable sigui del tipus que toque!!! sino el colsums altanto
  # TODO ficar potser algun filtre de minim delements (rows) a fer colsum
  
  # TODO check que el type sigui o "std" o "square" o "crossprod"
  if(type == "std"){
    colSums(table1)
  } else if(type == "square"){
    colSums(table1 ^ 2)
  } else if (type == "crossprod") {
    
    ids <- pheno$scanID
    sample.index <- which(getScanID(geno[[1]]) %in% ids)
    results <- Reduce(c, lapply(geno, function(x){
      GWASTools::resetIterator(x)
      crossprod <- colSums(table1 * t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>% 
                                        replace_na(as.list(means))))
      while(GWASTools::iterateFilter(x)){
        crossprod <- c(crossprod, colSums(table1 * t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>%
                                                       replace_na(as.list(means)))))
      }
      return(crossprod)
    }))
    return(results)
  } else if (type == "std_vect") {
    sum(table1^2)
  } else if (type == "square_vect") {
    sum(table1^2)
  } else if (type == "std_geno") {
    # Table 1 contains the GenoBlockIterator list
    # Table 2 contains the pheno data, grab the individuals 
    # (it has been filtered previously with fastGWAS_PHENO_removeNAindiv)
    ids <- table2$scanID
    sample.index <- which(getScanID(table1[[1]]) %in% ids)
    results <- Reduce(c, lapply(table1, function(x){
      GWASTools::resetIterator(x)
      sums <- colSums(t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>% replace_na(as.list(means))))
      while(GWASTools::iterateFilter(x)){
        sums <- c(sums, colSums(t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>% replace_na(as.list(means)))))
      }
      return(sums)
    }))
    return(results)
  } else if (type == "square_geno") {
    # Table 1 contains the GenoBlockIterator list
    # Table 2 contains the pheno data, grab the individuals 
    # (it has been filtered previously with fastGWAS_PHENO_removeNAindiv)
    ids <- table2$scanID
    sample.index <- which(getScanID(table1[[1]]) %in% ids)
    results <- Reduce(c, lapply(table1, function(x){
      GWASTools::resetIterator(x)
      sums <- colSums(t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>% replace_na(as.list(means))) ^ 2)
      while(GWASTools::iterateFilter(x)){
        sums <- c(sums, colSums(t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>% replace_na(as.list(means))) ^ 2))
      }
      return(sums)
    }))
    return(results)
  } else if (type == "GWASsums") {
    ids <- pheno$scanID
    sample.index <- which(getScanID(geno[[1]]) %in% ids)
    
    # n.cores <- parallel::detectCores() - 1
    # my.cluster <- parallel::makeCluster(
    #   n.cores, 
    #   type = "PSOCK"
    # )
    # doParallel::registerDoParallel(cl = my.cluster)
    # foreach::getDoParRegistered()
    # foreach::getDoParWorkers()
    # 
    # results <- Reduce(c, parLapply(my.cluster, geno, function(i){
    #   node <- GWASTools::GdsGenotypeReader(i@data@filename, allow.fork = TRUE)
    #   
    #   geno_data <- t(replace.NA(GWASTools::getGenotype(node)[,sample.index], means, F))
    #   
    #   crossprod <- Rfast::colsums(table1 * geno_data)
    #   sums <- Rfast::colsums(geno_data)
    #   squared_sums <- Rfast::colsums(geno_data ^ 2)
    #   
    #   return(list(crossprod = crossprod, sums = sums, squared_sums = squared_sums))
    # }))
    # 
    # stopCluster(my.cluster)
    
      results <- Reduce(c, lapply(geno, function(x){
      GWASTools::resetIterator(x)
      # genoData_temp <- t(as.tibble(getGenotypeSelection(x, scan=sample.index, order="selection")) %>% replace_na(as.list(means)))
      genoData_temp <- t(replace.NA(getGenotypeSelection(x, scan=sample.index, order="selection"), means, F))
      # crossprod <- Rfast::colsums(table1 * genoData_temp)
      crossprod <- colSums(table1 * genoData_temp)
      # sums <- Rfast::colsums(genoData_temp)
      sums <- colSums(genoData_temp)
      # squared_sums <- Rfast::colsums(genoData_temp ^ 2)
      squared_sums <- colSums(genoData_temp ^ 2)

      while(GWASTools::iterateFilter(x)){
        genoData_temp <- t(replace.NA(getGenotypeSelection(x, scan=sample.index, order="selection"), means, F))
        # crossprod <- c(crossprod, Rfast::colsums(table1 * genoData_temp))
        crossprod <- c(crossprod, colSums(table1 * genoData_temp))
        # sums <- c(sums, Rfast::colsums(genoData_temp))
        sums <- c(sums, colSums(genoData_temp))
        # squared_sums <- c(squared_sums, Rfast::colsums(genoData_temp ^ 2))
        squared_sums <- c(squared_sums, colSums(genoData_temp ^ 2))
      }

      return(list(crossprod = crossprod, sums = sums, squared_sums = squared_sums))
    }))
        
    results_combined <- tapply(unlist(results, use.names = FALSE), rep(names(results), lengths(results)), FUN = c)
    return(results_combined)
  } else {
    stop()
  }
  
}

#' Title
#'
#' @param pheno 
#' @param objective 
#' @param covars 
#'
#' @return
#' @export
#'
#' @examples
fastGWAS_PHENO_removeNAindiv <- function(pheno, objective, covars){
  # TODO ficar un inherits o algo aixi per assegurarnos que el pheno es de classe PhenoInfoTable
  return(pheno[complete.cases(pheno[, c(objective, covars)]), c(objective, covars, "scanID")])
}

#' Title
#'
#' @param geno 
#' @param pheno 
#'
#' @return
#' @export
#'
#' @examples
fastGWAS_S_impute <- function(geno, pheno){
  
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
  output <- rowSums(results) / n_snps
  return(output)

# 
#   # TODO ficar un inherits o algo aixi per assegurarnos que el pheno es de classe PhenoInfoTable i GenoInfoTable
#   # Remove IDs from geno that do not have the objective/covariables
#   browse
#   geno <- geno[rownames(geno) %in% pheno$scanID,]
#   # Impute geno with the sample mean of the observed genotypes
#   k <- which(is.na(geno), arr.ind=TRUE)
#   geno[k] <- rowMeans(geno, na.rm=TRUE)[k[,1]]
#   return(geno)
}
