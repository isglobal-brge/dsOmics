##' Performing a linear regression analysis on pooled data from multiple studies for every CpG site
##' 
##' Function fits a generalized linear model to methylation data for each feature (CpG site) 
##' in the data sets considered, using user specified phenotypic variables as covariates
##' Outputs a matrix containing a beta value, standard error and p-value for each CpG site
##' @title Linear regression analysis of pooled data for each CpG site in study
##' @param cpgs an optional parameter input as a vector of integer values which indicates the indices of specific 
##' CpG sites that should be analysed. By default all CpG sites are analysed
##' @param model list of phenotypic variables to use as covariates in the regression analysis in the form: 
##' "phenoVar1 ~ phenoVar2 + phenoVar3 ... + phenoVarN"
##' @param molecular.data name of the DataSHIELD object to which the methylation data set has been assigned
##' @param pheno.data name of the DataSHIELD object to which the phenotype data set has been assigned
##' @param datasources Opal object or list of opal objects denoting the opal server(s) information
##' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use during 
##' parallel processing. Argument can only be > 1 when the function is run on a linux machine
##' @param type.p.adj multiple comparison correction method. Default 'fdr' 
##' @param cellCountsAdjust optional TRUE/FALSE parameter which indicates whether or not the linear regression
##' models should be adjusted for the estimated cell counts by including the variables in the models.
##' NOTE: This assumes that the Opal pheno tables for every study include the necessary estimated cell count data 
##' originally computed when running the createOpalFiles function
##' @export
##' @examples
##' 

ds.lmFeature <- function(cpgs=NULL, model, molecular.data, pheno.data, datasources, mc.cores = 1, 
                         type.p.adj='fdr', cellCountsAdjust = FALSE){
  
  cellCountWarning = FALSE
  
  #Adding cell count variables to model if cellCountsAdjust argument has been specified
  if(isTRUE(cellCountsAdjust)){
    
    phenoVars = ds.colnames(pheno.data)
    
    #Checking that the cell counts have been computed for all studies
    if(length(phenoVars) > 1){
      for(i in 1:(length(phenoVars)-1)){
        if(!identical(grep("_Scaled", phenoVars[[i]], value = TRUE), 
                      grep("_Scaled", phenoVars[[i+1]], value = TRUE))){
          return(message(paste(" ", "ERROR: Cell count variables do not match between one or more of the studies or are missing.",
                               "       Set 'cellCountsAdjust' = FALSE and re-run analysis to avoid error.",
                        " ", sep="\n")))
        }
      }
    }
    
    #Adding cell count variables as covariates to model
    cellCountsVars = grep("_Scaled", phenoVars[[1]], value = TRUE)
    if(length(cellCountsVars) > 0){
      model = paste(model, paste0(cellCountsVars, collapse = "+"), sep="+")
    }else{
      cellCountWarning = TRUE
    }
    
  }
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  #Verifying that the number of features is the same for all studies
  studyCPGnums = ds.dim(molecular.data)
  if(length(studyCPGnums) > 1){
    for(i in 1:(length(studyCPGnums)-1)){
      if(!(studyCPGnums[[i]][1] == studyCPGnums[[i+1]][1])){
        return(message(paste(" ", "ERROR: The number of features (cpg sites) do not match across all studies",
                             " ", sep="\n")))
      }
    }
  }
  
  #Setting the number of cpgs to loop over as the total number of 
  #features in the studies if no indices are specified
  if(is.null(cpgs)){
    cpgs = 1:studyCPGnums[[1]][1]
  }
  
  
  lmFeatureGLM <- function(k, cpgs, model, molecular.data, pheno.data, datasources, vars){
    
    #Creating and executing function call which returns all methylation data related 
    #to cpg site cpgs[k] and the phenotypic data corresponding to variables specified in 'model' formula
    cally <- paste0("lmFeatureDS(cpgIndex=", cpgs[k], ", ", "model=", model, ", " ,"molecular.data=", molecular.data, ", ", "pheno.data=", pheno.data,")")
    datashield.assign(datasources, 'data.sel', as.symbol(cally))
    
    dep <- paste(vars, collapse="+")
    colnames <- ds.colnames('data.sel')
    feature <- colnames[[1]][1]
    mm <- as.formula(paste(feature, "~", dep))
    
    mod <- ds.glm.o(mm, family='gaussian', data='data.sel')
    
    metrics = as.data.frame(mod$coefficients[2, c(1,2,4)])
    names(metrics) = feature
    
    return(metrics)
    
  }
  
  ans = t(as.data.frame(mclapply(1:length(cpgs), lmFeatureGLM, 
                                 cpgs, model, molecular.data, pheno.data, datasources, vars, 
                                 mc.cores = mc.cores))) 
  
  if(isTRUE(cellCountWarning)){
    message(paste(" ", "WARNING: Cell count data not found in any of the studies.",
                  "         Linear regression analysis performed without cell count variables.",
                  " ", sep="\n"))
  }
  
  colnames(ans) <- c("beta", "s.e.", "p.value")
  
  #Ordering matrix by ascending p.value
  ans = ans[order(ans[, 'p.value']), ]
  
  #Calculating adjusted p-values
  ans = cbind(ans, p.adj = p.adjust(ans[,'p.value'],  method=type.p.adj))
  
  class(ans) <- c("dsMethy", class(ans))
  
  return(ans)
  
}