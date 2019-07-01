##' Performing a linear regression analysis on pooled data from multiple studies for every CpG site
##' 
##' Function fits a generalized linear model to methylation data for each feature (CpG site) 
##' in the data sets considered, using user specified phenotypic variables as covariates
##' Outputs a matrix containing a beta value, standard error and p-value for each CpG site
##' @title Linear regression analysis of pooled data for each CpG site in study
##' @param cpgs an optional parameter input as a vector of integer values which indicates the indices of specific 
##' CpG sites that should be analysed. By default all CpG sites are analysed
##' @param model list of phenotypic variables to use as covariates in the regression analysis in the form: 
##' "phenoVar1 ~ phenoVar2 ~ phenoVar3 ... ~ phenoVarN"
##' @param molecular.data name of the DataSHIELD object to which the methylation data set has been assigned
##' @param pheno.data name of the DataSHIELD object to which the phenotype data set has been assigned
##' @param datasources Opal object or list of opal objects denoting the opal server(s) information
##' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use during 
##' parallel processing. Argument can only be > 1 when the function is run on a linux machine
##' @export
##' @examples
##' 

ds.lmFeature <- function(cpgs=NULL, model, molecular.data, pheno.data, datasources, mc.cores = 1){
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  #Verifying that the number of features is the same for all studies
  studyCPGnums = ds.dim(molecular.data)
  for(i in 1:(length(studyCPGnums)-1)){
    if(!(studyCPGnums[[i]][1] == studyCPGnums[[i+1]][1])){
      return(message(paste(" ", "ERROR: The number of features (cpg sites) do not match across all studies",
                           " ", sep="\n")))
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
    cally <- paste0("lmFeatureDS(i=", cpgs[k], ", ", "model=", model, ", " ,"molecular.data=", molecular.data, ", ", "pheno.data=", pheno.data,")")
    datashield.assign(datasources, 'data.sel', as.symbol(cally))
    
    dep <- paste(vars, collapse="+")
    colnames <- ds.colnames('data.sel')
    feature <- colnames$study1[1]
    mm <- as.formula(paste(feature, "~", dep))
    
    mod <- ds.glm.o(mm, family='gaussian', data='data.sel')
    
    metrics = as.data.frame(mod$coefficients[2, c(1,2,4)])
    names(metrics) = feature
    
    return(metrics)
    
  }
  
  ans = t(as.data.frame(mclapply(1:length(cpgs), lmFeatureGLM, 
                                 cpgs, model, molecular.data, pheno.data, datasources, vars, 
                                 mc.cores = mc.cores))) 
  
  colnames(ans) <- c("beta", "s.e.", "p.value")
  
  #Ordering matrix by ascending p.value
  ans = ans[order(ans[, 'p.value']), ]
  
  class(ans) <- c("dsMethy", class(ans))
  
  return(ans)
  
}