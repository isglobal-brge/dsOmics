#'
#' @title Differential expression analysis on the server-side
#' @description Performs a differential expression analysis based on
#'  the negative binomial distribution. 
#' @details This function performs a differential expression analysis 
#' for \code{SummarizedExperiment} input using functions from 
#' \code{edgeR} Bioconductor's package. 
#' 
#' 
#' \code{edgeRDS} is a server-side aggregate function called by \code{ds.edgeR}.
#' 
#' @param set \code{RangedSummarizedExperiment} object
#' @param variable_names grouping variables used to perform differential expression analysis
#' @param intercept model intercept 
#' @param dispersion type of dispersion to estimate
#' @param normalization normalization method 
#' @param contrast character or numeric vector specifying the contrast
#' @param levels character or factor vector specifying the names of the contrast parameters 
#' @param test test for differential expression analysis 
#' @param coef a character or an integer indicating the coefficient to be tested equal to zero 
#' 
#' @author L. Abarrategui for DataSHILED development team
#' @export
#'

edgeRDS<-function(set, variable_names, intercept, dispersion, normalization, contrast, levels, test, coef)
{ 
  set<-eval(parse(text=set))
  
  counts<-SummarizedExperiment::assays(set)$counts
  pheno <- SummarizedExperiment::colData(set)
  
  #Group variable
  group<- paste("set",variable_names[1],sep = "$")
  group<-eval(parse(text=group))
  
  #design formula
  variable_names <- unlist(strsplit(variable_names, split=","))
  ff <- paste("~", paste(c(intercept,variable_names), collapse="+"))
  design <- model.matrix(stats::formula(ff), data = pheno)

  #Create DGEList object
  DGEList.object<-edgeR::DGEList(counts = counts, 
                                 samples= pheno, group = group)
  
  #Filtering
  keep <- edgeR::filterByExpr(DGEList.object)
  DGEList.object<- DGEList.object[keep,,keep.lib.sizes=FALSE]
  
  #Normalization
  DGEList.object <- edgeR::calcNormFactors(object = DGEList.object,
                                           method = normalization)
  
  
  
  #estimate the dispersion
  if(dispersion==1)
  {
   #common and tagwise dispersions
   DGEList.object <- edgeR::estimateDisp(y = DGEList.object,
                                         design = design)
  }
  
  if(dispersion==2 | dispersion==3)
  {
    #common dispersions
    DGEList.object  <- edgeR::estimateCommonDisp(y = DGEList.object,
                                                 design = design)
    
  }
  if(dispersion==3)
  {
    #tagwise dispersions
    DGEList.object <- edgeR::estimateTagwiseDisp(y = DGEList.object,
                                                 design = design)
    
    
  }
  
  
  
  if(!is.null(contrast)) 
  { 
    contrast <- unlist(strsplit(contrast, split=","))
    contrast.numeric<-as.numeric(contrast)
    if(is.null(contrast.numeric))
      {
       if(levels != "design")
         {
          levels <- unlist(strsplit(levels, split=",")) 
          colnames(design)<-levels 
          }
    contrast <-limma::makeContrasts(contrasts = contrast,levels = levels)
       }else{
      contrast<-contrast.numeric
    }
  }
  
  
  
  if(test==1)
  #Quasi-likelihood F-test
  {
    fit <- edgeR::glmQLFit(DGEList.object,design)
    if(is.null(contrast))
    {
    fit <- edgeR::glmQLFTest(fit, coef = coef)
    }else{
    fit <- edgeR::glmQLFTest(fit,coef = coef, contrast = contrast)
    }
  }
  
 if (test==2)
 #Likelihood ratio test
  {
   fit <- edgeR::glmFit(DGEList.object,design)
   if (is.null(contrast)){
     fit <- edgeR::glmLRT(fit,coef=coef) 
   }else{
     fit <- edgeR::glmLRT(fit,coef=coef,contrast = contrast) 
   }
    
  }
  
  #Result table
  results<-edgeR::topTags(fit)
  results<-results$table
  
  return(results)
}

#AGGREGATE FUNCTION
#edgeRDS