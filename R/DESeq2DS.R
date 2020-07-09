#'
#' @title Differential expression analysis on the server-side
#' @description Performs the Differential expression analysis based on the Negative Binomial distribution on the server-side
#' @details This function perform a differential expression analysis 
#' for \code{Count matrix}, \code{htseq-count} and \code{SummarizedExperiment} input 
#' @param vars variables for the model
#' @param set The data file 
#' @param test "Wald" or "LRT"
#' @param fitType "parametric", "local", or "mean"
#' @param sfType "ratio", "poscounts", or "iterate"
#' @param betaPrior TRUE or FALSE
#' @param reduced test="LRT"
#' @param contrast model matrix contrast
#' 
#' @author L. Abarrategui for DataSHILED development team
#' @export
#'

DESeq2DS<-function(vars, set,test, fitType, sfType, betaPrior, reduced, contrast)
{  
  set<-eval(parse(text=set))

   #design formula
   vars <- unlist(strsplit(vars, split=","))
   ff <- paste("~", paste(c(vars), collapse="+"))
   
   #convert data into DESeqDataSet
     dds <- DESeq2::DESeqDataSet(se = set, 
                                 design = stats::formula(ff))
     

  

  if(is.null(reduced))
    {
    dds<-DESeq2::DESeq(dds,test = test, fitType = fitType, sfType = sfType, betaPrior = betaPrior)  
  }else{
    reduced <- unlist(strsplit(reduced, split=","))
    reduced <- paste("~", paste(c(reduced), collapse="+"))
    dds<-DESeq2::DESeq(dds,test = test, fitType = fitType, sfType = sfType, betaPrior = betaPrior, reduced = stats::formula(reduced)) 
  }
   
  if(is.null(contrast)){
  res<-DESeq2::results(dds) 
  }else{
  contrast<-unlist(strsplit(contrast, split=","))
  res<-DESeq2::results(dds, contrast = contrast) 
  }
    

  return(res)
}

#AGGREGATE FUNCTION
# DESeq2DS