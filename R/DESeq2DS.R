#'
#' @title Differential expression analysis on the server-side
#' @description Performs the DGE based on
#'  the Negative Binomial distribution 
#' @details This function perform a DGE analysis 
#' for \code{SummarizedExperiment} input 
#' @param vars variables for the model
#' @param set \code{SummarizedExperiment} object
#' @param test "Wald" or "LRT"
#' @param fitType "parametric", "local", or "mean"
#' @param sfType "ratio", "poscounts", or "iterate"
#' @param reduced test="LRT"
#' @param contrast model matrix contrast
#' 
#' @import dplyr
#' 
#' @author L. Abarrategui for DataSHILED development team
#' @export
#'

DESeq2DS<-function(vars, set,test, fitType, sfType, reduced, contrast)
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
    dds<-DESeq2::DESeq(dds,test = test, fitType = fitType, sfType = sfType)  
  }else{
    reduced <- unlist(strsplit(reduced, split=","))
    reduced <- paste("~", paste(c(reduced), collapse="+"))
    dds <- DESeq2::DESeq(dds,test = test, fitType = fitType,
                         sfType = sfType, reduced = stats::formula(reduced)) 
  }
   
  if(is.null(contrast)){
    res <- DESeq2::results(dds) 
  } else{
    contrast<-unlist(strsplit(contrast, split=","))
    res <- DESeq2::results(dds, contrast = contrast) 
  }
# Gene names (rownames) to column to avoid loosing them when converting to tibble
  res <- tibble::rownames_to_column(data.frame(res), "gene")
  return(as_tibble(res)[1,])
}

#AGGREGATE FUNCTION
#DESeq2DS
