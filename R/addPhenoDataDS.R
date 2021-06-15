#' @title Add Phenotype data to ExpressionSet
#'
#' @param x \code{ExpressionSet} ExpressionSet to which add phenotype information
#' @param pheno \code{data.frame} Table with the new phenotypes
#' @param identifier \code{character} Name of the ID column on the phenotypes data.frame
#' @param complete_cases \code{bool} If \code{TRUE} only the matching individuals 
#' between the ExpressionSet and the phenotypes table will be included on the resulting ExpressionSet. If 
#' \code{FALSE} all the individuals on the input ExpressionSet will be on the output ExpressionSet
#' @param alternate_eset_id \code{character} Alternate ID of the eSet pheno data, by default the rownames 
#' of the eSet pheno data act as ID, use this argument if the ID to merge the individuals is on a column of the pheno data. 
#' Input NULL for the standard behaviour of using the rownames of the pheno data as ID.
#'
#' @return
#' @export
#'
#' @examples
addPhenoDataDS <- function(x, pheno, identifier, alternate_eset_id, complete_cases){
  
  if(!(any(identifier %in% colnames(pheno)))){
    stop("Identifier [", identifier, "] is not on the phenotypes table")
  }
  
  og_pheno <- Biobase::pData(x)
  og_pheno_md <- Biobase::varMetadata(x)
  
  new_variables <- colnames(pheno)[!(identifier == colnames(pheno))]
  old_variables <- colnames(og_pheno)
  
  if(is.null(alternate_eset_id)){
    og_individuals <- rownames(og_pheno)
  } else {
    og_individuals <- as.character(og_pheno[,alternate_eset_id])
  }
  
  if(!is.null(alternate_eset_id) & length(unique(og_individuals)) != nrow(og_pheno)){
    stop('The selectected alternate_eset_id[', alternate_eset_id, '] does not correspond to a unique identifier, there are repeated IDs in this column')
  }
  
  new_individuals <- as.character(unlist(pheno[,identifier]))
  common_individuals <- new_individuals %in% og_individuals
  
  if(all(common_individuals == FALSE)){
    stop('No common individuals between the ExpressionSet and the new pheno data')
  }
  
  new_pheno <- pheno[common_individuals,]
  new_pheno[,identifier] <- as.character(unlist(new_pheno[,identifier]))
  og_pheno <- cbind(og_pheno, og_individuals_id = og_individuals)
  og_pheno$og_individuals_id <- as.character(og_pheno$og_individuals_id)
  
  if(complete_cases == TRUE){
    rownames_new_pheno <- rownames(og_pheno[og_individuals %in% new_individuals[common_individuals],])
  } else{
    stop('complete cases FALSE not implemented')
  }
  
  if(complete_cases == TRUE){
    new_pheno <- dplyr::right_join(og_pheno, new_pheno, by = c("og_individuals_id" = identifier))
    assay_data <- Biobase::exprs(x)[,colnames(Biobase::exprs(x)) %in% rownames_new_pheno]
  }
  else{
    stop('complete cases FALSE not implemented')
    # new_pheno <- dplyr::left_join(og_pheno, new_pheno, by = c("og_individuals_id" = identifier))
    # assay_data <- Biobase::exprs(x)
  }
  
  rownames(new_pheno) <- rownames_new_pheno
  new_pheno$og_individuals_id <- NULL
  
  if(any(new_variables %in% old_variables)){stop("Variables conflict between ExpressionSet and new PhenoData")}
  for(i in new_variables){
    og_pheno_md <- eval(str2expression(paste0("rbind(og_pheno_md, ", i, " = NA)")))
  }
  
  new_pheno <- new("AnnotatedDataFrame", data=new_pheno, varMetadata=og_pheno_md)
  
  eset <- Biobase::ExpressionSet(assayData = assay_data,
                                 phenoData = new_pheno,
                                 featureData = Biobase::featureData(x),
                                 annotation = Biobase::annotation(x))
  
  return(eset)
  
}