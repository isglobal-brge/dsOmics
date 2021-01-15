#' @title Add Phenotype data to ExpressionSet
#'
#' @param x \code{ExpressionSet} ExpressionSet to which add phenotype information
#' @param pheno \code{data.frame} Table with the new phenotypes
#' @param identifier \code{character} Name of the ID column on the phenotypes data.frame
#'
#' @return
#' @export
#'
#' @examples
addPhenoDataDS <- function(x, pheno, identifier){
  
  if(!(any(identifier %in% colnames(pheno)))){
    stop("Identifier [", identifier, "] is not on the phenotypes table")
  }
  
  og_pheno <- Biobase::pData(x)
  og_pheno_md <- Biobase::varMetadata(x)
  
  new_variables <- colnames(pheno)[!(identifier == colnames(pheno))]
  old_variables <- colnames(og_pheno)
  
  og_individuals <- rownames(og_pheno)
  new_individuals <- pheno[,identifier]
  common_individuals <- new_individuals %in% og_individuals
  
  new_pheno <- pheno[common_individuals,]
  og_pheno <- cbind(og_pheno, og_individuals_id = og_individuals)
  
  new_pheno <- dplyr::left_join(og_pheno, new_pheno, by = c("og_individuals_id" = identifier))
  rownames(new_pheno) <- new_pheno$og_individuals_id
  new_pheno$og_individuals_id <- NULL
  
  if(any(new_variables %in% old_variables)){stop("Variables conflict between ExpressionSet and new PhenoData")}
  for(i in new_variables){
    og_pheno_md <- eval(str2expression(paste0("rbind(og_pheno_md, ", i, " = NA)")))
  }
  
  new_pheno <- new("AnnotatedDataFrame", data=new_pheno, varMetadata=og_pheno_md)
  
  eset <- Biobase::ExpressionSet(assayData = Biobase::exprs(x),
                                 phenoData = new_pheno,
                                 featureData = Biobase::featureData(x),
                                 annotation = Biobase::annotation(x))
  
  return(eset)
  
}