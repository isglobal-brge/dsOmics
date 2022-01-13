#' Title
#'
#' @param genoData 
#' @param snp_list 
#'
#' @return
#' @export
#'
#' @examples
subsetGenoDS <- function(genoData, snp_list){
  
  # TODO disclosure check for subset of snps dimensions
  genoData_og <- get(genoData, envir = parent.frame())
  
  # TODO check that genoData_og inherits something that makes sense
    
  if(snp_list == "ethnic_snps"){
    data("ethnic_snps")
    snp_subset <- ethnic_snps
    
  } else {
    # TODO cas en que el usuari pase la seva llista de snps
    # snp_subset <- ???
  }
  
  
    
  original_rs <- GWASTools::getVariable(genoData_og, "snp.rs.id")
  original_ids <- GWASTools::getVariable(genoData_og, "snp.id")
  
  ids_of_interest <- original_ids[which(original_rs %in% snp_subset)]
  
  GWASTools::close(genoData_og)
  
  new_f <- tempfile()
  
  # TODO, depenen si es un gds sol o un gdsData amb el fenotip, sha de canviar la seguent liniea
  # ia que el gds sol no te el slot @data

  gdsSubset(genoData_og@data@filename, new_f,
            sample.include=NULL, snp.include=ids_of_interest,
            sub.storage=NULL,
            compress="LZMA_RA",
            block.size=5000,
            verbose=TRUE,
            allow.fork=TRUE)

  # TODO, depenen si es un gds sol o un gdsData amb el fenotip, sha de canviar la seguent liniea
  # ia que el gds sol no te el slot @data
  
  genoData_og@data <- GWASTools::GdsGenotypeReader(genoData_og@data@filename, allow.fork = TRUE)
  assign(genoData, genoData_og, envir = parent.frame())
  
  new_gds <- GWASTools::GdsGenotypeReader(new_f, allow.fork = TRUE)
  
  return(new_gds)
  
}