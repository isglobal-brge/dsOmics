#' Title
#'
#' @param genoData 
#' @param snpBlock 
#'
#' @return
#' @export
#'
#' @examples
fastAlleleFrequencyDS <- function(genoData, snpBlock){
  
  if(inherits(genoData, "GenotypeData")){
    geno <- GWASTools::GenotypeBlockIterator(genoData, snpBlock=snpBlock)
  } else if (inherits(genoData, "GdsGenotypeReader")) {
    gds_gd <- GWASTools::GenotypeData(genoData)
    geno <- GWASTools::GenotypeBlockIterator(gds_gd, snpBlock=snpBlock)
  } else {
    stop()
  }
  
  gfile <- openfn.gds(geno@data@filename, allow.duplicate = T, allow.fork = T)
  n <- index.gdsn(gfile, "genotype")
 
  iterations <- length(geno@snpFilter)
  
  # sums <- matrix(data = NA, ncol = ncol(geno@scanAnnot@data), ncol = iterations)
  
  sums_col <- Reduce(rbind, lapply(1:iterations, function(x){
    snp.sel <- geno@snpFilter[[x]]
    geno_iteration <- readex.gdsn(n, list(NULL, snp.sel), .value = 3, .substitute = NA)
    return(cbind(colSums(geno_iteration, na.rm = T),
          colSums(is.na(geno_iteration))))
  }))
  
  rs <- GWASTools::getVariable(geno, "snp.rs.id")[sums_col[,1] != 0]

  sums_col_no_0 <- sums_col[sums_col[,1] != 0,]
  sums_col_no_0[,2] <- length(GWASTools::getVariable(geno, "sample.id")) - sums_col_no_0[,2]
  
  sums_col_tot <- sums_col_no_0[,1] / (2 * sums_col_no_0[,2])
  
  sums_col_tot[sums_col_tot > 0.5] <- 1 - sums_col_tot[sums_col_tot > 0.5]
  
  return(tibble(rs = rs, n = sums_col_no_0[,2], MAF = sums_col_tot))
  
}