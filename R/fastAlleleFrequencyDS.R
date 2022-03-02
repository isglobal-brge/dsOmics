#' @title Allelic frequency
#' 
#' @description Calculates the frequency of the A allele on a server side GenotypeData object.
#'
#' @param genoData \code{GenotypeData} or \code{GdsGenotypeReader} Genotype data
#' @param snpBlock \code{integer} number of SNPs to read at each iteration, 
#' tune this parameter to improve performance.
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
    stop('genoData object is not of class ["GenotypeData" or "GdsGenotypeReader"]')
  }
  
  #############################################################
  # MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  #nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  
  scan_number <- GWASTools::nscan(genoData)
  
  if(scan_number < nfilter.tab){
    stop("Disclosure issue: Not enough individuals to aggregate: N.individuals less than nfilter.tab")
  }
  
  gfile <- openfn.gds(geno@data@filename, allow.duplicate = T, allow.fork = T)
  n <- index.gdsn(gfile, "genotype")
 
  iterations <- length(geno@snpFilter)
  
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