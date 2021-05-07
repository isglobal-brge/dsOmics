#' @title Hardy-Weinberg Equilibrium testing
#' 
#' @description This function performs exact Hardy-Weinberg Equilibrium testing (using Fisher's Test) 
#' over a selection of SNPs. It also counts genotype, calculates allele frequencies, 
#' and calculates inbreeding coefficients.
#'
#' @param genoData \code{GenotypeData} object
#' @param chromosome \code{character} Chromosome to study. \code{"all"} to study all available chromosomes
#' @param geno.counts \code{bool} if \code{TRUE}, genotype counts are returned in the output data.frame
#' @param block.size \code{numeric}  number of SNPs to read in at once
#' @param permute \code{bool} logical indicator for whether to permute alleles before calculations
#' @param controls \code{bool} logical to calculate the HWE test only on the controls
#' @param controls_column \code{character} name of the case/controls column of the covariates 
#' used to create the GenotypeData object. Only used if \code{controls = TRUE}
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
exactHWEDS <- function(genoData, chromosome, geno.counts, block.size, permute, controls, 
                       controls_column){
  
  if(inherits(genoData, "GenotypeData")){
    chr <- getChromosome(genoData)
    
    if(chromosome != "all"){
      range_chr <- range(which(chr == chromosome))
      if(!(chromosome %in% autosomeCode(genoData))){
        warning(paste0("HWE is only available for autosomes, performing HWE with all available autosomes."))
        range_chr <- range(which(is.element(chr, 1:22)))
      }
      else if(any(is.infinite(range_chr))){
        warning(paste0("Couldn't find chromosome [", chromosome, "], performing HWE with all available autosomes."))
        range_chr <- range(which(is.element(chr, 1:22)))
      }
      snpStart <- range_chr[1]
      snpEnd <- range_chr[2]
    }
    else{
      range_chr <- range(which(is.element(chr, 1:22)))
      snpStart <- range_chr[1]
      snpEnd <- range_chr[2]
    }

    if(controls){
      if(missing(controls_column)){stop('"controls" is set to TRUE, but no "controls_column" was provided')}
      scan.exclude <- genoData@scanAnnot@data$scanID[genoData@scanAnnot@data[[controls_column]] == 1]
    } else{scan.exclude <- NULL}
    GWASTools::getScanID(genoData)
    ans <- GWASTools::exactHWE(genoData = genoData, 
                    geno.counts = geno.counts,
                    snpStart = snpStart,
                    snpEnd = snpEnd,
                    scan.exclude = scan.exclude,
                    block.size = block.size, 
                    verbose = FALSE,
                    permute = permute)
    rs <- data.frame(snpID = GWASTools::getVariable(genoData, "snp.id"),
                     rs = GWASTools::getVariable(genoData, "snp.rs.id"))
    ans <- dplyr::left_join(ans, rs, by = "snpID")
    return(tibble::as_tibble(ans) %>% mutate_at(c(3:6, 8:9), as.numeric) %>% 
             select(!c("snpID", "f")) %>% dplyr::relocate(rs))
  }
  else{
    stop(paste0("Object of incorrect type [", class(genoData), "] exactHWE requires object of type GenotypeData"))
  }
}