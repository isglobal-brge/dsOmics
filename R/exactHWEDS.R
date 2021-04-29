#' @title Hardy-Weinberg Equilibrium testing
#' 
#' @description This function performs exact Hardy-Weinberg Equilibrium testing (using Fisher's Test) 
#' over a selection of SNPs. It also counts genotype, calculates allele frequencies, 
#' and calculates inbreeding coefficients.
#'
#' @param genoData \code{GenotypeData} object
#' @param sexcol \code{character} Name of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object
#' @param male \code{character} Name of the male indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{male} than \code{Male})
#' @param female \code{character} Name of the female indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{female} than \code{Female})
#' @param chromosome \code{character} Chromosome to study. \code{"all"} to study all available chromosomes
#' @param geno.counts \code{bool} if \code{TRUE}, genotype counts are returned in the output data.frame
#' @param block.size \code{numeric}  number of SNPs to read in at once
#' @param permute \code{bool} logical indicator for whether to permute alleles before calculations
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
exactHWEDS <- function(genoData, sexcol, male, female, chromosome, geno.counts, block.size, permute){
  
  if(inherits(genoData, "GenotypeData")){
    if(sexcol %in% colnames(genoData@scanAnnot@data) == FALSE){
      stop(paste0("Selected sexcol [", sexcol, "] can't be found on the GenotypeData"))
    }
    if(!all(c(male, female) %in% unique(genoData@scanAnnot@data[,sexcol]))){
      stop(paste0("Incorrect male or female identifier [", male, "/", female,
                  "]. Available gender identifiers: ", paste0(unique(genoData@scanAnnot@data[,sexcol]), collapse = "/")))
    }

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
    genoData@scanAnnot@sexCol <- sexcol
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == male] <- "M"
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == female] <- "F"
    ans <- GWASTools::exactHWE(genoData = genoData, 
                    geno.counts = geno.counts,
                    snpStart = snpStart,
                    snpEnd = snpEnd,
                    block.size = block.size, 
                    verbose = FALSE,
                    permute = permute)
    rs <- GWASTools::getVariable(genoData, "snp.rs.id")
    return(tibble::as_tibble(ans) %>% mutate_at(c(3:6, 8:9), as.numeric) %>% 
             tibble::add_column(rs=rs) %>% select(!c("snpID")) %>% 
             dplyr::relocate(rs))
  }
  else{
    stop(paste0("Object of incorrect type [", class(genoData), "] exactHWE requires object of type GenotypeData"))
  }
}