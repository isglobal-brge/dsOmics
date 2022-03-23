#' @title Allelic frequency
#' 
#' @description Calculates the frequency of the A allele
#'
#' @param genoData \code{GenotypeData} \code{\link{GenotypeData}} object
#'
#' @return A matrix with a row for each SNP. Columns "M" for males, "F" for females, and "all" for all scans give 
#' frequencies of the A allele. Sample size for males, females, and all is returned as "n.M", "n.F", and "n", respectively.
#'  "MAF" is the minor allele frequency over all scans.
#' @export
#' 
#' @import dplyr

alleleFrequencyDS <- function(genoData){
  if(inherits(genoData, "GenotypeData")){
    
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
    
    ans <- GWASTools::alleleFrequency(genoData, verbose = FALSE)
    rs <- GWASTools::getVariable(genoData, "snp.rs.id")
    ans <- tibble::as_tibble(ans) %>% tibble::add_column(rs=rs) %>% dplyr::relocate(rs)
    
    # diffP
    #############################################################
    # CAPTURE THE diffP SETTINGS
    nfilter.diffP.epsilon <- getOption("default.nfilter.diffP.epsilon")
    nfilter.diffP.resampleN <- getOption("default.nfilter.diffP.resampleN")
    #############################################################
    
    if(!is.null(nfilter.diffP.epsilon) & !is.null(nfilter.diffP.resampleN)){
      # Resample `nfilter.diffP.resampleN` times the alleleFreq (removing a random ID each time)
      resamp <- do.call(rbind, lapply(1:nfilter.diffP.resampleN, function(x){
        # Select resample individuals
        individuals <- GWASTools::getVariable(genoData, "sample.id")
        individuals_resample <- individuals[-sample(1:length(individuals), 1)]
        # New temp file
        new_f <- tempfile()
        # Subset all SNPs all individuals - random one
        gdsSubset2(genoData@data@filename, new_f,
                   sample.include=individuals_resample, snp.include=NULL,
                   sub.storage=NULL,
                   compress="LZMA_RA",
                   verbose=TRUE,
                   allow.fork=TRUE)
        new_gds <- GWASTools::GdsGenotypeReader(new_f, allow.fork = TRUE)
        # Add pheno information to subset
        new_gds <- GWASTools::GenotypeData(new_gds, 
                                           scanAnnot = genoData@scanAnnot[genoData@scanAnnot@data$scanID %in% 
                                                                            individuals_resample])
        ans2 <- GWASTools::alleleFrequency(new_gds, verbose = FALSE)
        ans2 <- tibble::as_tibble(ans2) %>% tibble::add_column(rs=rs) %>% dplyr::relocate(rs)
        # Merge resample with original results
        merged_data <- merge(ans, ans2, by = "rs")
        # Get l1-sensitivity of freq
        freq_sens <- max(merged_data$MAF.x - merged_data$MAF.y)
        # Return l1-sensitivities
        return(data.frame(freq_sens = freq_sens))
      }))
      
      # Extract max l1-sensitivities
      freq_sens <- max(resamp$freq_sens)
      # Add Laplacian noise to  freq with mean mean 0 and
      # scale l1-sensitivity / nfilter.diffP.epsilon
      ans_diffP <- ans %>% mutate(MAF := MAF + Laplace_noise_generator(m = 0, 
                                                                       b = freq_sens/nfilter.diffP.epsilon, 
                                                                       n.noise = nrow(ans)))
      # MAF filter
      #############################################################
      # CAPTURE THE diffP SETTINGS
      default.nfilter.MAF <- getOption("default.nfilter.MAF")
      #############################################################
      if(!is.null(default.nfilter.MAF)){
        ans_diffP <- ans_diffP %>% filter(MAF > default.nfilter.MAF)
      }
      return(ans_diffP)
    } else {
      # MAF filter
      #############################################################
      # CAPTURE THE diffP SETTINGS
      default.nfilter.MAF <- getOption("default.nfilter.MAF")
      #############################################################
      if(!is.null(default.nfilter.MAF)){
        ans <- ans %>% filter(MAF > default.nfilter.MAF)
      }
      return(ans)
    }
    
  }
  else(stop(paste0("Object of incorrect type [", class(genoData), "] alleleFrequency requires object of type GenotypeData")))
}