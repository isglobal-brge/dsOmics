#' @title Genome-wide association analysis (GWAS) on the server-side
#' 
#' @description Performs GWAS using \code{GENESIS} Bioconductor package
#' @param genoData a \code{GenotypeData} object which is a container for storing genotype data
#' from a GWAS toghether with the metadata associated with the subjects (i.e. phenotypes and/or covariates)
#' and SNPs
#' @param outcome A character string specifying the name of the outcome variable in \code{genoData}
#' @param covars A vector of character strings specifying the names of the fixed effect covariates 
#' in \code{genoData}; an intercept term is automatically included. If NULL (default) the only fixed effect 
#' covariate is the intercept term
#' @param family A description of the generalized linear model used. The defatul is "binomial" that is defined
#' for case/control studies. Quantitative traits can be analyzed by using "gaussian". Other values are accepted.
#' @param snpBlock an integer specifying the number of SNPs in an iteration block. See \code{GenotypeIterator} 
#' function in GWASTools package.
#'  
#' @param ... other arguments of \code{fitNullModel} function in GENESIS package
#' @return a matrix with SNPs ordered by p-values
#' 
#' @author Gonzalez, JR.
#'
#' @import dplyr
#'
#' @export 
#' 
GWASDS <- function(genoData, outcome, covars=NULL, family="binomial", snpBlock, ...){
  if(!is.null(covars)){
    covars <- unlist(strsplit(covars, split=","))
  }
  nullmod <- GENESIS::fitNullModel(genoData, outcome = outcome, 
                                   covars = covars, 
                                   family = family, ...)
  genoIterator <- GWASTools::GenotypeBlockIterator(genoData, snpBlock=snpBlock)
  assoc <- GENESIS::assocTestSingle(genoIterator, null.model = nullmod)
  assoc$rs <- GWASTools::getVariable(genoData, "snp.rs.id")[assoc$variant.id]
  alleles <- GWASTools::getVariable(genoData, "snp.allele")[assoc$variant.id] # ref/alt
  assoc$ref_allele <- substring(alleles, 1, 1)
  assoc$alt_allele <- substring(alleles, 3, 3)
  ans <- assoc %>% as_tibble() %>%
    select(variant.id, rs, everything()) %>% 
    select(!c("Score", "Score.SE", "Score.Stat", "PVE", "MAC")) %>%
    dplyr::rename(p.value=Score.pval) %>% select(!c("variant.id")) %>%
    arrange(p.value) %>% dplyr::relocate(p.value, .after = pos)

  #############################################################
  # CAPTURE THE diffP SETTINGS
  nfilter.diffP.epsilon <- getOption("default.nfilter.diffP.epsilon")
  nfilter.diffP.resampleN <- getOption("default.nfilter.diffP.resampleN")
  #############################################################
  
  if(!is.null(nfilter.diffP.epsilon) & !is.null(nfilter.diffP.resampleN)){
    # Resample `nfilter.diffP.resampleN` times the GWAS data (removing a random ID each time)
    l1.sens <- lapply(1:nfilter.diffP.resampleN, function(x){
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
      # Fit the same model to the subset
      nullmod2 <- GENESIS::fitNullModel(new_gds, outcome = outcome, 
                                        covars = covars, 
                                        family = family, ...)
      genoIterator2 <- GWASTools::GenotypeBlockIterator(new_gds, snpBlock=snpBlock)
      assoc2 <- GENESIS::assocTestSingle(genoIterator2, null.model = nullmod2)
      assoc2$rs <- GWASTools::getVariable(new_gds, "snp.rs.id")[assoc2$variant.id]
      alleles2 <- GWASTools::getVariable(new_gds, "snp.allele")[assoc2$variant.id] # ref/alt
      assoc2$ref_allele <- substring(alleles2, 1, 1)
      assoc2$alt_allele <- substring(alleles2, 3, 3)
      ans2 <- assoc2 %>% as_tibble() %>%
        select(variant.id, rs, everything()) %>% 
        select(!c("Score", "Score.SE", "Score.Stat", "PVE", "MAC", "chr", "pos", 
                  "Score.pval", "n.obs", "ref_allele", "alt_allele", "variant.id", "freq"))
      results <- dplyr::left_join(ans %>% select(c("rs", "Est", "Est.SE")), ans2, by = "rs")
      results <- results %>% dplyr::mutate(est_sens = abs(Est.x - Est.y)) %>%
        mutate(estSE_sens = abs(Est.SE.x - Est.SE.y)) %>% select(c("rs", "est_sens", "estSE_sens")) %>%
        replace(is.na(.), 0)
      
      return(results)
    })
    
    # Extract max l1-sensitivities
    l1.sens <- Reduce(function(...){dplyr::left_join(..., by = "rs")}, l1.sens)

    cmd <- parse(text = paste0("c('",paste(colnames(l1.sens)
                                           [grepl("^est_sens", colnames(l1.sens))],
                                           collapse = "','"),"')"))
    est_sens <- matrixStats::rowMaxs(as.matrix(l1.sens[,eval(cmd)]))
    est_sens[which(est_sens == 0)] <- .Machine$double.xmin
    cmd <- parse(text = paste0("c('",paste(colnames(l1.sens)
                                           [grepl("^estSE_sens", colnames(l1.sens))],
                                           collapse = "','"),"')"))
    estSE_sens <- matrixStats::rowMaxs(as.matrix(l1.sens[,eval(cmd)]))
    estSE_sens[which(estSE_sens == 0)] <- .Machine$double.xmin
    
    
    laplace_noise <- Laplace_noise_generator(m = 0,
                                             b = est_sens/nfilter.diffP.epsilon,
                                             n.noise = nrow(ans))
    quants <- quantile(laplace_noise,c(0.05,0.95))
    laplace_noise[laplace_noise < quants[1]] <- quants[1]
    laplace_noise[laplace_noise > quants[2]] <- quants[2]
    ans_diffP <- ans %>% mutate(Est := Est + laplace_noise)
    
    
    laplace_noise2 <- Laplace_noise_generator(m = 0,
                                              b = estSE_sens/nfilter.diffP.epsilon,
                                              n.noise = nrow(ans))
    quants <- quantile(laplace_noise2,c(0.05,0.95))
    laplace_noise2[laplace_noise2 < quants[1]] <- quants[1]
    laplace_noise2[laplace_noise2 > quants[2]] <- quants[2]
    ans_diffP <- ans_diffP %>% mutate(Est.SE := Est.SE + laplace_noise2)
    
    ans_diffP <- ans_diffP %>% mutate(p.value := 2 * pnorm(-abs(Est / Est.SE)))
    
    # -abs(Est / Est.SE) = qnorm(p.value/2)
    # laplace_noise = abs(laplace_noise) / qnorm(ans$p.value/2)
    # if(!any(c(est_sens, estSE_sens) == 0)){ # TODO revisar aquest any no fa falta crec
    #   
    # }
    # TODO do not return freq!!!
    # else {ans_diffP <- ans}
    return(ans_diffP)
  } else {
    return(ans)
  }
}



