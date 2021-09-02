#' @title Get Ploygenic Risk Score
#' 
#' @details This function resolves a list of resources subsetting them by the 
#' SNPs of risk, this does not ensure that all the SNPs of risk will be found on the 
#' data. From all the found SNPs of risk, if an individual has less than 'snp_threshold' (percetage)
#' of SNPs with data, it will be dropped (SNP with no data is marked on the VCF as ./.). If an individual 
#' passes this threshold filter but still has SNPs with no data, those SNPs will be counted on the 
#' polygenic risk score as non-risk-alleles, to take this infomation into account, the number of SNPs 
#' with data for each individual is returned as 'n_snps'.
#'
#' @param resources \code{list} of all the VCF resources with biallelic genotype information. It is advised to 
#' have one VCF resource per chromosome, a big VCF file with all the information is always slower 
#' to use.
#' @param snp_threshold \code{numeric} (default \code{80}) Threshold to drop individuals. See details for 
#' further information.
#' @param ... Corresponds to the ROI table passed from the client. It is assembled if not NULL.
#'
#' @return \code{data.frame} were the rownames are the individuals. The columns found are: \cr
#' prs: Polygenic risk score per individual \cr
#' prs_nw: Polygenic risk score without weights (weight 1 for each risk allele) \cr
#' n_snps: Number of SNPs with information for each individual
#' @export
#'

PRSDS <- function(resources, snp_threshold, ...){
  
  ROI <- unlist(list(...))
  ROI <- data.frame(matrix(ROI[1:(length(ROI)-1)], ncol = as.numeric(ROI[length(ROI)])))
  if(ncol(ROI) == 5){
    colnames(ROI) <- c("chr_name", "start", "end", "effect_allele", "effect_weight")
    ROI$start <- as.numeric(ROI$start)
    ROI$end <- as.numeric(ROI$end)
    ROI$chr_name <- as.numeric(ROI$chr_name)
    ROI$effect_weight <- as.numeric(ROI$effect_weight)
  } else if(ncol(ROI) == 3){
    colnames(ROI) <- c("rsID", "effect_allele", "effect_weight")
    ROI$effect_weight <- as.numeric(ROI$effect_weight)
  } else{stop()}
  
  # Get the found SNPs, positions and alleles and merge into data frame
  found_rs <- do.call(c,lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.rs.id")
  }))
  # Check if no rs's have been found on the supplied vcf resources. Interrupt execution if so
  if(length(found_rs) == 0){
    stop('No SNPs from the PRS list found on the provided resources')
  }
  found_positions <- do.call(c,lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.position")
  }))
  found_alleles <- do.call(c,lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.allele")
  }))
  resource_name <- do.call(c,lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.rs.id")
  }))
  if(length(resources) == 1){
    individuals <- GWASTools::getVariable(resources[[1]], "sample.id")
  } else {
    individuals <- if(class(Reduce(function(x,y) if (identical(x,y)) x else FALSE, lapply(resources, function(x){
      GWASTools::getVariable(x, "sample.id")}))) != "logical"){
      GWASTools::getVariable(resources[[1]], "sample.id")
    } else {stop('Different individuals among the provided resources')}
  }
  
  # Get GDS parameters
  gds_info <- data.frame(rsID = found_rs, 
                         start = found_positions, 
                         reference_allele = do.call(c, lapply(found_alleles, function(x){substr(x, 1, 1)})), 
                         alternate_allele = do.call(c, lapply(found_alleles, function(x){substr(x, 3, 3)})),
                         resource_name = resource_name)
  if(colnames(ROI) == "rsID"){
    gds_with_risks <- dplyr::left_join(gds_info, ROI, by = "rsID")
  } else {
    gds_with_risks <- dplyr::left_join(gds_info, ROI, by = "start")
  }
  # Detect if effect_allele is on the reference_allele or alternate_allele (or none)
  #   'inverse' means that the reference_allele of the vcf is the effect_allele (risk allele)
  #   'correct' means that the alternate_allele of the vcf is the effect_allele
  #   'not_present' means that the effect_allele is not present on the alternate or reference allele
  allele_placement <- NULL
  for(i in seq(1, nrow(gds_with_risks))){
    if(gds_with_risks$effect_allele[i] %in% gds_with_risks$reference_allele[i]){
      allele_placement <- c(allele_placement, "inverse")
    } else if(gds_with_risks$effect_allele[i] %in% gds_with_risks$alternate_allele[i]) {
      allele_placement <- c(allele_placement, "correct")
    } else {
      allele_placement <- c(allele_placement, "not_present")
    }
  }
  gds_with_risks <- cbind(gds_with_risks, allele_placement)
  # Get genotype
  geno <- data.frame(do.call(cbind, lapply(resources, function(x){
    GWASTools::getVariable(x, "genotype")
  })))
  colnames(geno) <- found_rs
  # Get the percentage of SNPs available for every individual (where 100% is all the SNPs found
  # on the VCFs, not all the SNPs proposed by the ROI)
  # Also get the n_snps to return with the PRS results
  percentage_snps <- NULL
  n_snps <- NULL
  total_snps <- length(found_rs)
  for(i in seq(1, nrow(geno))){
    snps <- sum(geno[i,] != 3)
    percentage_snps <- c(percentage_snps, 100*(snps / total_snps))
    n_snps <- c(n_snps, snps)
  }
  names(n_snps) <- individuals
  # Remap geno; 0 means double alternate allele, 1 means single alternate allele, 2 means no alternate allele
  # remap to: 2 double alternate allele, 1 single alternate allele, 0 no alternate allele
  # the above is done taking into account if the reference and alternate allele are the effect allele
  # or not (calculated above)
  geno <- sapply(colnames(geno), function(x){
    if(gds_with_risks[which(gds_with_risks == x),]$allele_placement == "correct"){
      # Remove missing data (geno == 3) by setting as 2 which represent no alternate alleles
      geno[,x][geno[,x] == 3] <- 2
      plyr::mapvalues(geno[,x], from = c(0, 1, 2), to = c(2, 1, 0))
    } else if(gds_with_risks[which(gds_with_risks == x),]$allele_placement == "inverse") {
      # Remove missing data (geno == 3) by setting as 0 which represent no reference alleles
      geno[,x][geno[,x] == 3] <- 0
      geno[,x]
    } else if(gds_with_risks[which(gds_with_risks == x),]$allele_placement == "not_present") {
      NULL
    }
  })
  rownames(geno) <- individuals
  
  # Calculate PRS
  prs <- geno %*% gds_with_risks$effect_weight
  # Remove individuals that have less than snp_threshold percentage_snps
  prs[which(percentage_snps < snp_threshold)] <- NA
  # calculate PRS without weights
  prs_nw <- rowSums(geno)
  # TODO probability of prs_nw (snpassoc)
  p_prs_nw <- SNPassoc::pscore(prs_nw, colnames(geno))
  # TODO devolver todo como una tabla mejor??
  return(data.frame(prs = prs, prs_nw = prs_nw, p_prs_nw = p_prs_nw, n_snps = n_snps))
}
