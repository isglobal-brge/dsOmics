#' Title
#'
#' @param resources 
#' @param pgs_id 
#'
#' @return
#' @export
#'
#' @examples
PRSDS <- function(resources, pgs_id, snp_threshold){
  # Get PGS data to later calculate genetic risk score
  ROI <- .retrievePGS(pgs_id)
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
  individuals <- if(do.call(identical, lapply(resources, function(x){
    GWASTools::getVariable(x, "sample.id")}))){
    GWASTools::getVariable(resources[[1]], "sample.id")} else {
      stop('Different individuals among the provided resources')
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
  # If weight_type is present and is equal to "OR" or "HR", convert the effect_weight
  # to log(effect_weight) to get the beta
  # This is done row by row in case not all rows have the same weight_type
  if(!is.null(gds_with_risks$weight_type)){
    for(i in seq(1, nrow(gds_with_risks))){
      if(c("OR", "HR") %in% gds_with_risks$weight_type[i]){
        gds_with_risks$effect_weight[i] <- log(gds_with_risks$effect_weight[i])
      }
    }
  }
  # Calculate PRS
  prs <- geno %*% gds_with_risks$effect_weight
  # Remove individuals that have less than snp_threshold percentage_snps
  prs[which(percentage_snps < snp_threshold)] <- NA
  # calculate PRS without weights
  prs_nw <- rowSums(geno)
  # TODO probability of prs_nw (snpassoc)
  # p_prs_nw <- SNPassoc::pscore(prs_nw, colnames(geno))
  # TODO devolver todo como una tabla mejor??
  return(list(prs = prs, prs_nw = prs_nw, n_snps = n_snps))
}

#' Title
#'
#' @param pgs_id 
#'
#' @return
#'
#' @examples
.retrievePGS <- function(pgs_id){
  pgs <- httr::GET(paste0("https://www.pgscatalog.org/rest/score/", pgs_id))
  pgs_text <- httr::content(pgs, "text")
  pgs_scoring_file <- jsonlite::fromJSON(pgs_text, flatten = TRUE)$ftp_scoring_file
  if(is.null(pgs_scoring_file)){stop('[', pgs_id, '] Not found in www.pgscatalog.org')}
  
  # Download scoring file
  destination <- tempfile(fileext = ".txt.gz")
  download.file(pgs_scoring_file, destination)
  
  # Unzip file
  unziped_file <- tempfile(fileext = ".txt")
  R.utils::gunzip(destination, unziped_file)
  
  # Read file into R
  # Logic https://www.pgscatalog.org/downloads/
  scorings <- read.delim(unziped_file, comment.char = "#")
  if(c("chr_name", "chr_position") %in% colnames(scorings)){
    data <- data.frame(chr_name = scorings$chr_name,
                       start = scorings$chr_position,
                       end = scorings$chr_position,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    if("weight_type" %in% colnames(scorings)){
      data <- data %>% tibble::add_column(weight_type = scorings$weight_type)}
    return(data)
  } else {
    data <- data.frame(rsID = scorings$rsID,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    if("weight_type" %in% colnames(scorings)){
      data <- data %>% tibble::add_column(weight_type = scorings$weight_type)}
    return(data)
  }
}
