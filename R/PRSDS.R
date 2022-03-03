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
#' @param ... Corresponds to the prs_table table passed from the client. It is assembled if not NULL.
#'
#' @return \code{data.frame} were the rownames are the individuals. The columns found are: \cr
#' prs: Polygenic risk score per individual \cr
#' prs_nw: Polygenic risk score without weights (weight 1 for each risk allele) \cr
#' n_snps: Number of SNPs with information for each individual
#' @export
#'

PRSDS <- function(resources, snp_threshold, snp_assoc, pgs_id, ...){

  # Get table of SNPs and betas
  if(!is.null(pgs_id)){ # Retrieve pgs_id table
    prs_table <- .retrievePGS(pgs_id)
  } else { # Construct custom table
    prs_table <- unlist(list(...))
    prs_table <- data.frame(matrix(prs_table[1:(length(prs_table)-1)], ncol = as.numeric(prs_table[length(prs_table)])))
    if(ncol(prs_table) == 5){
      colnames(prs_table) <- c("chr_name", "start", "end", "effect_allele", "effect_weight")
      prs_table$start <- as.numeric(prs_table$start)
      prs_table$end <- as.numeric(prs_table$end)
      prs_table$chr_name <- as.numeric(prs_table$chr_name)
      prs_table$effect_weight <- as.numeric(prs_table$effect_weight)
    } else if(ncol(prs_table) == 3){
      colnames(prs_table) <- c("rsID", "effect_allele", "effect_weight")
      prs_table$effect_weight <- as.numeric(prs_table$effect_weight)
    }
  }
  
  # Get the  SNPs, positions and alleles and merge into data frame
  positions <- lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.position")
  })
  chromosomes <- lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.chromosome")
  })
  rs <- lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.rs.id")
  })
  alleles <- lapply(resources, function(x){
    GWASTools::getVariable(x, "snp.allele")
  })

  # Get the matching indexes
  if("rsID" %in% colnames(prs_table)){ # Get rs
    match_indexes <- lapply(rs, function(x){
      matches <- x %in% prs_table$rsID
      if(any(matches)){
        which(matches)
      } else {
        NULL
      }
    })
  } else { # Get positions and chr number
    match_indexes <- lapply(positions, function(x){
      matches <- x %in% prs_table$start
      if(any(matches)){
        which(matches)
      } else {
        NULL
      }
      
    })
  }
  
  # Stop execution if no match_indexes are found
  if(is.null(unlist(match_indexes))){
    stop('No matching SNPs found')
  }

  # Check that individuals are consistent
  if(length(resources) == 1){
    individuals <- GWASTools::getVariable(resources[[1]], "sample.id")
  } else {
    individuals <- if(class(Reduce(function(x,y) if (identical(x,y)) x else FALSE, lapply(resources, function(x){
      GWASTools::getVariable(x, "sample.id")}))) != "logical"){
      GWASTools::getVariable(resources[[1]], "sample.id")
    } else {stop('Different individuals among the provided resources')}
  }
  
  # Get GDS parameters
  gds_info <- do.call(rbind, lapply(1:length(match_indexes), function(x){
    if(is.null(match_indexes[[x]])){
      NULL
    } else {
      data.frame(rsID = rs[[x]][match_indexes[[x]]],
                 start = positions[[x]][match_indexes[[x]]],
                 reference_allele = unlist(lapply(alleles[[x]][match_indexes[[x]]], function(y){substr(y, 1, 1)})), 
                 alternate_allele = unlist(lapply(alleles[[x]][match_indexes[[x]]], function(y){substr(y, 3, 3)})),
                 chr_name = as.numeric(chromosomes[[x]][match_indexes[[x]]]))
    }
  }))
  
  # Merge gds_indo and prs_table for matching individuals
  if("rsID" %in% colnames(prs_table)){
    gds_with_risks <- dplyr::left_join(gds_info, prs_table, by = "rsID")
  } else {
    gds_with_risks <- dplyr::left_join(gds_info, prs_table, by = c("start", "chr_name"))
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
  geno <- do.call(cbind, lapply(1:length(resources), function(x){
    if(is.null(match_indexes[[x]])){
      NULL
    } else {
      snps <- do.call(cbind, lapply(1:length(match_indexes[[x]]), function(y){
        GWASTools::getVariable(resources[[x]], 
                               "genotype", 
                               start=c(1,match_indexes[[x]][y]), count=c(-1,1))
      }))
      return(snps)
    }
    
  }))
  colnames(geno) <- gds_info$rsID
  # Get the percentage of SNPs available for every individual (where 100% is all the SNPs found
  # on the VCFs, not all the SNPs proposed by the prs_table)
  # Also get the n_snps to return with the PRS results
  percentage_snps <- NULL
  n_snps <- NULL
  total_snps <- length(gds_info$rsID)
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
    if(gds_with_risks[which(gds_with_risks$rsID == x),]$allele_placement == "correct"){
      # Remove missing data (geno == 3) by setting as 2 which represent no alternate alleles
      geno[,x][geno[,x] == 3] <- 2
      plyr::mapvalues(geno[,x], from = c(0, 1, 2), to = c(2, 1, 0))
    } else if(gds_with_risks[which(gds_with_risks$rsID == x),]$allele_placement == "inverse") {
      # Remove missing data (geno == 3) by setting as 0 which represent no reference alleles
      geno[,x][geno[,x] == 3] <- 0
      geno[,x]
    } else if(gds_with_risks[which(gds_with_risks$rsID == x),]$allele_placement == "not_present") {
      # NULL
      rep(0, nrow(geno))
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
  if(snp_assoc){
    p_prs_nw <- tryCatch({
      SNPassoc::pscore(prs_nw, colnames(geno))
    }, error = function(w){
      "SNPassoc (>2.0-3) not available on the Opal"
    })
  } else {
    p_prs_nw <- rep(NA, length(prs))
  }

  # TODO devolver todo como una tabla mejor??
  return(data.frame(prs = prs, prs_nw = prs_nw, p_prs_nw = p_prs_nw, n_snps = n_snps))
}

#' @title Auxiliary function to merge the PRS results to a table by ID
#' @export
PRSDS_aux <- function(prs_results, prs_results_name, table, id){
  if(!(id %in% colnames(table))){
    stop('[',id,'] is not present on the table to be merged.')
  }
  if(is.na(prs_results$p_prs_nw)){
    prs_results <- prs_results %>% 
      select(prs, prs_nw, n_snps) %>% 
      rename(!!paste0('prs_', prs_results_name) := prs, 
             !!paste0('prs_nw_', prs_results_name) := prs_nw,
             !!paste0('n_snps_', prs_results_name) := n_snps)
  } else {
    prs_results <- prs_results %>% 
      select(prs, prs_nw, p_prs_nw, n_snps) %>% 
      rename(!!paste0('prs_', prs_results_name) := prs, 
             !!paste0('prs_nw_', prs_results_name) := prs_nw,
             !!paste0('p_prs_nw_', prs_results_name) := p_prs_nw,
             !!paste0('n_snps_', prs_results_name) := n_snps)
  }
  
  return(merge(table, prs_results,
        by.x = colnames(table)[which(names(table) == id)], by.y = "row.names")[,union(names(table), names(prs_results))])
}

#' @title Internal function: Get PGS catalog table of polygenic risks
#'
#' @param pgs_id \code{character} ID of the PGS catalog to be used to calculate the polygenic risk score. 
#' Polygenic Score ID & Name from https://www.pgscatalog.org/browse/scores/
#'
#' @return \code{data.frame} with the columns: \cr
#' - If chr_name and poition are found: \cr
#' + start \cr
#' + end \cr
#' + effect_allele \cr
#' + effect_weight \cr
#' + weight_type (if present on the catalog) \cr
#' - If rsID is found: \cr
#' + rsID \cr
#' + effect_allele \cr
#' + effect_weight \cr
#' + weight_type (if present on the catalog) \cr
#'

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
  } else {
    data <- data.frame(rsID = scorings$rsID,
                       # reference_allele = scorings$reference_allele,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    if("weight_type" %in% colnames(scorings)){
      data <- data %>% tibble::add_column(weight_type = scorings$weight_type)}
  }
  # If weight_type is present and is equal to "OR" or "HR", convert the effect_weight
  # to log(effect_weight) to get the beta
  # This is done row by row in case not all rows have the same weight_type
  if(!is.null(data$weight_type)){
    for(i in seq(1, nrow(data))){
      if(c("OR", "HR") %in% data$weight_type[i]){
        data$effect_weight[i] <- log(data$effect_weight[i])
      }
    }
    data$weight_type <- NULL
  }
  
  return(data)
  
}

.recodeprs_table <- function(scorings){
  if(c("chr_name", "chr_position") %in% colnames(scorings)){
    data <- data.frame(chr_name = scorings$chr_name,
                       start = scorings$chr_position,
                       end = scorings$chr_position,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    return(data)
  } else {
    data <- data.frame(rsID = scorings$rsID,
                       effect_allele = scorings$effect_allele,
                       effect_weight = scorings$effect_weight)
    return(data)
  }
}
