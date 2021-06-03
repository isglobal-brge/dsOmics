#' @title Class GenotypeData. 
#' 
#' @description Container for storing genotype data from a GWAS together with the 
#' metadata associated with the subjects (e.g. phenotypes and covariates) and SNPs 
#' involved in the study. This is a wrapper of \code{GenotypeData} function from GWASTools package.
#'
#' @param x a \code{GdsGenotypeReader} object (see GWASTools). It is the object that is
#' obtained in the opal servers
#' from a resource of type VCF2GDS
#' @param covars a data.frame or a tibble having the metadata of samples (i.e. phenotypes and/or covariates)
#' @param columnId \code{numeric} Column of the covars that contains the IDs
#' @param sexId \code{numeric} (default \code{NULL}) Column of the covars that contains the sex phenotype
#' @param male_encoding \code{character} (default \code{"male"}) String used to encode the male sex
#' phenotype on the covars table
#' @param female_encoding \code{character} (default \code{"female"}) String used to encode the female sex
#' phenotype on the covars table
#' @param case_control_column \code{character} (default \code{NULL}) Name of the column that holds the
#' case/control to relevel to 0/1
#' @param case \code{character} (default \code{NULL}) Encoding of the case of the \code{case_control_column}
#' @param control \code{character} (default \code{NULL}) Encoding of the control of the \code{case_control_column}
#' @param na_string \code{character vector} (default \code{NULL}) Encoding to be
#' put to NA of the \code{case_control_column}
#'
#' @return ...
#' @author Gonzalez, JR.
#'
#' @export

GenotypeDataDS <- function(x, covars, columnId, sexId, male_encoding, female_encoding,
                           case_control_column, case, control, ...){
  na_string <- unlist(list(...))
  names(covars)[columnId] <- "scanID"
  if(!is.null(sexId)){
    if(colnames(covars)[sexId] != "sex"){
      covars <- covars %>% tibble::add_column(sex = unlist(covars[, sexId]))
    }
    covars$sex[covars$sex %in% male_encoding] <- "M"
    covars$sex[covars$sex %in% female_encoding] <- "F"
  }
  
  if(!is.null(case_control_column)){
    covars[[case_control_column]][covars[[case_control_column]] %in% case] <- 1
    covars[[case_control_column]][covars[[case_control_column]] %in% control] <- 0
    if(!is.null(na_string)){
      covars[[case_control_column]][covars[[case_control_column]] %in% na_string] <- NA
    }
    covars[[case_control_column]] <- as.numeric(covars[[case_control_column]])
  }
  covars_id <- covars$scanID
  geno_id <- getScanID(x)
  if(length(geno_id) > length(covars_id)){
    stop('The covariates table is missing the individuals: ', 
         paste(geno_id[!(geno_id %in% covars_id)], collapse = ", "))
  }
  covars <- covars[covars_id %in% geno_id,]
  covars <- covars[match(geno_id, covars_id),]
  
  scanAnnot <- GWASTools::ScanAnnotationDataFrame(data.frame(covars))
  geno <- GWASTools::GenotypeData(x, scanAnnot = scanAnnot)
  geno
}