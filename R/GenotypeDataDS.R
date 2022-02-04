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
#' @param columnId \code{character} Column of the covars that contains the IDs
#' @param sexId \code{character} (default \code{NULL}) Column of the covars that contains the sex phenotype
#' @param male_encoding \code{character} (default \code{"male"}) String used to encode the male sex
#' phenotype on the covars table
#' @param female_encoding \code{character} (default \code{"female"}) String used to encode the female sex
#' phenotype on the covars table
#' @param case_control_column \code{character} (default \code{NULL}) Column that holds the
#' case/control to relevel to 0/1
#' @param case \code{character} (default \code{NULL}) Encoding of the case of the \code{case_control_column}
#' @param control \code{character} (default \code{NULL}) Encoding of the control of the \code{case_control_column}
#'
#' @return ...
#' @author Gonzalez, JR.
#'
#' @export

GenotypeDataDS <- function(x, covars, columnId, sexId, male_encoding, female_encoding,
                           case_control_column, case, control){
  # Decode hex variables (character at origin)
  male_encoding <- if(!is.null(male_encoding)){rawToChar((wkb::hex2raw(male_encoding)))}
  female_encoding <- if(!is.null(female_encoding)){rawToChar((wkb::hex2raw(female_encoding)))}
  case <- if(!is.null(case)){rawToChar((wkb::hex2raw(case)))}
  control <- if(!is.null(control)){rawToChar((wkb::hex2raw(control)))}
  columnId <- which(names(covars) == rawToChar((wkb::hex2raw(columnId))))
  sexId <- if(!is.null(sexId)){which(names(covars) == rawToChar((wkb::hex2raw(sexId))))}
  case_control_column <- if(!is.null(case_control_column)){
    val <- rawToChar(wkb::hex2raw(case_control_column))
    if(!(val %in% names(covars))){
      stop('The [case_control_column] supplied (', val, 
           ') is not present on the covars table')
      } else {val}
    } 
  
  names(covars)[columnId] <- "scanID"
  if(!is.null(sexId)){
    # TODO does the sexId supplied exist?
    if(colnames(covars)[sexId] != "sex"){
      covars <- covars %>% tibble::add_column(sex = unlist(covars[, sexId]))
    }
    # TODO does the male/female_encoding supplied exist?
    covars$sex[covars$sex %in% male_encoding] <- "M"
    covars$sex[covars$sex %in% female_encoding] <- "F"
  }
  # TODO check if the supplied `case` and `control` are even present on the dataset
  if(!is.null(case_control_column)){
    covars[[case_control_column]][covars[[case_control_column]] %in% case] <- 1
    covars[[case_control_column]][covars[[case_control_column]] %in% control] <- 0
    covars[[case_control_column]] <- as.numeric(covars[[case_control_column]])
  }

  covars_id <- covars$scanID
  geno_id <- getScanID(x)

  covars <- covars[covars_id %in% geno_id,]
  covars_id <- covars$scanID
  covars <- covars[match(geno_id, covars_id),]
  
  # TODO commented code corresponds to possible no.match argument to substitute NAs created by having more
  # geno individuals than pheno individuals. To be revised on the future.
  # if(!is.null(no.match)){
  #   # Check if no.match is numeric
  #   if(!is.na(as.numeric(no.match))){ # If can be converted successfully use as numeric
  #     no.match <- as.numeric(no.match)
  #   } # Otherwise use as it is
  #   rows_to_put_no_match <- is.na(covars$scanID)
  #   
  #   covars <- as_tibble(lapply(covars, function(x){
  #     tryCatch({
  #       x[rows_to_put_no_match] <- no.match
  #       return(x)
  #     }, error = function(w){return(x)})
  #   }))
  # }
  
  if(length(geno_id) > length(covars_id)){
    covars$scanID <- geno_id
  }
  covars <- covars[!is.na(covars$scanID),]
  
  scanAnnot <- GWASTools::ScanAnnotationDataFrame(data.frame(covars))
  geno <- GWASTools::GenotypeData(x, scanAnnot = scanAnnot)
  geno
}