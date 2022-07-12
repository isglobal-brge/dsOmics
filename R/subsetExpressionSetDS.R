#' @title Subset ExpressionSet
#' 
#' @description Subset ExpressionSet using a categorical variable of the covariates as filter. Can also subset by 
#' complete cases
#'
#' @param eSet \code{ExpressionSet} ExpressionSet to subset
#' @param objective_variable \code{character} (default \code{NULL}) Name of the covariate on the ExpressionSet to use as filter
#' @param objective_value \code{character} (default \code{NULL}) Name of the value from the \code{objective_variable} to filter. 
#' The resulting subset will be the individuals that match this value. 
#' To put in in code, it can be represented as: \code{subset <- expressionSet[expressionSet$objective_variable == objective_value,]}
#' @param complete_cases \code{bool} (default \code{FALSE}) If \code{TRUE} only the complete cases will be included on the subset. 
#' This option can be used with \code{objective_variable} and \code{objective_value} or without them, if those arguments are not 
#' present, the subset will be the complete cases of the whole ExpressionSet.
#'
#' @return Subseted \code{ExpressionSet}
#' @export

subsetExpressionSetDS <- function(eSet, objective_variable, objective_value, complete_cases){
  subset_indexes <- subset_indexes_complete_cases <- T
  
  # Extract the phenotype
  pheno <- eSet@phenoData@data
  
  # Check if objective_variable exists on the phenotypes
  if(!(objective_variable %in% colnames(pheno))){
    stop('The selected objective_variable [', objective_variable, 
         "] is not present on the ExpressionSet. Check which variables are available using the function ds.varLabels()")
  }
  
  if(!is.null(objective_variable) & !is.null(objective_value)){
    # Coerce variable into a character for consistency when using a numeric or character variable to subset
    objective_variable_asCharacter <- as.character(pheno[[objective_variable]])
    
    # Get the resulting subset size
    subset_indexes <- objective_variable_asCharacter == objective_value
    # The missings of the objective_variable turn out as NAs, we want them to be FALSE
    subset_indexes[is.na(subset_indexes)] <- FALSE
  } else if (!is.null(objective_variable) & is.null(objective_value) | is.null(objective_variable) & !is.null(objective_value)) {
    stop("If a subset by a phenotype is desired, make sure to introduce both arguments `objective_variable` and `objective_value`")
  }

  # Update the subset_indexes if complete_cases is TRUE
  if(complete_cases){
    subset_indexes_complete_cases <- complete.cases(pheno)
  }
  
  # Resulting subset
  subset_indexes <- as.logical(subset_indexes * subset_indexes_complete_cases)
  
  # Get subset size
  subset_size <- sum(subset_indexes, na.rm = T)
  
  # Disclosure control adapted from dsBase::dataFrameSubsetDS1 to ensure the resulting subset is not potentially disclosive
  #########################################################################
  # DataSHIELD MODULE: CAPTURE THE nfilter SETTINGS
  nfilter.subset <- as.numeric(getOption("nfilter.subset"))
  #########################################################################
  
  # CHECK SUBSET LENGTH IS CONSISTENT WITH nfilter FOR MINIMUM SUBSET SIZE
  
  if(subset_size < nfilter.subset){
    studysideMessage <- "Subset to be created is too small (<nfilter.subset)"
    stop(studysideMessage, call. = FALSE)
  }
  
  neweSet <- eSet[,subset_indexes]
  
  return(neweSet)

}
