returnDataByCondition <- function(cpg, phenoVar, molecular.data, pheno.data){
  
  #NB: Assign type function in DataShield
  
  errorList = list(phenoLevelsError = FALSE, cpgError = FALSE, phenoVarError = FALSE)
  
  #Checking if the phenotypic variable specified exists
  if((phenoVar %in% colnames(pheno.data)) == FALSE){
    errorList$phenoVarError = TRUE
  }else{
    #Determining the levels of the phenotype variable
    phenoLevels = levels(as.factor(pheno.data[, phenoVar]))
    phenoLevelCount = length(phenoLevels)
    
    #Checking if the number of phenotypic levels 
    #exceeds 10 in which case it can't be plotted
    if(phenoLevelCount > 10){
      errorList$phenoLevelsError = TRUE
    }
    
  }
  
  #Checking if the cpg site specified exists
  if((cpg %in% molecular.data$feature) == FALSE){
    errorList$cpgError = TRUE
  }
  
  #If any errors found in previous checks return error list
  #otherwise process data normally 
  if(any(sapply(errorList, isTRUE))){
    
    errorList = Filter(isTRUE, errorList)
    return(errorList)
    
  }else{
    
    #Creating a list to store the methylation data for each level of the phenotype variable
    methyl.data = vector("list", phenoLevelCount)
    names(methyl.data) = gsub(" |/|:|'", "_", phenoLevels)
    
    #For levels which are integers only add "Level_" prefix to level name in order to 
    #avoid errors on the client side when reading from data frame and referencing levels
    names(methyl.data) = lapply(names(methyl.data), function(x) ifelse(!grepl("\\D", x), paste0("Level_", x), x))
    
    for (k in 1:phenoLevelCount) {
      
      #List of phenotypic data row numbers corresponding to 
      #level[k] of the phenotypic variable of interest
      sampleRowNums = which(pheno.data[[phenoVar]]==phenoLevels[k])
      
      #Getting the list of sample names which correspond to each level[k] of the phenotypic variable
      phenoSampleNames = pheno.data[sampleRowNums, 'sample']
      
      #Retrieving methlyation values of each sample corresponding to specified 
      #level of phenotypic variable (only for denoted CpG site)
      methyl.data[[k]] = t(molecular.data[which(molecular.data$feature==cpg), 
                                          which(names(molecular.data) %in% phenoSampleNames)])
      
    }
    
    return(methyl.data)
    
  }
  
  
  
}

