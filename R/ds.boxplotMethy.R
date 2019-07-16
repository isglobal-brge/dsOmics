##' Generating boxplot summary of methylation values split by phenotypic variable levels 
##' 
##' Function generates a boxplot summary of methylation data by phenotypic condition 
##' for a specified CpG site. If multiple studies are analysed the boxplots provide representations
##' of the pooled data. Within a single chart a boxplot is generated for each level of the user 
##' specified phenotypic variable.   
##' @title Boxplot summary of methylation values split by phenotypic variable levels
##' @param cpg name of the CpG site the boxplot will be created for
##' @param phenoVar phenotypic variable, each level of which a boxplot will be created for
##' @param molecular.data name of the DataSHIELD object to which the methylation data set has been assigned
##' @param pheno.data name of the DataSHIELD object to which the phenotype data set has been assigned
##' @param datasources Opal object or list of opal objects denoting the opal server(s) information
##' @export
##' @examples
##' 

ds.boxplotMethy <- function(cpg, phenoVar, molecular.data, pheno.data, datasources){
  
  #Calling 'histDataByCondition' server side function which returns a list containing a 
  #vector of methylation values for each level of the phenotypic variable in 'methyl.data'
  cally <- paste0("returnDataByCondition(cpg=\'", cpg, "\', ", "phenoVar=\'", phenoVar, "\', " ,"molecular.data=", molecular.data, ", ", "pheno.data=", pheno.data,")")
  datashield.assign(datasources, 'methyl.data', as.symbol(cally))
  
  #Retrieving names of phenotypic variable levels for each study
  studyPhenoLevels = ds.names(x = 'methyl.data')
  
  #Function which tests for errors returned from server side function in each study
  errorCheck <- function(x){
    return(which(c("phenoLevelsError", "cpgError", "phenoVarError") %in%  x))
  }
  
  #Checking if a server side error was returned and printing corresponding message to console
  if(any(unlist(lapply(studyPhenoLevels, errorCheck)))){
    
    errorMessages=list(paste0("The number of levels associated with the phenotypic variable exceeds what can be plot.\n", 
                           "Max number of phenotypic levels permitted = 10", sep=''),
                    "The cpg site specified does not exist in one or more of the studies", 
                    "The phenotypic variable specified does not exist in one or more of the studies")
    
    message(paste(" ", "ERROR(S): ", 
                  paste0(errorMessages[unique((unlist(lapply(studyPhenoLevels, errorCheck))))], 
                         collapse = "\n"),
                  " ", sep="\n"))
    
  }else{
    
    #Checking that the phenotypic levels returned for each study match
    if(length(studyPhenoLevels)>1){
      for(i in 1:(length(studyPhenoLevels)-1)){
        if(!identical(studyPhenoLevels[[i]], studyPhenoLevels[[i+1]])){
          return(message(paste0(" ", 
                                paste0("ERROR: ", "Phenotypic levels for variable '", phenoVar, 
                                       "' do not match across all studies", sep = ""), 
                                " ", sep="\n")))
        }
      }
    }
    
    phenoLevels = studyPhenoLevels[[1]]
    
    #Creating a list to store summary data for each phenotypic level
    plotList = vector("list", length(phenoLevels))
    
    counter=1
    
    #Looping over levels of phenotypic variables and storing summary data
    for (i in phenoLevels) {
      
      ds.assign(toAssign = paste('methyl.data$', i, sep = ''), newobj = 'phenoLeveliData')
      ds.asNumeric(x='phenoLeveliData', newobj = 'numPhenoLeveliData')
      
      #NB:ds.quantileMean returns an error if length(x) < 5
      if(ds.length(x='numPhenoLeveliData') >= 5){
        
        quantData = ds.quantileMean(x='numPhenoLeveliData', type = 'combine')
        
        #Extracting relevant quantile data for boxplot
        plotList[[counter]] = list(stats = matrix(c(quantData["5%"], quantData["25%"], quantData["50%"], 
                                                    quantData["75%"], quantData["95%"]) ,5,1), n=10)
      }
      
      counter = counter + 1
    }
    
    
    #Grouping summary data from each level (only for non null entries in  plotList)
    groupbxp <- do.call(mapply, c(cbind, Filter(Negate(is.null), plotList)))
    nonNullLevels = phenoLevels[which(sapply(plotList, Negate(is.null)))]
    
    #Generating boxplot
    bxp(groupbxp, main = paste0("Boxplot of Methylation Values for ", cpg, sep=""), show.names = FALSE)
    axis(1, at=1:length(nonNullLevels), labels=nonNullLevels )
    
    #Generate warning which outputs all phenotypic levels not plotted (<5 observations)
    if(length(nonNullLevels) != length(phenoLevels)){
      message(paste(" ", "NOTE: The following levels of the phenotypic variable were excluded",
                    "      from the boxplot as they had fewer than 5 observations:", " ",
                    paste0("      ", paste0(phenoLevels[which(sapply(plotList, is.null))], collapse = ", "), sep=""),
                    " ", sep="\n"))
    }
    
    if(length(phenoLevels) == 1){
      message(paste(" ", "NOTE: Only 1 level detected for phenotypic variable specified and ",
                    "      therefore boxplot generated for all samples in data set", " ", sep="\n"))
      
    }
    
  }
  
}
