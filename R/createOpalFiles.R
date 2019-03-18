##' Creating genetic expression and phenotype dictionary and data files 
##' 
##' Function processes an Expression Set object or csv/txt files with
##' genetic expression and related phenotype data to create xls dictionary
##' and csv data files used for Opal data table creation in DataShield. 
##' @title Creating genetic expression and phenotype dictionary and data files
##' @param geneExpData an Expression Set object or name of genetic expression csv/txt file
##' @param phenoFile parameter used with csv/txt files to specify name of phenotype file
##' @param inputDIR optional parameter used with csv/txt files to specify directory containing data files
##' @export
##' @examples
##' 
##' 
##' 

createOpalFiles <- function(geneExpData, phenoFile, inputDIR = getwd()){

  require(Biobase)
  require(xlsx)
  
  if(class(geneExpData) == "ExpressionSet"){
    
    #-----Variables for creating Gene Expression Files-----#
    SampleNum <- length(sampleNames(geneExpData)) 
    gSamples <- sampleNames(geneExpData) 
    
    gDataSet <- cbind(feature=featureNames(geneExpData), exprs(geneExpData))
    rownames(gDataSet) <-  1:length(featureNames(geneExpData))
    
    #-----Variables for creating Phenotype Files-----#
    pDataSet <- cbind(sample=gSamples, pData(geneExpData))
    rownames(pDataSet) <-  1:length(sampleNames(geneExpData))
    pVars <- c("sample", sub(":", "_", varLabels(geneExpData)))
    
    
  } else if(class(geneExpData) == "character"){
    
    if(grepl("\\.csv", geneExpData) == TRUE){
      
      #Reading in gene expression data set
      gDataSet <- read.csv(paste(inputDIR, geneExpData, sep = "/"))
      
      #Reading in phenotype data set
      pDataSet <- read.csv(paste(inputDIR, phenoFile, sep = "/"))
      
      
    } else if(grepl("\\.txt", geneExpData) == TRUE){
      
      #Reading in gene expression data set
      gDataSet <- read.delim(paste(inputDIR, geneExpData, sep = "/"))
      
      #Reading in phenotype data set
      pDataSet <- read.delim(paste(inputDIR, phenoFile, sep = "/"))
    } else {
      
      stop("Invalid File Name Provided")
    }
    
    #-----Variables for creating Gene Expression Files-----#
    
    gSamples <- colnames(gDataSet)[-1] #List of sample names
    SampleNum <- length(gSamples) #Number of samples
    colnames(gDataSet)[1] <- "feature"
    
    #-----Variables for creating Phenotype Files-----#
    
    pDataSet <- cbind(sample=gSamples, pDataSet) #Data frame with phenotype data
    pVars <- sub(":", "_", colnames(pDataSet)) #List of phenotype variable names
    
  }
  
  #----------Gene Expression Dictionary and Data file----------#
  
  gExpNumRows <- SampleNum + 1 #Number of rows in gene expression dictionary file
  
  #Set values for gene expression dictionary variables
  gtableName <- rep("geneExp", gExpNumRows) 
  gValueType <- rep("decimal", SampleNum)
  gEntType <- rep("Participant", gExpNumRows)
  glabels <- rep("sample", SampleNum)
  gRepeatable <- rep(0, gExpNumRows)
  
  #Gene expression dictionary file data frame
  geneExpDict <- data.frame(table=gtableName,
                            name=c("feature", gSamples), 
                            valueType=c("text", gValueType), 
                            entityType=gEntType,
                            referencedEntityType=character(gExpNumRows),
                            mimeType=character(gExpNumRows),
                            unit=character(gExpNumRows),
                            repeatable=gRepeatable,
                            occuranceGroup=character(gExpNumRows))
  
  `label:en` <- c("feature", glabels)
  geneExpDict <- cbind(geneExpDict, `label:en`)
  
  #Data frame with values for categories tab in dictionary file
  geneExpCategories <- data.frame(table=character(),
                                  variable=character(),
                                  name=character(),
                                  code=character(),
                                  missing=integer())
  
  `label:en` <- character()
  geneExpCategories <- cbind(geneExpCategories, `label:en`)
  
  #Creating dictionary and data files
  write.xlsx2(geneExpDict, "geneExp-dictionary.xls", sheetName = "Variables", 
              col.names = TRUE, row.names = FALSE, append = FALSE)
  write.xlsx2(geneExpCategories, "geneExp-dictionary.xls", sheetName = "Categories", 
              col.names = TRUE, row.names = FALSE, append = TRUE)
  
  write.csv(gDataSet, "geneExp.csv")
  
  #------------------------------------------------------------#
  
  #----------Phenotype Dictionary and Data file----------#
  
  pVarNum <- length(pVars) #Number of phenotype variables  
  ptableName <- rep("pheno", pVarNum) #Set phenotype table name
  pEntType <- rep("Participant", pVarNum) #Set value for entityType
  pRepeatable <- rep(0, pVarNum) #Set value for repeatable
  
  #Function used to set the valueType for each phenotype variable
  setPhenoVarType <- function(x){
    if(x == 'numeric'){
      vType = 'decimal'
    } else if(x == 'integer'){
      vType = 'integer'
    } else{
      vType = 'text'
    }
    return(vType)
  }
  
  #Determining the data type of each phenotype variable
  pVarTypes <- unlist(lapply(lapply(pDataSet, class), setPhenoVarType))
  
  #Phenotype dictionary file data frame
  phenoDict <- data.frame(table=ptableName, 
                          name=pVars, 
                          valueType=pVarTypes, 
                          entityType=pEntType,
                          referencedEntityType=character(pVarNum),
                          mimeType=character(pVarNum),
                          unit=character(pVarNum),
                          repeatable=pRepeatable,
                          occuranceGroup=character(pVarNum))
  
  `label:en`<- pVars
  phenoDict <- cbind(phenoDict, `label:en`)
  
  #Data frame with values for "categories" tab in dictionary file
  phenoCategories <- data.frame(table=character(),
                                variable=character(),
                                name=character(),
                                code=character(),
                                missing=integer())
  
  `label:en`<-character()
  phenoCategories <- cbind(phenoCategories, `label:en`)
  
  colnames(pDataSet) <- sub(":", "_", colnames(pDataSet))
  
  
  #Creating dictionary and data files
  write.xlsx2(phenoDict, "pheno-dictionary.xls", sheetName = "Variables", 
              col.names = TRUE, row.names = FALSE, append = FALSE)
  write.xlsx2(phenoCategories, "pheno-dictionary.xls", sheetName = "Categories", 
              col.names = TRUE, row.names = FALSE, append = TRUE)
  
  write.csv(pDataSet, "pheno.csv")
  
  #------------------------------------------------------#
  
}
