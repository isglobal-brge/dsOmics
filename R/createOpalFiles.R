##' Creating genetic expression and phenotype dictionary and data files 
##' 
##' Function processes an Expression Set object or csv/txt files with
##' genetic expression and related phenotype data to create xls dictionary
##' and csv data files used for Opal data table creation in DataShield. 
##' @title Creating genetic expression and phenotype dictionary and data files
##' @param geneExpData an Expression Set object or name of genetic expression csv/txt file
##' @param phenoFile parameter used with csv/txt files to specify name of phenotype file
##' @param inputDIR optional parameter used with csv/txt files to specify directory containing data files
##' @param estimateCellCounts optional TRUE/FALSE parameter which determines whether epigenomic data is calculated and added to phenotypic data file
##' @param cellTypeRef cell type reference used to estimate cell counts
##' @param geneExpTableName optional parameter which sets the gene expression table name in Opal
##' @param phenoTableName optional parameter which sets the phenotype data table name in Opal
##' @export
##' @examples
##' 
##' 
##' 

createOpalFiles <- function(geneExpData, phenoFile, inputDIR = getwd(), estimateCellCounts = FALSE, 
                            cellTypeRef = "blood gse35069 complete", geneExpTableName = "geneExp",
                            phenoTableName = "pheno"){

  require(Biobase)
  require(xlsx)
  require(meffil)
  
  if(class(geneExpData) == "ExpressionSet"){
    
    #-----Variables for creating Gene Expression Files-----#
    SampleNum <- length(sampleNames(geneExpData)) 
    gSamples <- sampleNames(geneExpData) 
    
    gDataSet <- cbind(feature=featureNames(geneExpData), exprs(geneExpData))
    rownames(gDataSet) <-  1:length(featureNames(geneExpData))
    
    #-----Variables for creating Phenotype Files-----#
    pDataSet <- cbind(sample=gSamples, pData(geneExpData))
    rownames(pDataSet) <-  1:length(sampleNames(geneExpData))
    
    #Removing special characters from phenotype variable names 
    pVars <- c("sample", gsub(" |/|:", "_", varLabels(geneExpData)))
    
    #Variable for calculating epigenomic data
    cpgs <- exprs(geneExpData)
    
    
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
    pVars <- gsub(" |/|:", "_", colnames(pDataSet)) #List of phenotype variable names
    
    #Epigenomic data variables
    cpgs <- gDataSet[,-1]
    rownames(cpgs) <- gDataSet[,1]
    cpgs <- as.matrix(cpgs)
  }
  
  #Calculating epigenomic variables
  if(isTRUE(estimateCellCounts)){
    
    cellCounts <- tryCatch(
      meffil.estimate.cell.counts.from.betas(cpgs, cell.type.reference = cellTypeRef, verbose = TRUE), 
      
      #Error handling function will return NULL to cellCounts object if an error is encountered
      error=function(e){
        message(paste(" ", "WARNING: Dataset provided cannot be used to estimate cell counts", 
                      "Phenotype csv data file will be created without epigenomic variables", 
                      "Set input argument estimateCellCounts = FALSE to avoid this warning", " ", sep="\n"))
      })
    
    if(!is.null(cellCounts) ){
      #Standardizing all epigenomic variables
      cellCounts_scaled = scale(cellCounts)
      colnames(cellCounts_scaled) = paste(colnames(cellCounts_scaled), "Scaled", sep="_")
      
      #Add epigenomic variables to phenotype data.frame
      pDataSet <- cbind(pDataSet, cellCounts, cellCounts_scaled)
      pVars <- append(pVars, c(colnames(cellCounts), colnames(cellCounts_scaled)))
    }  
    
  }
  
  #----------Gene Expression Dictionary and Data file----------#
  
  gExpNumRows <- SampleNum + 1 #Number of rows in gene expression dictionary file
  
  #Set values for gene expression dictionary variables
  gtableName <- rep(geneExpTableName, gExpNumRows) 
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
                            occurrenceGroup=character(gExpNumRows))
  
  `label:en` <- c("feature", glabels)
  geneExpDict <- cbind(geneExpDict, `label:en`)
  
  #Data frame with values for categories tab in dictionary file
  geneExpCategories <- data.frame(table=NA,
                                  variable=NA,
                                  name=NA,
                                  code=NA,
                                  missing=NA)
  
  `label:en` <- NA
  geneExpCategories <- cbind(geneExpCategories, `label:en`)
  
  #Creating dictionary and data files
  write.xlsx2(geneExpDict, paste0(geneExpTableName, "-dictionary.xls", sep=""), sheetName = "Variables", 
              col.names = TRUE, row.names = FALSE, append = FALSE)
  write.xlsx2(geneExpCategories, paste0(geneExpTableName, "-dictionary.xls", sep=""), sheetName = "Categories", 
              col.names = TRUE, row.names = FALSE, append = TRUE)
  
  write.csv(gDataSet, paste0(geneExpTableName, ".csv", sep=""), na = "")
  
  #------------------------------------------------------------#
  
  #----------Phenotype Dictionary and Data file----------#
  
  pVarNum <- length(pVars) #Number of phenotype variables  
  ptableName <- rep(phenoTableName, pVarNum) #Set phenotype table name
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
  
  #Removing special characters from phenotype variable names 
  colnames(pDataSet) <- gsub(" |/|:", "_", colnames(pDataSet))
  
  #Phenotype dictionary file data frame
  phenoDict <- data.frame(table=ptableName, 
                          name=pVars, 
                          valueType=pVarTypes, 
                          entityType=pEntType,
                          referencedEntityType=character(pVarNum),
                          mimeType=character(pVarNum),
                          unit=character(pVarNum),
                          repeatable=pRepeatable,
                          occurrenceGroup=character(pVarNum))
  
  `label:en`<- pVars
  phenoDict <- cbind(phenoDict, `label:en`)
  
  #Data frame with values for "categories" tab in dictionary file
  phenoCategories <- data.frame(table=NA,
                                variable=NA,
                                name=NA,
                                code=NA,
                                missing=NA)
  
  `label:en`<- NA
  phenoCategories <- cbind(phenoCategories, `label:en`)
  
  #Creating dictionary and data files
  write.xlsx2(phenoDict, paste0(phenoTableName, "-dictionary.xls"), sheetName = "Variables", 
              col.names = TRUE, row.names = FALSE, append = FALSE)
  write.xlsx2(phenoCategories, paste0(phenoTableName, "-dictionary.xls"), sheetName = "Categories", 
              col.names = TRUE, row.names = FALSE, append = TRUE)
  
  write.csv(pDataSet, paste0(phenoTableName, ".csv"), na = "")
  
  #------------------------------------------------------#
  
}
