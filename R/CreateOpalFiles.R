#Script to create phenotype and gene expression data csv files with associated xls dictionaries

library(Biobase)
library(xlsx)

#Input directory for ExpressenSet data file and output directory for csv & dictionary files
#By default both set to working directory
inputDIR <- getwd()
otuputDIR <- getwd()

#Getting ExpressionSet data file names
Data_Files <- list.files(path = inputDIR, pattern = "*.Rdata")

for (ES_File in Data_Files) {
  
  #Load data from ExpressionSet file
  tmp_env <- new.env()
  tmp_obj <- load(paste(inputDIR, ES_File, sep = "/"), tmp_env)
  ES_Data <- get(tmp_obj, tmp_env)
  rm(tmp_env)
  
  #Data set label used in table and file names
  DataLabel <-  sub(".Rdata", "", ES_File)
  
  #----------Phenotype Dictionary and Data file----------#
  
  pVars <- sub(":", "_", varLabels(ES_Data))
  pVarNum <- length(pVars) #Number of phenotype variables  
  ptableName <- rep(paste("pheno", DataLabel, sep = ""), pVarNum) #Set phenotype table name
  pEntType <- rep("Participant", pVarNum) #Set value for entityType
  pRepeatable <- rep(0, pVarNum) #Set value for repeatable
  
  #Determining the data type of each phenotype variable
  pVarTypes <- 
    unlist(lapply(lapply(pData(ES_Data), class), function(x) ifelse(x == "numeric", "integer","text")))
  
  #Phenotype dictionary file data frame
  phenoDict <- data.frame(table=ptableName, 
                          name=pVars, 
                          valueType=pVarTypes, 
                          entityType=pEntType,
                          referencedEntityType=character(pVarNum),
                          mimeType=character(pVarNum),
                          unit=character(pVarNum),
                          repeatable=pRepeatable,
                          occuranceGroup=character(pVarNum),
                          `label:en`= pVars)
  
  #Creating data frame with pheno data 
  ES_PhenoData <- cbind(sample=sampleNames(ES_Data), pData(ES_Data))
  rownames(ES_PhenoData) <-  1:length(sampleNames(ES_Data))
  colnames(ES_PhenoData) <- sub(":", "_", colnames(ES_PhenoData))
  
  #Creating dictionary and data files
  write.xlsx2(phenoDict, paste("pheno", DataLabel, "-dictionary.xls", sep = ""), 
              sheetName = "Variables", col.names = TRUE, row.names = FALSE, append = FALSE)
  
  write.csv(ES_PhenoData, paste("pheno", DataLabel, ".csv", sep = ""))
  
  #------------------------------------------------------#
  
  
  #----------Gene Expression Dictionary and Data file----------#
  
  SampleNum <- length(sampleNames(ES_Data)) #Number of samples
  gExpNumRows <- SampleNum + 1 #Number of rows in gene expression dictionary file
  
  #Set values for gene expression dictionary variables
  gtableName <- rep(paste("geneExp", DataLabel, sep = ""), gExpNumRows) 
  gValueType <- rep("decimal", SampleNum)
  gEntType <- rep("Participant", gExpNumRows)
  glabels <- rep("sample", SampleNum)
  gRepeatable <- rep(0, gExpNumRows)
  
  #Gene expression dictionary file data frame
  geneExpDict <- data.frame(table=gtableName,
                            name=c("feature", sampleNames(ES_Data)), 
                            valueType=c("text", gValueType), 
                            entityType=gEntType,
                            referencedEntityType=character(gExpNumRows),
                            mimeType=character(gExpNumRows),
                            unit=character(gExpNumRows),
                            repeatable=gRepeatable,
                            occuranceGroup=character(gExpNumRows),
                            `label:en`= c("feature", glabels))
  
  #Creating data frame with gene expression data
  ES_GeneExpData <- cbind(feature=featureNames(ES_Data), exprs(ES_Data))
  rownames(ES_GeneExpData) <-  1:length(featureNames(ES_Data))
  
  #Creating dictionary and data files
  write.xlsx2(geneExpDict, paste("geneExp", DataLabel, "-dictionary.xls", sep = ""), 
              sheetName = "Variables", col.names = TRUE, row.names = FALSE, append = FALSE)
  
  write.csv(ES_GeneExpData, paste("geneExp", DataLabel, ".csv", sep = ""))
  
  #------------------------------------------------------------#
}