# load required libraries
library(Biobase)

# load funtion doing linear model
source("R/lmFeature.R")

# get data
load("data/GSE19711.Rdata")
gse66351.sel

#
# GET THE REQUIRED DATASETS (molecular and phenotypes)
# 


# 1. molecular data -> features in rows, samples in columns

molecular <- exprs(gse66351.sel)


# 2. phenotypic data -> samples in rows, variables in columns
pheno <- pData(gse66351.sel)

# condition variable can be accessed as
head(pheno$casecon)
table(pheno$casecon)

# covariates can be accessed as
head(pheno$Sex)

# assess linear regression for a given CpG (cg05343755)
ii <- which(rownames(molecular)=="cg05343755")
lmFeature(ii, "casecon ~ Sex",
          molecular.data=molecular, 
          pheno.data=pheno)



# assess linear regression to the first two CpGs
ans <- lapply(1:10, lmFeature, "casecon ~ Sex", 
              molecular.data=molecular, 
              pheno.data=pheno)
out <- do.call(rbind, ans)
out


