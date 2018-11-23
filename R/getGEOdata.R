# Illustrative data are obtained from public repository (GEO)
# This script is used to harmonized the data
# and to select some features to avoid computational issues

library(Biobase)

gg <- GEOquery::getGEO("GSE66351")
gse66351 <- gg[[1]]
names(pData(gse66351))[c(41, 45, 49)] <- c("age", "casecon", "Sex")
pData(gse66351)$age <- as.numeric(pData(gse66351)$age)

cpgs <- featureNames(gse66351)
sel <- sample(cpgs, 500)
gse66351.sel <- gse66351[sel,]


gg <- GEOquery::getGEO("GSE80970")
gse80970 <- gg[[1]]
names(pData(gse80970))[c(39, 41, 43)] <- c("age", "casecon", "Sex")
pData(gse80970)$age <- as.numeric(pData(gse80970)$age)
pData(gse80970)$casecon[pData(gse80970)$casecon == "control"] <- "CTRL"
pData(gse80970)$casecon[pData(gse80970)$casecon == "Alzheimer's disease"] <- "AD"

temp <- sample(cpgs, 20)
sel2 <- unique(c(sel, temp))
gse80970.sel <- gse80970[sel2,]


save(gse66351.sel, file="data/GSE19711.Rdata")
save(gse80970.sel, file="data/GSE80970.Rdata")

