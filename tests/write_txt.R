library(Biobase)
setwd("c:/juan/CREAL/GitHub/dsOmics/data/")
load("GSE66351.Rdata")
load("GSE80970.Rdata")

ss <- exprs(gse66351.sel)
geneExpr <- data.frame(ProbeId=rownames(ss), ss)
pheno <- pData(gse66351.sel)
pheno$id <- rownames(pheno)

write.csv(geneExpr, file="gse66351-gene.csv", quote=FALSE, 
            row.names = FALSE)
write.csv(pheno, file="gse66351-pheno.csv", quote=FALSE, 
            row.names = FALSE)

ss <- exprs(gse80970.sel)
geneExpr <- data.frame(ProbeId=rownames(ss), ss)
pheno <- pData(gse80970.sel)
pheno$id <- rownames(pheno)

write.csv(geneExpr, file="gse80970-gene.csv", quote=FALSE, 
            row.names = FALSE)
write.csv(pheno, file="gse80970-pheno.csv", quote=FALSE, 
            row.names = FALSE)
