library(resourcer)
library(DSOpal)
library(dsBaseClient)
library(DSLite)
library(SNPRelate)
source("tests/getClientVCF.R")


################   GWAS data analysis

dslite.server <- newDSLiteServer(resources = list(
  vcf = newResource(name = "vcf", url = "https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/1000G.vcf", format = "VCF2GDS"),
  covars = newResource(name = "covars", url = "https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/covars_1000G.txt", format = "ssv")
))

resourcer::registerResourceResolver(GDSFileResourceResolver$new())


builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "dslite.server", resource = "vcf", driver = "DSLiteDriver")
logindata <- builder$build()
logindata
# login and assign resources
conns <- datashield.login(logins = logindata, assign = TRUE, symbol = "vcf.res")

datashield.assign.expr(conns, symbol = "gds", expr = quote(as.resource.object(vcf.res)))
# assign covars resource
datashield.assign.resource(conns, symbol = "covars.res", resource = list(study1 = "covars"))
datashield.assign.expr(conns, symbol = "covars", expr = quote(as.resource.data.frame(covars.res)))

ds.ls()
ds.class("gds")
ds.class("covars")# inspect data in the first study (DSLite only!)
dslite.server$getSessionData(conns[[1]]@sid, "gds")
dslite.server$getSessionData(conns[[1]]@sid, "covars")


dslite.server$aggregateMethod("indexGdsnDS", "gdsfmt::index.gdsn")
dslite.server$aggregateMethod("readGdsnDS", "gdsfmt::read.gdsn")
dslite.server$aggregateMethod("snpgdsGetGenoDS", "SNPrelate::snpgdsGetGeno")


cally <- "readGdsnDS(indexGdsnDS(gds, 'sample.id'))"
ids <- datashield.aggregate(conns,  cally)
if (!identical(ids,vars[,1]))
 stop("There VCF ids are not identical to those in the first column of 'covars' ")


cally <- "snpgdsGetGenoDS(gds, snp.id=1)"
datashield.aggregate(conns, as.symbol(cally))

dslite.server$aggregateMethod("selSNP", function(x, snp, vars) {
  g <- as.numeric(snpgdsGetGeno(gds.f, snp.id = snps.id[1]))
  
  data.frame(x[feature, ])[, c(feature, vars)]
})




# .... R
gds.f<-SNPRelate::snpgdsOpen("c:/tmp/gds.f.gds", readonly = FALSE)
covars <- read.delim("c:/juan/CREAL/GitHub/brgeUtils/data/covars_1000G.txt", sep="")

sample.id <- read.gdsn(index.gdsn(gds.f, "sample.id"))
snps.id <- read.gdsn(index.gdsn(gds.f, "snp.id"))
snps <- read.gdsn(index.gdsn(gds.f, "snp.rs.id"))
g <- as.numeric(snpgdsGetGeno(gds.f, snp.id = snps.id[1]))
g



datashield.logout(conns)



##############  First try with Yannick's function

source("tests/getClientVCF.R")


# register resolver so that magic happens

resourcer::registerResourceResolver(GDSFileResourceResolver$new())



# use it (method and snpfirstdim are optional)

res <- resourcer::newResource(url = "https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/1000G.vcf", format = "VCF2GDS")

client <- resourcer::newResourceClient(res)


gds.f <- client$getValue()

sample.id <- read.gdsn(index.gdsn(gds.f, "sample.id"))
snps.id <- read.gdsn(index.gdsn(gds.f, "snp.id"))
snps <- read.gdsn(index.gdsn(gds.f, "snp.rs.id"))
g <- as.numeric(snpgdsGetGeno(gds.f, snp.id = snps.id[1]))
g

client$close()

