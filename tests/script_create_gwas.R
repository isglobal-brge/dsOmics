#### Create new
library(resourcer)
library(DSLite)
library(dsBaseClient)
library(dsBase)


# make a DSLite server with resources inside
dslite.server <- newDSLiteServer(resources = list(
  gwas1 = resourcer::newResource(name = "gwas1", 
                                 url = "https://github.com/isglobal-brge/brgeUtils/raw/master/data/gwas1.rda",
                                 format = "data.frame")))


dslite.server$assignMethods()

# build login details
builder <- DSI::newDSLoginBuilder()
builder$append(server = "gwas1", url = "dslite.server", 
               resource = "gwas1", driver = "DSLiteDriver")
logindata <- builder$build()# login and assign resources
conns <- datashield.login(logins = logindata, assign = TRUE, symbol = "res")# R data file resource
datashield.assign.expr(conns, symbol = "ES", expr = quote(as.resource.object(res)))


ff <- function(x, snp) {
  ans <- data.frame(as.numeric(snpMatrix1[,snp]), pheno1)
  names(ans)[1] <- snp
  ans
}


dslite.server$assignMethod("ds.snpMatrix", ff)
cally <- "ds.snpMatrix(ES, snp='rs10915297')"
datashield.assign(conns, 'dd', as.symbol(cally))
ds.colnames('dd')

glmSNP('snp', c("affected",  "sex"), 'dd')

model <- "obesity ~ gender"
snps <- "rs10915297"
ds.glmSNP('rs10915297', model = "obesity~gender", data = 'dd')
