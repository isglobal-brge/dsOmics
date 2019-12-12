#### Create new
library(resourcer)
library(DSLite)
library(dsBaseClient)
library(dsBase)
library(dsOmics)


# make a DSLite server with resources inside
dslite.server <- newDSLiteServer(resources = list(
  GSE66351 = resourcer::newResource(name = "GSE66351", url = "https://github.com/epigeny/dsOmics/raw/master/data/GSE66351.Rdata", format = "ExpressionSet"),
  GSE80970 = resourcer::newResource(name = "GSE80970", url = "https://github.com/epigeny/dsOmics/raw/master/data/GSE80970.Rdata", format = "ExpressionSet")
))

dslite.server$assignMethods()

# build login details
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "dslite.server", resource = "GSE66351", driver = "DSLiteDriver")
builder$append(server = "study2", url = "dslite.server", resource = "GSE80970", driver = "DSLiteDriver")
logindata <- builder$build()# login and assign resources
conns <- datashield.login(logins = logindata, assign = TRUE, symbol = "res")# R data file resource


datashield.assign.expr(conns, symbol = "ES", expr = quote(as.resource.object(res)))

dslite.server$assignMethods()


dslite.server$assignMethod("selFeature", selFeature)
dslite.server$assignMethod("cellCounts", cellCounts)
dslite.server$assignMethod("exprsDS", function(x) Biobase::exprs(x))
dslite.server$assignMethod("pDataDS", function(x) Biobase::pData(x))

dslite.server$aggregateMethod("featureNamesDS", function(x) Biobase::featureNames(x))



feature <- "cg21477232"
features <- c("cg21477232", "cg21477232")
vars <- c("casecon", "Sex")

lmFeature(feature, vars, "ES")
out <- ds.lmFeature(features, model=casecon~Sex, eSets="ES")



cellCounts <- function(x, cellTypeRef="blood gse35069 complete"){
  try( meffil::meffil.estimate.cell.counts.from.betas(Biobase::exprs(x),
                                                      cell.type.reference = cellTypeRef, verbose = FALSE), TRUE)
}





# ---------------------



ds.colnames("D")




datashield.aggregate(conns, "ds.exprs(ES)")

datashield.assign(conns, 'dd', "ds.exprs(ES)")
ds.dim("dd")

dslite.server$aggregateMethod("selFeature", function(x, feature, vars) {
  data.frame(x[feature, ])[, c(feature, vars)]
})

datashield.aggregate(conns, as.symbol('selFeature(ES, feature="cg21477232", vars="casecon")'))

cally <- "selFeature(ES, feature='cg21477232')"
datashield.assign(conns, 'dd', as.symbol(cally))
ds.dim('dd')

model <- "casecon ~ Sex"
mt <- as.formula(model)
vars <- all.vars(mt)

data <- 'ES'

lmFeature("cg21477232", vars, 'ES')





# and use it
datashield.aggregate(conns, "selData(ES, feature='cg21477232')")

cally <- "selData(ES, feature='cg21477232')"
datashield.assign(conns, 'dd', as.symbol(cally))
ds.dim("dd")

mod <- ds.glm(cg21477232~casecon, family='gaussian', data='dd', datasources = conns)



ds.ls(conns)



##########


builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", 
               user = "dsuser", password = "password", 
               resource = "omicBioC.gwas1", driver = "OpalDriver")
logindata <- builder$build()

conns <- DSI::datashield.login(logins = logindata, 
                               assign = TRUE, symbol = "res")


datashield.assign.expr(conns, symbol = "ES", 
                       expr = quote(as.resource.object(res)))

ds.dataFrame("ES", newobj = "snps")

datashield.assign(conns, 'pheno', as.symbol("res[[1]]"))

datashield.logout(conns)

