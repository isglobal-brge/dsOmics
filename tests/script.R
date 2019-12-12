library(DSOpal)
library(dsBaseClient)
library(dsOmics)
library(parallel)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", user = "dsuser", password = "password", resource = "test.GSE66351", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-test.obiba.org", user = "dsuser", password = "password", resource = "test.GSE80970", driver = "OpalDriver")
logindata <- builder$build()

# login and assign resources
conns <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "res")
# R data file resource
ds.class("res")


datashield.assign.expr(conns, symbol = "ES", expr = quote(as.resource.object(res)))

feature <- "cg21477232"
out <- ds.lmFeature(feature, model=casecon~Sex, eSets="ES")


datashield.logout(conns)


model <- "casecon ~ Sex"
mt <- as.formula(model)
vars <- all.vars(mt)

dep <- paste(vars, collapse="+")

feature <- "cg21477232"
mm <- as.formula(paste(feature, "~", dep))
mod <- ds.glm(mm, family='gaussian', data='D')

lmFeature(feature, model, 'D', vars)

ans <- t(as.data.frame(lapply(features, lmFeatureGLM, model=model, 
                                eSet.data='D', vars=vars))) 



datashield.assign.expr(conns, symbol = "ES", 
                       expr = quote(as.resource.object(res)))

out <- ds.lmFeature("cg21477232", model=model, eSet.data = "D", eSet="ES")





###########

library(resourcer)
res2 <- resourcer::newResource(name = "GSE66351", url = "https://github.com/epigeny/dsOmics/raw/master/data/GSE66351.Rdata", format = "ExpressionSet")
client <- resourcer::newResourceClient(res2)
# coerce to a data.frame
dim(client$asDataFrame())
# coerce to a dplyr's tbl
dim(client$asTbl())
# the raw ExpressionSet object
dim(client$getValue())

datashield.logout(client)



res <- resourcer::newResource(name = "exposome", 
                              url = "https://github.com/isglobal-brge/rexposome/blob/master/data/exposome.rda", 
                              format = "ExposomeSet")
client <- resourcer::newResourceClient(res)

