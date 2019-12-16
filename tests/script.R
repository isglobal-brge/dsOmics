library(DSOpal)
library(dsBaseClient)
library(dsOmics)

builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "https://opal-test.obiba.org", user = "dsuser", password = "password", resource = "test.GSE66351", driver = "OpalDriver")
builder$append(server = "study2", url = "https://opal-test.obiba.org", user = "dsuser", password = "password", resource = "test.GSE80970", driver = "OpalDriver")
logindata <- builder$build()

# login and assign resources
conns <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "res")

ds.ls(conns)


datashield.assign.expr(conns, symbol = "ES", expr = quote(as.resource.object(res)))

feature <- "cg21477232"
out <- ds.lmFeature(feature, model=casecon~Sex, eSet="ES",
                    connections=conns)


ans <- ds.limma(model=casecon~Sex, eSet = "ES",
                connections = conns)

datashield.logout(conns)


#
# limma
#


#
# VCF
#
