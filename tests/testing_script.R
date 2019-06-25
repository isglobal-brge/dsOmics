library(opal)
library(dsBaseClient)
library(dsStatsClient)
library(dsGraphicsClient)
library(dsModellingClient)
library(dsBetaTestClient)

server <- c("study1","study2")
url <- c("http://192.168.56.100:8080","http://192.168.56.101:8080")
table.m <- c("OMICS.molecular1","OMICS.molecular2")

logindata <- data.frame(server,url,user="administrator",password="datashield_test&",table=table.m)

# login and assign the whole dataset
opals <- datashield.login(logins=logindata,assign=TRUE, symbol='M')
datashield.assign(opals, "P", "OMICS.pheno")

ds.ls()
ds.dim('M')
ds.dim('P')

ds.lmFeature(i=1:10, model='casecon ~ Sex', molecular.data='M', pheno.data='P', datasources=opals)

ds.ls()


