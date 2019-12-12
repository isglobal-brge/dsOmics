source("tests/getClientVCF.R")


# register resolver so that magic happens

resourcer::registerResourceResolver(GDSFileResourceResolver$new())



# use it (method and snpfirstdim are optional)

res <- resourcer::newResource(url = "https://raw.githubusercontent.com/isglobal-brge/scoreInvHap/master/inst/extdata/example.vcf?method=biallelic.only&snpfirstdim=TRUE", format = "VCF2GDS")

client <- resourcer::newResourceClient(res)

gds.f <- client$getValue()

n <- index.gdsn(gds.f, "genotype")
read.gdsn(n, start=c(1, 1), count=c(1, 31))


client$close()
