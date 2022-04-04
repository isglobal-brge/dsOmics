#' @title Genome-wide association analysis (GWAS) using PLINK 
#' 
#' @description Performs GWAS using PLINK using shell command lines
#' 
#' @param client \code{ssh} SSH resource
#' @param ... \code{character vector} Arguments to pass to the PLINK software
#' 
#' @author Gonzalez, JR.
#'
#' @import readr
#'
#' @export 
#' 
plinkDS <- function(client, ...){
  
  dots <- list(...)
  
  dashedNames <- unlist(lapply(names(dots), function(n) {
    paste0("--", n)
  }))

  plink.command <- gsub("-- ", "--", paste(dashedNames, dots, collapse = " ", sep=" "))
  plink.command <- unlist(strsplit(plink.command, " "))
  
  plink.command <- c(plink.command, "--noweb", "--out")
  tempDir <- client$tempDir()
  command <- c(plink.command, paste0(tempDir, '/out'))
  
  plink <- client$exec('plink1', command)
  
  outs <- client$exec('ls', tempDir)$output
  outs <- outs[-grep(".hh$|.log$|.nof$", outs)]
  
  client$downloadFile(paste0(tempDir, '/', outs))
  
  if (length(outs)==0){
    ans <- plink$error
  }
  else {
    if (length(outs)==1) {
      results <- readr::read_table(outs)
      }
    else {
      results <- c("There are more than 1 table as output")
    }
    ans <- list(results=results, plink.out = plink)
  }
    
  client$removeTempDir()
  client$close()
  
  return(ans)
}
