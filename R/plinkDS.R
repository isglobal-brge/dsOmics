#' @title Genome-wide association analysis (GWAS) using PLINK 
#' 
#' @description Performs GWAS using PLINK using shell command lines
#' 
#' @param client ...
#' @param ... 
#' 
#' @author Gonzalez, JR.
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
  
  client$downloadFile(paste0(tempDir, '/out.*'))  
  
  outs <- client$exec('ls', tempDir)$output
  outs <- outs[-grep(".hh$|.log$|.nof$", outs)]
  
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
