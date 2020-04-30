#' @title Genome-wide association analysis (GWAS) using PLINK 
#' 
#' @description Performs GWAS using PLINK using shell command lines
#' 
#' @param client ...
#' @param ... 
#' 
#' @author Gonzalez, JR.
#'
#' @import readr
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
  
  else{
    nn <- sapply(strsplit(outs, "\\."), "[", 2)
    results <- list()
    for (i in 1:length(outs)){
      results[[i]] <- readr::read_table(outs[i])
    }
    names(results) <- nn
    ans <- list(results=results, plink.out = plink)
  }
  
  client$removeTempDir()
  
  return(ans)
}
