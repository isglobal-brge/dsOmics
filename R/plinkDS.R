#' @title Genome-wide association analysis (GWAS) using PLINK 
#' 
#' @description Performs GWAS using PLINK using shell command lines
#' 
#' @param 
#' @author Gonzalez, JR.
#'
#' @import readr
#' @export 
#' 
plinkDS <- function(client, plink.command,   ...){
  
  plink.command <- unlist(strsplit(plink.command, " "))
  plink.command <- c(plink.command, "--noweb", "--out")
  base.out <- basename(tempdir())
  file.created <- paste0('/tmp/', base.out, '-out')
  command <- c(plink.command, file.created)
  plink <- client$exec('plink1', command)
  
  client$downloadFile('/tmp/*-out.*')  
  
  outs <- client$exec('ls', '/tmp')$output
  outs <- outs[grep(base.out, outs)]
  outs <- outs[-grep(".hh|.log|.nof", outs)]
  
  nn <- sapply(strsplit(outs, "\\."), "[", 2)
  results <- list()
  for (i in 1:length(outs)){
    results[[i]] <- readr::read_table(outs[i])
  }
  names(results) <- nn
  
  client$exec('rm', paste0('/tmp/',  base.out, '-out.*'))
  
  ans <- list(results=results, plink.out = plink)
  return(ans)
}
