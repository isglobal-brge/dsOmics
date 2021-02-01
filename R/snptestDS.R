#' @title Interface to run SNPtest commands on a ssh connection (with SNPtest 2.5.2 installed)
#'
#' @param client \code{ssh} SSH resource
#' @param ... \code{character vector} Arguments to pass to the snptest_v2.5.2 software
#'
#' @return List containing: \cr
#' - Results: Table of results (typical SNPtest output file) \cr
#' - Output: Console output of the ssh query
#' @export
#'

snptestDS <- function(client, ...){
  
  dots <- list(...)
  arguments <- unlist(dots)
  
  # Add temp folder and set it for output
  tempDir <- client$tempDir()
  arguments <- c(arguments, "-o", paste0(tempDir, "/ex.out"))
  
  # snptest execution command
  snptest.command <- "snptest_v2.5.2"
  
  # Run snptest
  snptest <- client$exec(snptest.command, arguments)
  
  # Get output file
  client$downloadFile(paste0(tempDir, '/ex.out'))  
  
  # Checkoutputs
  outs <- client$exec('ls', tempDir)$output
  
  if (length(outs)==0){ # No outputs means errors
    ans <- snptest$error
  }
  else {
    if (length(outs)==1) {
      results <- readr::read_table(outs)
    }
    else {
      results <- c("There are more than 1 table as output")
    }
    ans <- list(results=results, plink.out = snptest)
  }
  
  client$removeTempDir()
  client$close()
  
  return(ans)
}
  