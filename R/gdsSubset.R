#' @title Write a subset of data in a GDS file to a new GDS file
#' 
#' @description Function addapted from GWASTools::gdsSubset. Changed the line 31 (of the source file) to have the argument 
#' \code{allow.duplicate = TRUE}
#' 
#' @author Adrienne Stilp
#'
#' @param parent.gds Name of the parent GDS file
#' @param sub.gds Name of the subset GDS file
#' @param sample.include Vector of sampleIDs to include in sub.gds
#' @param snp.include Vector of snpIDs to include in sub.gds
#' @param sub.storage storage type for the subset file; defaults to original storage type
#' @param compress The compression level for variables in a GDS file (see \link{GWASTools::add.gdsn} for options)
#' @param block.size for GDS files stored with scan,snp dimensions, the number of SNPs to read from the parent 
#' file at a time. Ignored for snp,scan dimensions.
#' @param verbose Logical value specifying whether to show progress information.
#' @param allow.fork Logical value specifying whether to enable multiple forks to access the gds file simultaneously.

gdsSubset2 <- function(parent.gds,
                      sub.gds,
                      sample.include=NULL,
                      snp.include=NULL,
                      sub.storage=NULL,
                      compress="LZMA_RA",
                      block.size=5000,
                      verbose=TRUE,
                      allow.fork=FALSE){
  
  # this function only works for gds files having up to two dimensions, named "snp" and "sample"
  gds <- openfn.gds(parent.gds, allow.fork=allow.fork, allow.duplicate = TRUE)
  
  # check that sample.include are all elements of sample.id
  sampID <- read.gdsn(index.gdsn(gds, "sample.id"))
  if (is.null(sample.include)) {
    sample.include <- sampID
  }
  chk <- all(sample.include %in% sampID)
  if (!chk) stop("sample.include elements are not all members of GDS sample.id")
  
  # logical vector for selecting samples
  sampsel <- is.element(sampID, sample.include)
  
  # check that snp.include are all elements of "snp"
  snpID <- read.gdsn(index.gdsn(gds, "snp.id"))
  if (is.null(snp.include)){
    snp.include <- snpID
  }
  chk <- all(snp.include %in% snpID)
  if (!chk) stop("snp.include elements are not all members of GDS snp.id")
  
  # logical vector for selecting snps
  snpsel <- is.element(snpID, snp.include)
  
  # new gds file
  gds.sub <- createfn.gds(sub.gds)
  
  ## TO DO: copy attributes for all
  
  # add sample.id
  node <- add.gdsn(gds.sub, "sample.id", sampID[sampsel], compress=compress, closezip=TRUE, check=TRUE)
  
  # add snp info
  add.gdsn(gds.sub, "snp.id", snpID[snpsel], compress=compress, closezip=TRUE, check=TRUE)
  
  add.gdsn(gds.sub, "snp.position", read.gdsn(index.gdsn(gds, "snp.position"))[snpsel], compress=compress, closezip=TRUE, check=TRUE)
  
  node.parent <- index.gdsn(gds, "snp.chromosome")
  node.sub <- add.gdsn(gds.sub, "snp.chromosome", read.gdsn(node.parent)[snpsel], compress=compress, closezip=TRUE, check=TRUE, storage="uint8")
  attributes <- get.attr.gdsn(node.parent)
  if (length(attributes) > 0){
    for (attribute in names(attributes)){
      put.attr.gdsn(node.sub, attribute, attributes[[attribute]])
    }
  }
  
  dimnames.parent <- ls.gdsn(gds)
  
  # alleles if they exist?
  if ("snp.allele" %in% dimnames.parent){
    add.gdsn(gds.sub, "snp.allele", read.gdsn(index.gdsn(gds, "snp.allele"))[snpsel], compress=compress, closezip=TRUE, check=TRUE)
  }
  # rsID if it exists
  if ("snp.rs.id" %in% dimnames.parent){
    add.gdsn(gds.sub, "snp.rs.id", read.gdsn(index.gdsn(gds, "snp.rs.id"))[snpsel], compress=compress, closezip=TRUE, check=TRUE)
  }
  
  # other dimensions
  dimnames <- dimnames.parent[!(dimnames.parent %in% c("sample.id", "snp.id", "snp.position", "snp.chromosome", "snp.allele", "snp.rs.id"))]
  
  sampID.parent <- sampID
  sampID.sub <- sampID[sampsel]
  
  snpID.parent <- snpID
  snpID.sub <- snpID[snpsel]
  
  # TO DO: check snp.order vs scan.order for snp,scan or scan,snp files
  for (dimname in dimnames){
    
    # node in parent file
    node.parent <- index.gdsn(gds, dimname)
    
    # check dimensions
    desc <- objdesp.gdsn(node.parent)
    
    if (length(desc$dim) != 2) next
    
    if (verbose) message(paste("working on", dimname))
    
    # check attributes of parent node to get array dimensions
    attributes <- get.attr.gdsn(node.parent)
    node.dim <- objdesp.gdsn(node.parent)$dim
    dimType <- NULL
    if ("sample.order" %in% names(attributes)) {
      
      dimType <- "scan,snp" 
      if (!isTRUE(all.equal(node.dim, c(length(sampID.parent), length(snpID.parent))))) stop(paste("dimensions of", dimname, "do not match sample.order"))
      
    } else if ("snp.order" %in% names(attributes)) {
      
      dimType <- "snp,scan"
      if (!isTRUE(all.equal(node.dim, c(length(snpID.parent), length(sampID.parent))))) stop(paste("dimensions of", dimname, "do not match sample.order"))
      
    } else {
      stop("gds node", dimname, "must have 'snp.order' or 'sample.order' as an attribute")
    }
    
    # subset node
    if (dimType == "scan,snp"){
      valdim <- c(sum(sampsel), sum(snpsel))
    } else {
      valdim <- c(sum(snpsel), sum(sampsel))
    }
    
    if (is.null(sub.storage)) {
      storage <- objdesp.gdsn(node.parent)$storage
    } else {
      storage <- sub.storage
    }
    
    node.sub <- add.gdsn(gds.sub, dimname, valdim=valdim, storage=storage)
    
    if (length(attributes) > 0){
      for (attribute in names(attributes)){
        if (attribute == "missing.value" & !is.null(sub.storage) && sub.storage == "bit2"){
          # hard coded missing value
          put.attr.gdsn(node.sub, attribute, 3)
        } else {
          put.attr.gdsn(node.sub, attribute, attributes[[attribute]])
        }
      }
    }
    
    # if genotypeDim is snp,scan then add data sample by sample (possibly with blocks eventually)
    # if genotypeDim is scan,snp then add data snp by snp in blocks of block.size
    if (dimType == "snp,scan"){
      
      for (i in 1:length(sampID.sub)){
        
        if (verbose & i %% 10==0) message(paste(dimname, "- sample", i, "of", sum(sampsel)))
        
        sample <- sampID.sub[i]
        i.parent <- which(sampID.parent %in% sample)
        i.sub <- which(sampID.sub %in% sample)
        
        start.parent <- c(1, i.parent)
        start.sub <- c(1, i.sub)
        count <- c(-1, 1)
        
        dat <- read.gdsn(node.parent, start=start.parent, count=count)
        write.gdsn(node.sub, dat[snpsel], start=start.sub, count=count)
        
      }
    } else if (dimType == "scan,snp"){
      
      nblocks <- ceiling(length(snpID.parent) / block.size)
      i.sub <- 1 # keep track of index in the subset file
      for (i in 1:nblocks){
        if (verbose) message(paste(dimname, "- block", i, "of", nblocks))
        
        # parent snps in this block:
        snps <- snpID.parent[(i-1)*block.size + (1:block.size)]
        snps <- snps[!is.na(snps)] # get rid of NAs from indexing errors
        # snpsel for this block:
        snpsel.block <- which(snps %in% snpID.sub)
        # if there are no snps included, go to the next block
        if (length(which(snps %in% snpID.sub)) == 0) next
        
        # otherwise, continue on: get the parent index to start at
        start.parent <- c(1, (i-1)*block.size + 1)
        count.parent <- c(-1, length(snps))
        
        # read in the block of genotypes
        dat <- read.gdsn(node.parent, start=start.parent, count=count.parent, simplify="none")
        # subset the genotypes for the subset file
        dat <- dat[, snpsel.block, drop=F]
        
        # write them out
        start.sub <- c(1, i.sub)
        count <- c(-1, ncol(dat))
        write.gdsn(node.sub, dat[sampsel, ], start=start.sub, count=count)
        
        # update parameters
        #cnt <- cnt + bsize
        i.sub <- i.sub + ncol(dat)
      }
    }
    
    sync.gds(gds.sub)
    
  }
  
  closefn.gds(gds)
  closefn.gds(gds.sub)
  cleanup.gds(sub.gds, verbose=verbose)
}