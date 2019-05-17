suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(pracma)
  library(doParallel)
  library(foreach)
})


ComputeNormalizationFactors <- function(reads)  {

  # Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome 
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  coverageWeight <- lapply(chrLabel, function(chr) 1/mean(occ[[chr]]))
  names(coverageWeight) <- chrLabel
  return(coverageWeight)
}


noCores <- switch(Sys.info()[['sysname']],
                  Windows = 1, # the parallel functions do not work in Windows...
                  Linux   = floor(detectCores(logical = FALSE) / 2), # logical = FALSE does not work in linux..., so we divide by 2
                  Darwin  = detectCores(logical = FALSE)) # Mac

if(noCores > 1) {
  noCores <- min(noCores - 1, 8) # Do not use more than 8 cores (to prevent memory issues)
}

# Parallelized version, much faster (~100x faster on my macbook pro)
ComputeCoverageMatrix <- function(lMin, lMax, beforeRef, afterRef, reads, 
                                  coverageWeight, referenceGRanges, readLength, referenceGenome)
{
  # For human and mouse data, limit the number of cores to 2 (increase this # according to the available memory that is available on your system)
  if (referenceGenome %in% c('mm9', 'mm10', 'hg18', 'hg19', 'hg38')){
    noCores = min(c(2, noCores)) 
  }
  
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  noChr <- length(chrLabel)
  
  registerDoParallel(cores=noCores)
  
  # For each fragment size, compute the corresponding coverage/occupancy
  occMatrixTranspose <- foreach(l=lMin:lMax, .combine='cbind') %dopar% {
    # Keep only the reads with the specific length l
    goodReadsInd <- (readLength == l)
    goodReads <- reads[goodReadsInd]
    
    # Compute average occupancy
    occ <- coverage(goodReads, weight=coverageWeight)
    alignedRegionsTranspose <- AlignRegionsTranspose(occ, referenceGRanges)
    rowMeans(alignedRegionsTranspose)
  }
  occMatrix <- t(occMatrixTranspose)
  rownames(occMatrix) <- c()
  
  stopImplicitCluster()
  #registerDoSEQ()
  
  return(occMatrix)
}

AlignRegionsTranspose <- function(profile, referenceGRanges)
{
  # Create Views with all the referenceGRanges
  chrName <- unique(as.character(seqnames(referenceGRanges)))
  myViews <- Views(profile[chrName], as(referenceGRanges, "IntegerRangesList")[chrName])
  
  alignedProfilesList <- lapply(myViews, function(gr) viewApply(gr, as.vector))
  alignedProfiles <- do.call("cbind", alignedProfilesList)
  
  ## Get the index of referenceGRanges, which were reorganized by as(referenceGRanges, "IntegerRangesList")
  listInd <- split(1:length(referenceGRanges), as.factor(seqnames(referenceGRanges)))
  idx <- do.call("c", listInd)
  
  alignedProfiles <- alignedProfiles[, order(idx)]
  
  ## Flip regions from the Crick strand
  crickInd <- which(strand(referenceGRanges) == "-")
  alignedProfiles[ , crickInd] <- flipud(alignedProfiles[ , crickInd])
  
  return(alignedProfiles)
}

ConstructReferenceGRanges <- function(optSites, annotat, selectedReference, beforeRef, afterRef, genome, optAlign) {
  
  chrLen <- annotat$chrLen

  annotations <- annotat$annotations
  
  if (is.null(optSites)) {
    
    if(selectedReference == "Plus1" & genome != 'sacCer3') {
      print_help(opt_parser)
      stop("Plus1 annotations are available only for the sacCer3 genome. For other genomes, the available options for --reference are: TSS (default), TTS", call.=FALSE)
    } 
    
    referencePos <- annotations[[selectedReference]] 
    referenceChr <- annotations$Chr
    refStrand <- annotations$Strand
    
    validValues <- (! is.nan(referencePos))
    referencePos <- referencePos[validValues]
    referenceChr <- referenceChr[validValues]
    refStrand <- refStrand[validValues]
    
    watson <- (refStrand == 1)
    
  } else {
    
    referenceFilename <- optSites
    referenceFilePath <- file.path(annotationsBasePath, referenceFilename)
    sites <- import.bed(referenceFilePath)
    
    # Make sure the chromosome names are following the 'UCSC' convention
    seqlevelsStyle(sites) <- 'UCSC'
    
    switch(optAlign, 
           "center"    ={ fixAlign = "center" },
           "fivePrime" ={ fixAlign = "start" },
           "threePrime"={ fixAlign = "end" })
    
    referenceGRanges <- resize(sites, width = 1, fix = fixAlign)
    
    referencePos <- start(referenceGRanges)
    referenceChr <- seqnames(referenceGRanges)
    refStrand <- strand(referenceGRanges)
    
    watson <- as.vector(strand(refStrand) != '-')      # ('+' or '*')
  }
  
  leftEdge <- referencePos - beforeRef
  rightEdge <- referencePos + afterRef
  
  leftEdgeCrick <- referencePos - afterRef
  rightEdgeCrick <- referencePos + beforeRef
  
  leftEdge[!watson] <- leftEdgeCrick[!watson]
  rightEdge[!watson] <- rightEdgeCrick[!watson]
  
  wholeChr <- GRanges(seqnames=names(chrLen), IRanges(start=rep(1, length(chrLen)), end=chrLen))  
  
  referenceGRanges <- GRanges(seqnames=referenceChr, IRanges(start=leftEdge, end=rightEdge), strand=refStrand)
  result <- subsetByOverlaps(referenceGRanges, wholeChr, type="within", ignore.strand=TRUE)
  
  return(result)
  
}


CalculatePlotData <- function(params, reads, referenceGRanges) {
  
  # Compute the histogram
  readLength <- width(reads) 
  h <- hist(readLength, breaks=seq(from = 0.5, to = max(readLength)+0.5, by = 1), plot=FALSE)
  lengthHist <- 100 * h$density[params$lMin:params$lMax]
  
  totalNoReads <- length(reads)
  
  resized_reads <- switch(params$plotType, 
                      "OCC"            = reads,
                      "DYADS"          = resize(reads, width = 1, fix = "center"),
                      "FIVEPRIME_ENDS" = resize(reads, width = 1, fix = "start"),
                      "THREEPRIME_ENDS"= resize(reads, width = 1, fix = "end"))
  
  coverageWeight <- ComputeNormalizationFactors(resized_reads)
  
  occMatrix <- ComputeCoverageMatrix(params$lMin, params$lMax, params$beforeRef, params$afterRef, 
                         resized_reads, coverageWeight, referenceGRanges, readLength, params$genome)
  
  outputFilePath <- GetOutputMatrixFilePath(params$plotType, params$referencePointsBed, 
                                            params$reference, params$siteLabel, 
                                            params$lMin, params$lMax, params$sampleName)
  
  # Parameters to be saved:
  lMin <- params$lMin
  lMax <- params$lMax
  beforeRef <- params$beforeRef
  afterRef <- params$afterRef
  sampleName <- params$sampleName
  siteLabel <- params$siteLabel
  save(occMatrix, lMin, lMax, beforeRef, afterRef, totalNoReads, lengthHist, sampleName, siteLabel, file=outputFilePath)
  
}

