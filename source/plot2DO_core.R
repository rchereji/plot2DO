suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

# NUM_CORES <- detectCores(all.tests = FALSE, logical = TRUE) - 1
# print(NUM_CORES)

ComputeNormalizationFactors <- function(reads)  {

  # Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome 
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  coverageWeight <- lapply(chrLabel, function(chr) 1/mean(occ[[chr]]))
  # coverageWeight <- mclapply(chrLabel, function(chr) 1/mean(occ[[chr]]), mc.cores = NUM_CORES)
  names(coverageWeight) <- chrLabel
  return(coverageWeight)
}


# debug:
# lMin <- params$lMin
# lMax <- params$lMax
# beforeRef <- params$beforeRef
# afterRef <- params$afterRef
# coverageWeight <- ComputeNormalizationFactors(reads)
# readLength <- width(reads) 
ComputeCoverageMatrix <- function(lMin, lMax, beforeRef, afterRef, reads, 
                      coverageWeight, referenceGRanges, readLength)
{
  #resized reads to fix mm9 and hg19:
  # TODO: ver con Razvan si este filtro esta bien para resolver el problema de los levels
  # referenceGRangesChrNames <- as.character(seqnames(referenceGRanges))
  # resizedReads <- reads[as.character(seqnames(reads)) %in% referenceGRangesChrNames]
  # seqlevels(resizedReads) <- seqlevels(referenceGRanges)
  
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  noChr <- length(chrLabel)
  
  rowCount <- lMax - lMin + 1
  colCount <- 1 + beforeRef + afterRef
  occMatrix <- matrix(data = 0, nrow = rowCount, ncol = colCount)
  
  # For each fragment size, compute the corresponding coverage/occupancy
  for(l in lMin:lMax) {
    # Keep only the reads with the specific length l
    goodReadsInd <- (readLength == l)
    goodReads <- reads[goodReadsInd]
    
    # Compute average occupancy
    occ <- coverage(goodReads, weight=coverageWeight)
    alignedRegions <- AlignRegions(occ, referenceGRanges)
    occMatrix[l-lMin+1,] <- colMeans(alignedRegions)
  }
  
  return(occMatrix)
}


# profile <- occ

# Something is fucked up in this function! I will replace it with the old version

# AlignRegions = function(profile, referenceGRanges)
# {
#   # Obtain Views for all GRanges that we wish to align
#   myViews = Views(profile, referenceGRanges)
#   
#   # Convert the RleViewsList (myViews) into a matrix
#   alignedProfilesList = lapply(myViews, function(gr) t(viewApply(gr, as.vector)))
#   alignedProfiles = do.call("rbind", alignedProfilesList)
#   
#   ## Flip the rows corresponding to GRanges from the Crick strand
#   CrickInd = which(as.character(strand(referenceGRanges)) == "-")
#   alignedProfiles[CrickInd, ] = alignedProfiles[CrickInd, ncol(alignedProfiles):1]
#   
#   return(alignedProfiles)
# }

AlignRegions = function(profile, referenceGRanges)
{
  # Create Views with all the referenceGRanges
  chrName = unique(as.character(seqnames(referenceGRanges)))
  myViews = Views(profile[chrName], as(referenceGRanges, "IntegerRangesList")[chrName])
  alignedProfilesList = lapply(myViews, function(gr) t(viewApply(gr, as.vector)))
  alignedProfiles = do.call("rbind", alignedProfilesList)
  
  ## Get the index of referenceGRanges, which were reorganized by as(referenceGRanges, "IntegerRangesList")
  listInd = split(1:length(referenceGRanges), as.factor(seqnames(referenceGRanges)))
  idx = do.call("c", listInd)
  
  rownames(alignedProfiles) = idx
  alignedProfiles = alignedProfiles[order(idx),]
  
  ## Flip regions from the Crick strand
  CrickInd = which(as.character(strand(referenceGRanges)) == "-")
  alignedProfiles[CrickInd,] = alignedProfiles[CrickInd, ncol(alignedProfiles):1]
  return(alignedProfiles)
}


# optSites <- opt$sites
# annotat <- annotations
# optAlign <- opt$align
ConstructReferenceGRanges <- function(optSites, annotat, selectedReference, beforeRef, afterRef, genome, optAlign) {
  
  chrLen <- annotat$chrLen
  # annotations.tmp <- annotat$annotations
  # armar esto solo para los cromosomas de interes!!
  # filter with chr names - fix levels problem:
  # TODO paula: descomentar esto
  # annotations <- annotations.tmp[annotations.tmp$Chr %in% names(chrLen), ]
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
    
    referenceFilename = optSites
    referenceFilePath <- file.path(annotationsBasePath, referenceFilename)
    sites = import.bed(referenceFilePath)
    
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
  h <- hist(readLength, breaks=seq(from = 0.5, to = 1000.5, by = 1), plot=FALSE)
  lengthHist <- 100 * h$density[params$lMin:params$lMax]
  
  totalNoReads <- length(reads)
  
  resized_reads <- switch(params$plotType, 
                      "OCC"            = reads,
                      "DYADS"          = resize(reads, width = 1, fix = "center"),
                      "FIVEPRIME_ENDS" = resize(reads, width = 1, fix = "start"),
                      "THREEPRIME_ENDS"= resize(reads, width = 1, fix = "end"))
  
  coverageWeight <- ComputeNormalizationFactors(resized_reads)
  
  occMatrix <- ComputeCoverageMatrix(params$lMin, params$lMax, params$beforeRef, params$afterRef, 
                         resized_reads, coverageWeight, referenceGRanges, readLength)
  
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

