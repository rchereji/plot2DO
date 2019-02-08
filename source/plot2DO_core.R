suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

# NUM_CORES <- detectCores(all.tests = FALSE, logical = TRUE) - 1
# print(NUM_CORES)

Normalization <- function(reads)  {
  #################
  # Normalization #
  #################
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
Occupancy <- function(lMin, lMax, beforeRef, afterRef, reads, 
                      coverageWeight, referenceGRanges, readLength)
{
  #resized reads to fix mm9 and hg19:
  # TODO: ver con Razvan si este filtro esta bien para resolver el problema de los levels
  # referenceGRangesChrNames <- as.character(seqnames(referenceGRanges))
  # resizedReads <- reads[as.character(seqnames(reads)) %in% referenceGRangesChrNames]
  # seqlevels(resizedReads) <- seqlevels(referenceGRanges)
  #seqnames(reads.sub)
  
  occ <- coverage(reads)
  chrLabel <- seqlevels(occ)
  noChr <- length(chrLabel)
  
  rowCount <- lMax - lMin + 1
  colCount <- 1 + beforeRef + afterRef
  occMatrix <- matrix(data = 0, nrow = rowCount, ncol = colCount)
  
  # For each fragment size, compute the corresponding coverage/occupancy
  for(l in lMin:lMax) {
    # print(l)
    # Keep only the reads with the specific length l
    goodReadsInd <- (readLength == l)
    goodReadsInd[1:noChr] <- TRUE # extra 1bp reads at the ends all chromosomes
    goodReads <- reads[goodReadsInd]
    
    if (length(goodReads) > noChr) {
      # Compute average occupancy
      occ <- coverage(goodReads, weight=coverageWeight)
      # data <- PrepareAlignRegions(occ, referenceGRanges)
      alignedRegions <- AlignRegions(occ, referenceGRanges)
      # alignedRegions <- AlignRegions(data$profile, data$referenceGRanges)
      occMatrix[l-lMin+1,] <- colMeans(alignedRegions)
    }
  }
  
  return(occMatrix)
}

# profile <- occ
AlignRegions <- function(profile, referenceGRanges)
{
  # problema: referenceGRanges tiene 16 cromosomas y
  # profile tiene 16 + M
  # deben tener los mismos ???
  #referenceGRangesSeqInfo <- seqinfo(referenceGRanges)
  #referenceGRangesChrs <- seqnames(referenceGRangesSeqInfo)

  # TODO: ver esto con razvan. esto es para resolver el problema que aparece en views
  
  chrNames <- seqlevels(referenceGRanges)
  chrNames.util <- intersect(chrNames, names(profile))

  profileChr <- profile[chrNames.util]
  seqlevels(profileChr) <- chrNames.util #esto arregla el error de views que aparece cuando no coinciden los levels

  # TODO: paula limpiar referenceGRanges antes de pasarlo como argumento porque esta en un ciclo y es siempre igual
  referenceGRangesChr <- referenceGRanges[seqnames(referenceGRanges) %in% chrNames.util]
  seqlevels(referenceGRangesChr) <- chrNames.util

  # TODO: fin de ver esto con razvan
  
  # Obtain Views for all GRanges that we wish to align
  myViews <- Views(profileChr, referenceGRangesChr)
  # myViews <- Views(profile, referenceGRanges)

  # Convert the RleViewsList (myViews) into a matrix
  alignedProfilesList <- lapply(myViews, function(gr) {
  # alignedProfilesList <- mclapply(myViews, function(gr) {
    #limpio los que no tienen datos, dan error en rbind:
    if(length(gr) > 0){
      result <- t(viewApply(gr, as.vector))
    } else {
      # print("else")
      result <- NULL
    }
    return(result)
  })
  # }, mc.cores = NUM_CORES)

  alignedProfiles <- do.call("rbind", alignedProfilesList)
  ## Flip the rows corresponding to GRanges from the Crick strand
  crickInd <- which(as.character(strand(referenceGRanges)) == "-")
  if(length(crickInd) > 0) {
    alignedProfiles[crickInd,] <- alignedProfiles[crickInd, ncol(alignedProfiles):1]
  }
  return(alignedProfiles)

}

# optSites <- opt$sites
# annotat <- annotations
# optAlign <- opt$align
Align <- function(optSites, annotat, selectedReference, beforeRef, afterRef, genome, optAlign) {
  
  chrLen <- annotat$chrLen
  # annotations.tmp <- annotat$annotations
  # armar esto solo para los cromosomas de interes!!
  # filter with chr names - fix levels problem:
  # TODO paula: descomentar esto
  # annotations <- annotations.tmp[annotations.tmp$Chr %in% names(chrLen), ]
  annotations <- annotat$annotations
  
  if (is.null(optSites)) {
    
    # exit program ?
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
           "center"={ fixAlign = "center" },
           "fivePrime"={ fixAlign = "start" },
           "threePrime"={ fixAlign = "end" })
    
    referenceGRanges <- resize(sites, width = 1, fix = fixAlign)
    
    watson <- as.vector(strand(referenceGRanges) != '-')      # ('+' or '*')
    
    referencePos <- start(referenceGRanges)
    
    referenceChr <- seqnames(referenceGRanges)
    
    refStrand <- strand(referenceGRanges)
    
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
  # readLength is not recomputed for any type of plot (including those with resized reads)
  h <- hist(readLength, breaks=seq(from = 0.5, to = 1000.5, by = 1), plot=FALSE)
  lengthHist <- 100 * h$density[params$lMin:params$lMax]
  
  totalNoReads <- length(reads)
  
  if ("OCC" %in% params$plotType){
    resized_reads <- reads
  }
  if ("DYADS" %in% params$plotType){
    resized_reads <- resize(reads, width = 1, fix = "center")
  }
  if ("FIVEPRIME_ENDS" %in% params$plotType){
    resized_reads = resize(reads, width = 1, fix = "start")
  }
  if ("THREEPRIME_ENDS" %in% params$plotType){
    resized_reads = resize(reads, width = 1, fix = "end")
  }  
  
  coverageWeight <- Normalization(resized_reads)
  occMatrix <- Occupancy(params$lMin, params$lMax, params$beforeRef, params$afterRef, 
                         resized_reads, coverageWeight, referenceGRanges, readLength)
  outputFilePath <- GetOutputMatrixFilePath(params$plotType, params$referencePointsBed, 
                                            params$reference, params$siteLabel, 
                                            params$lMin, params$lMax, params$sampleName)
  
  # required to be saved:
  lMin <- params$lMin
  lMax <- params$lMax
  beforeRef <- params$beforeRef
  afterRef <- params$afterRef
  sampleName <- params$sampleName
  siteLabel <- params$siteLabel
  save(occMatrix, lMin, lMax, beforeRef, afterRef, totalNoReads, lengthHist, sampleName, siteLabel, file=outputFilePath)
  
}

