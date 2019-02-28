suppressPackageStartupMessages({
  library(tools)
  library(AnnotationHub)
  library(yaml)
  library(biomaRt)
  library(Rsamtools)
  library(rtracklayer)
})

LoadReads <- function(inputFilename, genome, annotations){
  
  inputType <- toupper(file_ext(inputFilename))
  
  switch(inputType, 
         BED={
           reads <- BuildReadsGRangesFromBED(inputFilename, genome, annotations)
         },
         BAM={
           reads <- BuildReadsGRangesFromBAM(inputFilename, genome, annotations)
         },
         {
           print_help(opt_parser)
           stop("File type not supported! Please provide a BED or BAM file as input.", call.=FALSE)
         }
  )
  
  return(reads)
  
}

BuildReadsGRangesFromBED <- function(bedFilename, genome, annotations){
  reads <- import.bed(bedFilename)
  genomeSeqInfo <- Seqinfo(seqnames = names(annotations$chrLen),
                           seqlengths = as.numeric(annotations$chrLen),
                           isCircular = rep(NA, length(annotations$chrLen)),
                           genome = rep(genome, length(annotations$chrLen)))
  
  good.seq.levels <- intersect(seqlevels(reads), seqlevels(genomeSeqInfo))
  
  result <- keepSeqlevels(reads, good.seq.levels, pruning.mode="coarse")
  seqlevels(result) <- good.seq.levels
  seqlevels(genomeSeqInfo) <- good.seq.levels
  seqinfo(result) <- genomeSeqInfo
  return(result)
}

BuildReadsGRangesFromBAM <- function(bamFilename, genome, annotations){
  all_fields <- c("rname", "pos", "isize")
  param <- ScanBamParam(what = all_fields, 
                        flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, 
                                           isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                                           isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                                           isNotPassingQualityControls = FALSE))

  bam <- scanBam(bamFilename, param=param)
  aln <- bam[[1]]
  rm(bam)
  
  posStrandReads <- aln$isize > 0
  reads <- GRanges(seqnames=Rle(aln$rname[posStrandReads]),
                   ranges = IRanges(start = aln$pos[posStrandReads], 
                                    width = aln$isize[posStrandReads]))
  
  genomeSeqInfo <- Seqinfo(seqnames = names(annotations$chrLen),
                           seqlengths = as.numeric(annotations$chrLen),
                           isCircular = rep(NA, length(annotations$chrLen)),
                           genome = rep(genome, length(annotations$chrLen)))
  
  good.seq.levels <- intersect(seqlevels(reads), seqlevels(genomeSeqInfo))
  
  result <- keepSeqlevels(reads, good.seq.levels, pruning.mode="coarse")
  seqlevels(result) <- good.seq.levels
  seqlevels(genomeSeqInfo) <- good.seq.levels
  seqinfo(result) <- genomeSeqInfo
  return(result)
}


CleanReads <- function(reads, chrLen, lMin, lMax){

  # Size selection; keep reads with the length between lMin and lMax
  readLength <- width(reads)
  goodReadsInd <- ((readLength >= lMin) & (readLength <= lMax))
  goodReadsLen <- reads[goodReadsInd]
  
  # Compute the raw coverage of these reads
  rawOcc <- coverage(goodReadsLen)

  # Discard regions where the raw coverage is very high -- sequencing artifacts
  # Compute the threshold occupancy for each chromosome
  thresholdFactor <- 10
  chrLabel <- seqlevels(rawOcc)
  chrCoverageThreshold <- lapply(chrLabel, function(chr) thresholdFactor * mean(rawOcc[[chr]]))
  names(chrCoverageThreshold) <- chrLabel

  badRanges <- sapply(names(chrCoverageThreshold),
                      function(chr) {
                        if(is.na(chrCoverageThreshold[[chr]])) {
                          result <- rawOcc[[chr]]
                        } else {
                          result <- reduce(slice(rawOcc[[chr]], lower=chrCoverageThreshold[[chr]], rangesOnly=TRUE) + 500)
                        }
                        return(result)
                      })

  allSeqnames <- names(badRanges)
  badRegions <- sapply(allSeqnames, function(x) {
      badRange <- badRanges[[x]]
      if(length(badRange) == 0) {
        result <- NULL
      } else {
        result <- GRanges(seqnames = x, ranges = badRange, strand = '*')
        seqlevels(result) <- allSeqnames
      }
      return(result)
    })

  badRegionsDF <- lapply(badRegions, function(gr) as.data.frame(gr))
  badRegionsDFAll <- do.call(rbind, badRegionsDF)
  if (ncol(badRegionsDFAll) > 0 & nrow(badRegionsDFAll) > 0) {
    badRegionsAll <- makeGRangesFromDataFrame(badRegionsDFAll)
    badReadIndex <- overlapsAny(goodReadsLen, badRegionsAll, ignore.strand=TRUE)
  }  else {
    # all chr are good; no regions with very high coverage
    badReadIndex <- rep(FALSE, length(goodReadsLen))
  }

  goodReadsLenCov <- goodReadsLen[!badReadIndex]
  
  return(reads)
}

LoadGenomeAnnotation <- function(genome){
  
  chrData <- GetChromosomeLength(genome)
  chrLen <- chrData[, 2]
  names(chrLen) <- chrData[, 1]

  annotations <- GetAnnotations(genome)
  
  result = list(chrLen = chrLen, annotations = annotations)
  return(result)

}

GetChromosomeLength <- function(genome) {
  
  configFilePath <- file.path(configBasePath, "annotation_config.yaml")
  config <- yaml.load_file(configFilePath)
  
  genomeConfig <- config$annotations[[genome]]
  
  url <- genomeConfig$chromosome_size$url
  configFilePath <- genomeConfig$chromosome_size$path
  chrFilter <- as.data.frame(genomeConfig$chromosome_size$filter)
  
  if (! is.null(configFilePath)) {
    filePath <- file.path(annotationsBasePath, configFilePath)
    if(! file.exists(filePath)) {
      if(is.null(url)) {
        # error: no url from where to download chrom.sizes
      } else {
        fileExtension <- file_ext(url)
        destFileName <- paste0("tmp.", fileExtension)
        destfile <- file.path(annotationsBasePath, destFileName)
        content <- download.file(url, destfile, method = "wget", quiet = TRUE)
        compressed <- (fileExtension == "gz")
        if(compressed){
          chrData <- read.table(gzfile(destfile), header = FALSE)
        } else {
          chrData <- read.table(destfile, header = FALSE)
        }
        
        #filter by chr list:
        chrDataToSave <- chrData[order(chrData$V1), c(1:2)]
        if (nrow(chrFilter) > 0 & ncol(chrFilter) > 0) {
          chrDataToSaveFilter <- merge(chrDataToSave, chrFilter, by=c(1))
        } else {
          chrDataToSaveFilter <- chrDataToSave
        }
      
        annotationsFolder <- dirname(filePath)       
        if (! file.exists(annotationsFolder)){           
            dir.create(annotationsFolder)
        }  
          
        # save downloaded data
        write.table(chrDataToSaveFilter, filePath, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
        # clean
        file.remove(destfile)
      }
    }
    result <- read.table(filePath, stringsAsFactors = FALSE, sep="\t")
  } else {
    # error: genomeConfig$chromosome_size$path is empty
  }  
  return(result)
}

GetAnnotations <- function(genome) {
  
  configFilePath <- file.path(configBasePath, "annotation_config.yaml")
  config <- yaml.load_file(configFilePath)
  
  genomeConfig <- config$annotations[[genome]]
  
  mart <- genomeConfig$mart
  file <- genomeConfig$file
  
  if (! is.null(mart)) {
    
    if(is.null(mart$host))  {
      ensembl <- useMart(biomart = mart$biomart, dataset = mart$dataset)  
    } else {
      ensembl <- useMart(host = mart$host,biomart = mart$biomart, dataset = mart$dataset)
    }
    
    annot <- getBM(mart$attributes, mart=ensembl)
    
    # exclude annotations with attr_id == ''
    annot <- annot[annot[mart$attributeId] != '', ]
    
    strand <- annot[[mart$vars$strand]]
    txStart <- annot[[mart$vars$txStart]]
    txEnd <- annot[[mart$vars$txEnd]]
    chrName <- annot[[mart$vars$chrName]]
    
    # check if name begins with "chr", fix it if not:
    someChrName <- tolower(chrName[1])
    if(substr(someChrName, 1, 3) != "chr"){
      chrName <- paste0("chr", chrName)
    } 
    
    TSS <- txStart
    TTS <- txEnd
    TSS[strand == -1] <- txEnd[strand == -1]
    TTS[strand == -1] <- txStart[strand == -1]
    
    annotations <- data.frame(Transcript = annot[mart$attributeId], 
                              Chr = chrName,
                              Strand = strand,
                              TSS = TSS,
                              TTS = TTS,
                              stringsAsFactors = FALSE) 
    
  } else if (! is.null(file)) {
    
    annotationsFilePath <- file.path(annotationsBasePath, file$path)
    annotations <- read.csv(annotationsFilePath, header=TRUE, stringsAsFactors=FALSE)
    
  }  
  
  return(annotations)
}
