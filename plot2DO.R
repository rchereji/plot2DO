#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Dataset file name [options: BED or BAM format]"),
  make_option(c("-r", "--reference"), type="character", default="TSS", 
              help="Reference points to align [options: TSS (default), TTS, Plus1]"),
  make_option(c("-l", "--minLength"), type="integer", default=50, 
              help="The smallest DNA fragment to be considered [default = %default]"),
  make_option(c("-L", "--maxLength"), type="integer", default=200, 
              help="The largest DNA fragment to be considered [default = %default]"),
  make_option(c("-u", "--upstream"), type="integer", default=1000, 
              help="Length of the upstream region that will be plotted [default = %default]"),
  make_option(c("-d", "--downstream"), type="integer", default=1000, 
              help="Length of the downstream region that will be plotted [default = %default]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least the dataset file name must be supplied.", call.=FALSE)
}

##################
# Initialization #
##################
# Type of alignments, i.e. reference points
selectedReference = opt$reference

# Size selection parameters: specify the interval of lengths to be analyzed
Lmin = opt$minLength   # the smallest DNA fragment to be considered
Lmax = opt$maxLength  # the largest DNA fragment to be considered

# Window selection parameters
beforeRef = opt$upstream  # length of the upstream region that will be plotted
afterRef = opt$downstream   # length of the downstream region that will be plotted

# Load the necessary R packages
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(caTools))
suppressPackageStartupMessages(library(colorRamps))

#########################################
# Import the paired-end sequencing data #
#########################################
# Data file name
inputFilename = opt$file
library(tools)
inputType = toupper(file_ext(inputFilename))
switch(inputType, 
       BED={
         sample.name = sub(".bed", "", inputFilename)
         reads = import.bed(inputFilename)
       },
       BAM={
         sample.name = sub(".bam", "", inputFilename)
         suppressPackageStartupMessages(library(GenomicAlignments))
         param = ScanBamParam(flag=scanBamFlag(isProperPair=TRUE))
         ga = readGAlignmentPairs(inputFilename, use.names=FALSE, param=param)
         reads = granges(ga)
       },
       {
         print_help(opt_parser)
         stop("File type not supported! Please provide a BED or BAM file as input.", call.=FALSE)
       }
)

# Discard the reads from the rDNA region, chrXII:451000-469000
rDNAregion = GRanges(seqnames = "chrXII",
                     ranges = IRanges(start=451000, end=469000),
                     strand = "*")
rDNAInd = overlapsAny(reads, rDNAregion, ignore.strand=TRUE)
reads = reads[!rDNAInd]

# Discard the reads from the yeast mitochondrial DNA
reads = reads[seqnames(reads) != 'chrM']


#########################################
# Compute the fragment length histogram #
#########################################
# Compute the histogram
readLength = width(reads)
h = hist(readLength, breaks=seq(from = 0.5, to = 1000.5, by = 1), plot=FALSE)
LengthHist = 100*h$density[Lmin:Lmax]


##################
# Size selection #
##################
# Eliminate the reads that are shorter than Lmin or longer than Lmax
goodReadsInd = ((readLength >= Lmin) & (readLength <= Lmax))
reads = reads[goodReadsInd]
readLength = width(reads)
TotalNoReads = length(reads)

# Make sure that each chromosome has at least a read on it (add an extra read of 1bp at the end of each chromosome)
chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
noChr = length(chrLen)
extraReads = GRanges(seqnames = names(chrLen),
                     ranges = IRanges(start = chrLen, width = 1),
                     strand = rep("+",noChr),
                     seqlengths = chrLen)
reads = c(reads, extraReads)


#################
# Normalization #
#################
# Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome 
Occ = coverage(reads)
chrLabel = seqlevels(Occ)
noChr = length(chrLabel)
coverageWeight = list()
for(chr in chrLabel)
{
  coverageWeight[[chr]] = 1/mean(Occ[[chr]])
}


##################################
## Compute 2D occupancy matrices #
##################################
dir.create("2D_occupancy", showWarnings = FALSE, recursive = TRUE)

# Load annotations
sacCer3transcripts = read.csv("sacCer3_annotations.csv", header=TRUE, stringsAsFactors=FALSE)

## Initialize the 2D occ. matrix
Occ_matrix = matrix(data = 0, nrow = Lmax-Lmin+1, ncol = 1+beforeRef+afterRef)

switch(selectedReference, 
       TSS={
         ReferencePos = sacCer3transcripts$TSS
         ReferenceChr = sacCer3transcripts$Chr
         RefStrand = sacCer3transcripts$Strand
         
         Watson = (RefStrand == 1)
         leftEdge = ReferencePos-beforeRef
         rightEdge = ReferencePos+afterRef
         
         leftEdgeCrick = ReferencePos-afterRef
         rightEdgeCrick = ReferencePos+beforeRef
         
         leftEdge[!Watson] = leftEdgeCrick[!Watson]
         rightEdge[!Watson] = rightEdgeCrick[!Watson]
         
         ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                                     IRanges(start=leftEdge,
                                             end=rightEdge),
                                     strand=RefStrand)
         
         # Construct GRanges for the entire chromosomes...
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         wholeChr = GRanges(seqnames=names(chrLen),
                            IRanges(start=rep(1, length(chrLen)),
                                    end=chrLen))
         # ... and remove the windows that fall outside of the chromosome edges
         ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                                 type="within", ignore.strand=TRUE)
       },
       TTS={
         ReferencePos = sacCer3transcripts$TTS
         ReferenceChr = sacCer3transcripts$Chr
         RefStrand = sacCer3transcripts$Strand
         # Keep only TTSs that are away from other TSSs
         validValues = (! is.nan(ReferencePos)) & ((sacCer3transcripts$NDR3pLength > 1000) | (sacCer3transcripts$NDR3pQualifier == "convergent"))
         ReferencePos = ReferencePos[validValues]
         ReferenceChr = ReferenceChr[validValues]
         RefStrand = RefStrand[validValues]
         
         Watson = (RefStrand == 1)
         leftEdge = ReferencePos-beforeRef
         rightEdge = ReferencePos+afterRef
         
         leftEdgeCrick = ReferencePos-afterRef
         rightEdgeCrick = ReferencePos+beforeRef
         
         leftEdge[!Watson] = leftEdgeCrick[!Watson]
         rightEdge[!Watson] = rightEdgeCrick[!Watson]
         
         ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                                     IRanges(start=leftEdge,
                                             end=rightEdge),
                                     strand=RefStrand)
         
         # Construct GRanges for the entire chromosomes...
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         wholeChr = GRanges(seqnames=names(chrLen),
                            IRanges(start=rep(1, length(chrLen)),
                                    end=chrLen))
         # ... and remove the windows that fall outside of the chromosome edges
         ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                                 type="within", ignore.strand=TRUE)  
       },
       Plus1={
         ReferencePos = sacCer3transcripts$Plus1
         ReferenceChr = sacCer3transcripts$Chr
         RefStrand = sacCer3transcripts$Strand
         validValues = ! is.nan(ReferencePos)
         ReferencePos = ReferencePos[validValues]
         ReferenceChr = ReferenceChr[validValues]
         RefStrand = RefStrand[validValues]
         
         Watson = (RefStrand == 1)
         leftEdge = ReferencePos-beforeRef
         rightEdge = ReferencePos+afterRef
         
         leftEdgeCrick = ReferencePos-afterRef
         rightEdgeCrick = ReferencePos+beforeRef
         
         leftEdge[!Watson] = leftEdgeCrick[!Watson]
         rightEdge[!Watson] = rightEdgeCrick[!Watson]
         
         ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                                     IRanges(start=leftEdge,
                                             end=rightEdge),
                                     strand=RefStrand)
         
         # Construct GRanges for the entire chromosomes...
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         wholeChr = GRanges(seqnames=names(chrLen),
                            IRanges(start=rep(1, length(chrLen)),
                                    end=chrLen))
         # ... and remove the windows that fall outside of the chromosome edges
         ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                                   type="within", ignore.strand=TRUE)   
       }
)

######################################
# Function to align multiple regions #
######################################
# Input arguments:
#   Profile           - The profile that needs to be aligned, e.g. coverage/occupancy profile
#   ReferenceGRanges  - GRanges with the windows centered on the reference points
#
# Output:
#   AlignedProfiles   - a matrix with each row corresponding to an aligned locus

AlignRegions = function(Profile, ReferenceGRanges)
{
  # Create Views with all the ReferenceGRanges
  chrName = unique(as.character(seqnames(ReferenceGRanges)))
  myViews = Views(Profile[chrName], as(ReferenceGRanges, "RangesList")[chrName])
  AlignedProfilesList = lapply(myViews, function(gr) t(viewApply(gr, as.vector)))
  AlignedProfiles = do.call("rbind", AlignedProfilesList)
  
  ## Get the index of ReferenceGRanges, which were reorganized by as(ReferenceGRanges, "RangesList")
  listInd = split(1:length(ReferenceGRanges), as.factor(seqnames(ReferenceGRanges)))
  idx = do.call("c", listInd)
  
  rownames(AlignedProfiles) = idx
  AlignedProfiles = AlignedProfiles[order(idx),]
  
  ## Flip regions from the Crick strand
  CrickInd = which(as.character(strand(ReferenceGRanges)) == "-")
  AlignedProfiles[CrickInd,] = AlignedProfiles[CrickInd, ncol(AlignedProfiles):1]
  return(AlignedProfiles)
}

# For each fragment size, compute the corresponding coverage/occupancy
for(L in Lmin:Lmax) {
  
  # Keep only the reads with the specific length L
  goodReadsInd = (readLength == L)
  goodReadsInd[1:noChr] = TRUE # extra 1bp reads at the ends all chromosomes
  goodReads = reads[goodReadsInd]
  
  if (length(goodReads) > noChr) {
    # Compute average occupancy
    Occ = coverage(goodReads, weight=coverageWeight)
    Occ_matrix[L-Lmin+1,] = colMeans(AlignRegions(Occ, ReferenceGRanges))
  }
}

# Save the 2D Occupancy matrix
switch(selectedReference, 
       TSS={
         save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, file=paste("2D_occupancy/Occ_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
       },
       TTS={
         save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, file=paste("2D_occupancy/Occ_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
       },
       Plus1={
         save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, file=paste("2D_occupancy/Occ_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
       }
)

################################
# Plot the 2D Occupancy matrix #
################################
# Function for drawing the color scale
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("Must have one more break than colour!")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3) # adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

# Plot the figure
switch(selectedReference, 
       TSS={
         Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors
         
         pdf(paste("2D_occupancy/Occ_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
         layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
         
         # Empty plot
         plot(0, type='n', axes=FALSE, ann=FALSE)
         
         # Add scale
         par(mar=c(5,6,2,0.5))
         image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
         axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
         mtext("Relative coverage (%)", cex=1, side=2, padj=-3, las=0)
         box()
         
         # 1D Occ
         par(mar=c(5,5,2,2))
         AvgOcc = colSums(Occ_matrix)
         plot(-beforeRef:afterRef, AvgOcc, axes=FALSE, typ='l', ann=FALSE, col = "blue",
              xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*max(AvgOcc)),
              panel.first = c(abline(v = seq(-500, 500, 500), col = "lightgray", lty = "dashed")))
         box()
         # x axis
         par(tcl= -0.2)
         axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(1, at=seq(-1000, 1000, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
         # y axis
         axis(2, at=seq(0,5,0.2), cex.axis=1.25)
         title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Average occupancy", cex.lab=1.4)
         
         # 2D Occ
         par(mar=c(5,5,2,2))
         image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
         # x axis
         par(tcl= -0.2)
         axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
         # y axis
         par(tcl= -0.2)
         axis(4, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(4, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
         
         abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
         box()
         
         # Empty plot
         plot(0, type='n', axes=FALSE, ann=FALSE)
         
         # Length histogram
         par(mar=c(5,5,2,1))
         plot(LengthHist, Lmin:Lmax, axes=FALSE, typ='l', ann=FALSE, col = "blue",
              ylim=c(Lmin, Lmax), xaxs="i", yaxs="i", xlim = c(0, 1.05*max(LengthHist)),
              panel.first = c(abline(h = c(100, 150), col = "lightgray", lty = "dashed")))
         box()
         # c axis
         axis(1, at=seq(0,1.05*max(LengthHist),1), cex.axis=1.25)
         # y axis
         par(tcl= -0.2)
         axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
         title(xlab="Percentage (%)", ylab="Fragment length (bp)", cex.lab=1.4)
         
         garbage = dev.off()
       },
       TTS={
         Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors
         
         pdf(paste("2D_occupancy/Occ_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
         layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
         
         # Empty plot
         plot(0, type='n', axes=FALSE, ann=FALSE)
         
         # Add scale
         par(mar=c(5,6,2,0.5))
         image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
         axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
         mtext("Relative coverage (%)", cex=1, side=2, padj=-3, las=0)
         box()
         
         # 1D Occ
         par(mar=c(5,5,2,2))
         AvgOcc = colSums(Occ_matrix)
         plot(-beforeRef:afterRef, AvgOcc, axes=FALSE, typ='l', ann=FALSE, col = "blue",
              xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*max(AvgOcc)),
              panel.first = c(abline(v = seq(-500, 500, 500), col = "lightgray", lty = "dashed")))
         box()
         # x axis
         par(tcl= -0.2)
         axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(1, at=seq(-1000, 1000, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
         # y axis
         axis(2, at=seq(0,5,0.2), cex.axis=1.25)
         title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Average occupancy", cex.lab=1.4)
         
         # 2D Occ
         par(mar=c(5,5,2,2))
         image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
         # x axis
         par(tcl= -0.2)
         axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
         # y axis
         par(tcl= -0.2)
         axis(4, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(4, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
         
         abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
         box()
         
         # Empty plot
         plot(0, type='n', axes=FALSE, ann=FALSE)
         
         # Length histogram
         par(mar=c(5,5,2,1))
         plot(LengthHist, Lmin:Lmax, axes=FALSE, typ='l', ann=FALSE, col = "blue",
              ylim=c(Lmin, Lmax), xaxs="i", yaxs="i", xlim = c(0, 1.05*max(LengthHist)),
              panel.first = c(abline(h = c(100, 150), col = "lightgray", lty = "dashed")))
         box()
         # c axis
         axis(1, at=seq(0,1.05*max(LengthHist),1), cex.axis=1.25)
         # y axis
         par(tcl= -0.2)
         axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
         title(xlab="Percentage (%)", ylab="Fragment length (bp)", cex.lab=1.4)
         
         garbage = dev.off()
       },
       Plus1={
         Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors, which should be really zeros
         
         pdf(paste("2D_occupancy/Occ_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
         layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
         
         # Empty plot
         plot(0, type='n', axes=FALSE, ann=FALSE)
         
         # Add scale
         par(mar=c(5,6,2,0.5))
         image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
         axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
         mtext("Relative coverage (%)", cex=1, side=2, padj=-3, las=0)
         box()
         
         # 1D Occ
         par(mar=c(5,5,2,2))
         AvgOcc = colSums(Occ_matrix)
         plot(-beforeRef:afterRef, AvgOcc, axes=FALSE, typ='l', ann=FALSE, col = "blue",
              xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*max(AvgOcc)),
              panel.first = c(abline(v = seq(-500, 500, 500), col = "lightgray", lty = "dashed")))
         box()
         # x axis
         par(tcl= -0.2)
         axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(1, at=seq(-1000, 1000, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
         # y axis
         axis(2, at=seq(0,5,0.2), cex.axis=1.25)
         title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Average occupancy", cex.lab=1.4)
         
         # 2D Occ
         par(mar=c(5,5,2,2))
         image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
         # x axis
         par(tcl= -0.2)
         axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
         # y axis
         par(tcl= -0.2)
         axis(4, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(4, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
         
         abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
         box()
         
         # Empty plot
         plot(0, type='n', axes=FALSE, ann=FALSE)
         
         # Length histogram
         par(mar=c(5,5,2,1))
         plot(LengthHist, Lmin:Lmax, axes=FALSE, typ='l', ann=FALSE, col = "blue",
              ylim=c(Lmin, Lmax), xaxs="i", yaxs="i", xlim = c(0, 1.05*max(LengthHist)),
              panel.first = c(abline(h = c(100, 150), col = "lightgray", lty = "dashed")))
         box()
         # c axis
         axis(1, at=seq(0,1.05*max(LengthHist),1), cex.axis=1.25)
         # y axis
         par(tcl= -0.2)
         axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
         par(tcl= -0.5)
         axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
         title(xlab="Percentage (%)", ylab="Fragment length (bp)", cex.lab=1.4)
         
         garbage = dev.off()
       }
)
