#!/usr/bin/env Rscript

library("optparse")

options = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Dataset file name [options: BED or BAM format]"),
  make_option(c("-t", "--type"), type="character", default="occ", 
              help="Types of distribution to plot [options: occ, dyads, fivePrime_ends, threePrime_ends; default = %default]"),
  make_option(c("-o", "--organism"), type="character", default="sacCer3", 
              help="Genome version [options: sacCer3 (default), dm3, dm6, ce10, ce11, mm9, mm10, hg18, hg19]"),
  make_option(c("-r", "--reference"), type="character", default="TSS", 
              help="Reference points to align [options: TSS (default), TTS, Plus1]"),
  make_option(c("-s", "--sites"), type="character", default=NULL, 
              help="Reference points in BED format"),
  make_option(c("-a", "--align"), type="character", default="center", 
              help="What points of the provided intervals to align? [options: center (default), fivePrime, threePrime]"),
  make_option(c("--siteLabel"), type="character", default="Sites", 
              help="Label for the aligned sites [default = %default]"),
  make_option(c("-l", "--minLength"), type="integer", default=50, 
              help="The smallest DNA fragment to be considered [default = %default]"),
  make_option(c("-L", "--maxLength"), type="integer", default=200, 
              help="The largest DNA fragment to be considered [default = %default]"),
  make_option(c("-u", "--upstream"), type="integer", default=1000, 
              help="Length of the upstream region that will be plotted [default = %default]"),
  make_option(c("-d", "--downstream"), type="integer", default=1000, 
              help="Length of the downstream region that will be plotted [default = %default]"),
  make_option(c("-m", "--colorScaleMax"), type="double", default=NULL, 
              help="Maximum value on the color scale (e.g. 0.02)"),
  make_option(c("--simplifyPlot"), type="character", default="off", 
              help="Simplify the plot (show only the 2D Occ.) [options: on, off (default)]"),
  make_option(c("--squeezePlot"), type="character", default="off", 
              help="Simplify the plot (show only the 2D Occ.) and squeeze the heat map [options: on, off (default)]")
)

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least the dataset file name must be supplied.", call.=FALSE)
}

if (opt$squeezePlot == "on"){
  opt$simplifyPlot = "on"
}

##################
# Initialization #
##################
# Type of plot
plot.type = toupper(opt$type)
plot.type = strsplit(plot.type, ',')[[1]]

# Genome
genome = opt$organism

# Data file name
inputFilename = opt$file

# Reference label
siteLabel = opt$siteLabel

if (! is.null(opt$colorScaleMax)) {
  colorScaleMax = opt$colorScaleMax
}

switch(genome, 
       sacCer3={
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
         
         # Load sacCer3 annotations
         annotations = read.csv("Annotations/sacCer3_annotations.csv", header=TRUE, stringsAsFactors=FALSE)
       },
       dm3={
         chrLen = c(23011544,21146708,24543557,27905053,1351857,22422827)
         names(chrLen) = c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
         
         # Get fly annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host='dec2014.archive.ensembl.org') # use this host for getting the annotations for dm3=BDGP5 (Fruitfly release 78)
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get BDGP5 annotations
         ensembl = useMart(host='dec2014.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='dmelanogaster_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       dm6={
         chrLen = c(23513712,25286936,28110227,32079331,1348131,23542271)
         names(chrLen) = c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
         
         # Get fly annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host='may2017.archive.ensembl.org') # use this host for getting the annotations for dm6=BDGP6
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get BDGP6 annotations
         ensembl = useMart(host='may2017.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='dmelanogaster_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       ce10={
         chrLen = c(15072423,15279345,13783700,17493793,20924149,17718866)
         names(chrLen) = c("chrI","chrII","chrIII","chrIV","chrV","chrX")
         
         # Get fly annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host=may2012.archive.ensembl.org') # use this host for getting the annotations for ce10=WS220
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get ce10 annotations
         ensembl = useMart(host='may2012.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='celegans_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       ce11={
         chrLen = c(15072434,15279421,13783801,17493829,20924180,17718942)
         names(chrLen) = c("chrI","chrII","chrIII","chrIV","chrV","chrX")
         
         # Get fly annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host=may2017.archive.ensembl.org') # use this host for getting the annotations for ce11=WBcel235
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get ce11 annotations
         ensembl = useMart(host='may2017.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='celegans_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       mm9={
         chrLen = c(197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555)
         names(chrLen) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
         
         # Get mouse annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host='may2012.archive.ensembl.org') # use this host for getting the annotations for mm9, hg19
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get mm9 annotations
         ensembl = useMart(host='may2012.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='mmusculus_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       mm10={
         chrLen = c(195471971,182113224,160039680,156508116,151834684,149736546,145441459,129401213,124595110,130694993,122082543,120129022,120421639,124902244,104043685,98207768,94987271,90702639,61431566,171031299,91744698)
         names(chrLen) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
         
         # Get mouse annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # Get mm9 annotations
         ensembl = useMart(biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='mmusculus_gene_ensembl')
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       hg18={
         chrLen = c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954)
         names(chrLen) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
         
         # Get human annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host='may2009.archive.ensembl.org') # use this host for getting the annotations for hg18
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get hg18 annotations
         ensembl = useMart(host='may2009.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='hsapiens_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       hg19={
         chrLen = c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
         names(chrLen) = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
         
         # Get human annotations
         suppressPackageStartupMessages(library(biomaRt))
         
         # listMarts(host='feb2014.archive.ensembl.org') # use this host for getting the annotations for hg19
         # other Ensembl archives: http://www.ensembl.org/info/website/archives/index.html
         
         # Get hg19 annotations
         ensembl = useMart(host='feb2014.archive.ensembl.org', 
                           biomart='ENSEMBL_MART_ENSEMBL', 
                           dataset='hsapiens_gene_ensembl')
         # listFilters(mart=ensembl)
         
         annot = getBM(c("refseq_mrna", "chromosome_name", "strand", "transcript_start", "transcript_end"), mart=ensembl)
         annot = annot[! (annot$refseq_mrna == ''), ]
         
         strand = annot$strand
         txStart = annot$transcript_start
         txEnd = annot$transcript_end
         
         TSS = txStart
         TTS = txEnd
         TSS[strand == -1] = txEnd[strand == -1]
         TTS[strand == -1] = txStart[strand == -1]
         
         annotations = data.frame(Transcript = annot$refseq_mrna, 
                                  Chr = paste("chr", annot$chromosome_name, sep=""), 
                                  Strand = strand,
                                  TSS = TSS,
                                  TTS = TTS)
       },
       {
         print_help(opt_parser)
         stop(paste("Genome ", genome, " is not supported yet. Supported options: sacCer3, dm3, dm6, ce10, ce11, mm9, mm10, hg18, hg19.", sep=""), call.=FALSE)
       }
)

# Type of alignments, i.e. reference points
selectedReference = opt$reference

# Size selection parameters: specify the interval of lengths to be analyzed
Lmin = opt$minLength   # the smallest DNA fragment to be considered
Lmax = opt$maxLength  # the largest DNA fragment to be considered

# Window selection parameters
beforeRef = opt$upstream  # length of the upstream region that will be plotted
afterRef = opt$downstream   # length of the downstream region that will be plotted

# Load the necessary R packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
})

#########################################
# Import the paired-end sequencing data #
#########################################
library(tools)
inputType = toupper(file_ext(inputFilename))
switch(inputType, 
       BED={
         sample.name = sub(".bed", "", inputFilename)
         reads = import.bed(inputFilename)
       },
       BAM={
         sample.name = sub(".bam", "", inputFilename)
         suppressPackageStartupMessages(library(Rsamtools))
         all_fields = c("rname", "pos", "isize")
         param = ScanBamParam(what = all_fields, 
                              flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, 
                                                 isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                                                 isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                                                 isNotPassingQualityControls = FALSE))
         bam = scanBam(inputFilename, param=param)
         
         # Keep only the proper reads, with the length > 0
         posStrandReads = (bam[[1]]$isize > 0)
         
         reads = GRanges(seqnames=Rle(bam[[1]]$rname[posStrandReads]),
                         ranges = IRanges(start=bam[[1]]$pos[posStrandReads], width=bam[[1]]$isize[posStrandReads]),
                         strand = "*")
         rm(bam)
       },
       {
         print_help(opt_parser)
         stop("File type not supported! Please provide a BED or BAM file as input.", call.=FALSE)
       }
)

# Remove problematic regions (e.g. repeated rDNA in the case of yeast)
if (genome == "sacCer3"){
  # Discard the reads from the rDNA region, chrXII:451000-469000
  rDNAregion = GRanges(seqnames = "chrXII",
                       ranges = IRanges(start=451000, end=469000),
                       strand = "*")
  rDNAInd = overlapsAny(reads, rDNAregion, ignore.strand=TRUE)
  reads = reads[!rDNAInd]
  
  # Discard the reads from the yeast mitochondrial DNA
  reads = reads[seqnames(reads) != 'chrM']
}


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
noChr = length(chrLen)
extraReads = GRanges(seqnames = names(chrLen),
                     ranges = IRanges(start = chrLen, width = 1),
                     strand = rep("+",noChr),
                     seqlengths = chrLen)
reads = c(extraReads, reads)


#########################
# Windows to be aligned #
#########################
if (is.null(opt$sites)) {
  switch(selectedReference, 
         TSS={
           ReferencePos = annotations$TSS
           ReferenceChr = annotations$Chr
           RefStrand = annotations$Strand
           
           validValues = (! is.nan(ReferencePos))
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
           wholeChr = GRanges(seqnames=names(chrLen),
                              IRanges(start=rep(1, length(chrLen)),
                                      end=chrLen))
           # ... and remove the windows that fall outside of the chromosome edges
           ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                               type="within", ignore.strand=TRUE)
         },
         TTS={
           ReferencePos = annotations$TTS
           ReferenceChr = annotations$Chr
           RefStrand = annotations$Strand
           
           validValues = (! is.nan(ReferencePos))
           # if (genome == "sacCer3"){    # in our publication we have used an additional filter to consider only the TTSs that are away from other TSSs, between convergent genes
           #                              # Uncomment this code block to use the same filter
           #   # Keep only TTSs that are away from other TSSs
           #   validValues = (! is.nan(ReferencePos)) & ((annotations$NDR3pLength > 1000) | (annotations$NDR3pQualifier == "convergent"))
           # }
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
           wholeChr = GRanges(seqnames=names(chrLen),
                              IRanges(start=rep(1, length(chrLen)),
                                      end=chrLen))
           # ... and remove the windows that fall outside of the chromosome edges
           ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                               type="within", ignore.strand=TRUE)  
         },
         Plus1={
           
           if (genome == 'sacCer3'){
             ReferencePos = annotations$Plus1
             ReferenceChr = annotations$Chr
             RefStrand = annotations$Strand
             
             validValues = (! is.nan(ReferencePos))
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
             wholeChr = GRanges(seqnames=names(chrLen),
                                IRanges(start=rep(1, length(chrLen)),
                                        end=chrLen))
             # ... and remove the windows that fall outside of the chromosome edges
             ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                                 type="within", ignore.strand=TRUE)
           } else {
             print_help(opt_parser)
             stop("Plus1 annotations are available only for the sacCer3 genome. For other genomes, the available options for --reference are: TSS (default), TTS", call.=FALSE)
           }
         }
  ) 
} else {
  referenceFilename = opt$sites
  sites = import.bed(referenceFilename)
  
  switch(opt$align, 
         "center"={
           ReferenceGRanges = resize(sites, width = 1, fix="center")
         },
         "fivePrime"={
           ReferenceGRanges = resize(sites, width = 1, fix="start")
         },
         "threePrime"={
           ReferenceGRanges = resize(sites, width = 1, fix="end")
         })
  
  Watson = as.vector(strand(ReferenceGRanges) != '-')      # ('+' or '*')
  ReferencePos = start(ReferenceGRanges)
  leftEdge = ReferencePos-beforeRef
  rightEdge = ReferencePos+afterRef
  
  leftEdgeCrick = ReferencePos-afterRef
  rightEdgeCrick = ReferencePos+beforeRef
  
  leftEdge[!Watson] = leftEdgeCrick[!Watson]
  rightEdge[!Watson] = rightEdgeCrick[!Watson]
  
  ReferenceGRanges = GRanges(seqnames=seqnames(ReferenceGRanges),
                             IRanges(start=leftEdge,
                                     end=rightEdge),
                             strand=strand(ReferenceGRanges))
  
  
  # Construct GRanges for the entire chromosomes...
  wholeChr = GRanges(seqnames=names(chrLen),
                     IRanges(start=rep(1, length(chrLen)),
                             end=chrLen))
  # ... and remove the windows that fall outside of the chromosome edges
  ReferenceGRanges = subsetByOverlaps(ReferenceGRanges, wholeChr,
                                      type="within", ignore.strand=TRUE)
}

# Create folder where figures are saved
if (is.null(opt$sites)) {
  switch(selectedReference, 
         TSS={
           dir.create("2D_occupancy_TSS", showWarnings = FALSE, recursive = TRUE)
         },
         TTS={
           dir.create("2D_occupancy_TTS", showWarnings = FALSE, recursive = TRUE)
         },
         Plus1={
           dir.create("2D_occupancy_Plus1", showWarnings = FALSE, recursive = TRUE)
         }
  )
} else {
  dir.create(paste("2D_occupancy_", opt$siteLabel, sep=""), showWarnings = FALSE, recursive = TRUE)
}

# Some functions used below

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





if ("OCC" %in% plot.type){
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
  ## Initialize the 2D occ. matrix
  Occ_matrix = matrix(data = 0, nrow = Lmax-Lmin+1, ncol = 1+beforeRef+afterRef)
  
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
  if (is.null(opt$sites)) {
    switch(selectedReference, 
           TSS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TSS/Occ_matrix.TSS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           TTS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TTS/Occ_matrix.TTS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           Plus1={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_Plus1/Occ_matrix.Plus1.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           }
    )
  } else {
    save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, siteLabel, file=paste("2D_occupancy_", opt$siteLabel, "/Occ_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
  }
  
  ################################
  # Plot the 2D Occupancy matrix #
  ################################
  # Plot the figure
  Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors
  if (! is.null(opt$colorScaleMax)) {
    Occ_matrix[Occ_matrix >= opt$colorScaleMax] = opt$colorScaleMax - 1e-10   # set a maximum threshold for the 2D Occ matrix
  }
  
  if (opt$squeezePlot == "on") { #Plot simplified & squeezed plot
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/Occ_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             TTS={pdf(paste("2D_occupancy_TTS/Occ_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             Plus1={pdf(paste("2D_occupancy_Plus1/Occ_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/Occ_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)
    }
  } else {
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/Occ_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             TTS={pdf(paste("2D_occupancy_TTS/Occ_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             Plus1={pdf(paste("2D_occupancy_Plus1/Occ_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/Occ_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
    }
  }
  
  if (opt$simplifyPlot == "on") { #Plot simplified plot
    
    if (opt$squeezePlot == "on") {
      layout(matrix(c(1,2), ncol=2), widths=c(4,1))
    } else {
      layout(matrix(c(1,2), ncol=2), widths=c(5,1))
    }
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
    # y axis
    par(tcl= -0.2)
    axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
    par(tcl= -0.5)
    axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    
    abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Fragment length (bp)", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)})}
    box()
    
    # Add scale
    if (opt$squeezePlot == "on") {
      par(mar=c(25,2,2,1))
    } else {
      par(mar=c(15,3,2,1))
    }
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    if (opt$squeezePlot == "off") {
      mtext("Relative coverage (%)", cex=1.4, side=2, padj=-3.5, las=0)
    }
    box()
    
    garbage = dev.off()
    
  } else { #Plot multi-panel figure
    layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
    
    # Empty plot
    plot(0, type='n', axes=FALSE, ann=FALSE)
    
    # Add scale
    par(mar=c(15,6,2,0))
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    mtext("Relative coverage (%)", cex.lab=1.4, side=2, padj=-3.5, las=0)
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
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Average occupancy", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Average occupancy", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Average occupancy", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Average occupancy", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Average occupancy", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Average occupancy", cex.lab=1.4)})}
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
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
}





if ("DYADS" %in% plot.type){
  resized_reads = resize(reads, width = 1, fix = "center")
  
  #################
  # Normalization #
  #################
  # Compute rescaling coefficients, such that after rescaling the average dyad density is 1, for each chromosome 
  Occ = coverage(resized_reads)
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
  ## Initialize the 2D occ. matrix
  Occ_matrix = matrix(data = 0, nrow = Lmax-Lmin+1, ncol = 1+beforeRef+afterRef)
  
  # For each fragment size, compute the corresponding coverage/occupancy
  for(L in Lmin:Lmax) {
    
    # Keep only the reads with the specific length L
    goodReadsInd = (readLength == L)
    goodReadsInd[1:noChr] = TRUE # extra 1bp reads at the ends all chromosomes
    goodReads = resized_reads[goodReadsInd]
    
    if (length(goodReads) > noChr) {
      # Compute average occupancy
      Occ = coverage(goodReads, weight=coverageWeight)
      Occ_matrix[L-Lmin+1,] = colMeans(AlignRegions(Occ, ReferenceGRanges))
    }
  }
  
  # Save the 2D Occupancy matrix
  if (is.null(opt$sites)) {
    switch(selectedReference, 
           TSS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TSS/Dyads_matrix.TSS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           TTS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TTS/Dyads_matrix.TTS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           Plus1={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_Plus1/Dyads_matrix.Plus1.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           }
    )
  } else {
    save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, siteLabel, file=paste("2D_occupancy_", opt$siteLabel, "/Dyads_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
  }
  
  
  ################################
  # Plot the 2D Occupancy matrix #
  ################################
  # Plot the figure
  Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors
  if (! is.null(opt$colorScaleMax)) {
    Occ_matrix[Occ_matrix >= opt$colorScaleMax] = opt$colorScaleMax - 1e-10   # set a maximum threshold for the 2D Occ matrix
  }
  
  if (opt$squeezePlot == "on") { #Plot simplified & squeezed plot
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/Dyads_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             TTS={pdf(paste("2D_occupancy_TTS/Dyads_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             Plus1={pdf(paste("2D_occupancy_Plus1/Dyads_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/Dyads_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)
    }
  } else {
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/Dyads_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             TTS={pdf(paste("2D_occupancy_TTS/Dyads_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             Plus1={pdf(paste("2D_occupancy_Plus1/Dyads_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/Dyads_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
    }
  }
  
  if (opt$simplifyPlot == "on") { #Plot simplified plot
    
    if (opt$squeezePlot == "on") {
      layout(matrix(c(1,2), ncol=2), widths=c(4,1))
    } else {
      layout(matrix(c(1,2), ncol=2), widths=c(5,1))
    }
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
    # y axis
    par(tcl= -0.2)
    axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
    par(tcl= -0.5)
    axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    
    abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Fragment length (bp)", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)})}
    box()
    
    # Add scale
    if (opt$squeezePlot == "on") {
      par(mar=c(25,2,2,1))
    } else {
      par(mar=c(15,3,2,1))
    }
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    if (opt$squeezePlot == "off") {
      mtext("Relative coverage (%)", cex=1.4, side=2, padj=-3.5, las=0)
    }
    box()
    
    garbage = dev.off()
    
  } else { #Plot multi-panel figure
    layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
    
    # Empty plot
    plot(0, type='n', axes=FALSE, ann=FALSE)
    
    # Add scale
    par(mar=c(15,6,2,0))
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    mtext("Relative coverage (%)", cex.lab=1.4, side=2, padj=-3.5, las=0)
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
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Average dyad density", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Average dyad density", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Average dyad density", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Average dyad density", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Average dyad density", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Average dyad density", cex.lab=1.4)})}
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
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
}





if ("FIVEPRIME_ENDS" %in% plot.type){
  resized_reads = resize(reads, width = 1, fix = "start")
  
  #################
  # Normalization #
  #################
  # Compute rescaling coefficients, such that after rescaling the average dyad density is 1, for each chromosome 
  Occ = coverage(resized_reads)
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
  ## Initialize the 2D occ. matrix
  Occ_matrix = matrix(data = 0, nrow = Lmax-Lmin+1, ncol = 1+beforeRef+afterRef)

  # For each fragment size, compute the corresponding coverage/occupancy
  for(L in Lmin:Lmax) {
    
    # Keep only the reads with the specific length L
    goodReadsInd = (readLength == L)
    goodReadsInd[1:noChr] = TRUE # extra 1bp reads at the ends all chromosomes
    goodReads = resized_reads[goodReadsInd]
    
    if (length(goodReads) > noChr) {
      # Compute average occupancy
      Occ = coverage(goodReads, weight=coverageWeight)
      Occ_matrix[L-Lmin+1,] = colMeans(AlignRegions(Occ, ReferenceGRanges))
    }
  }
  
  # Save the 2D Occupancy matrix
  if (is.null(opt$sites)) {
    switch(selectedReference, 
           TSS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TSS/fivePrime_ends_matrix.TSS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           TTS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TTS/fivePrime_ends_matrix.TTS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           Plus1={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_Plus1/fivePrime_ends_matrix.Plus1.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           }
    )
  } else {
    save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, siteLabel, file=paste("2D_occupancy_", opt$siteLabel, "/fivePrime_ends_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
  }
  
  
  ################################
  # Plot the 2D Occupancy matrix #
  ################################
  # Plot the figure
  Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors
  if (! is.null(opt$colorScaleMax)) {
    Occ_matrix[Occ_matrix >= opt$colorScaleMax] = opt$colorScaleMax - 1e-10   # set a maximum threshold for the 2D Occ matrix
  }
  
  if (opt$squeezePlot == "on") { #Plot simplified & squeezed plot
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/fivePrime_ends_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             TTS={pdf(paste("2D_occupancy_TTS/fivePrime_ends_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             Plus1={pdf(paste("2D_occupancy_Plus1/fivePrime_ends_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/fivePrime_ends_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)
    }
  } else {
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/fivePrime_ends_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             TTS={pdf(paste("2D_occupancy_TTS/fivePrime_ends_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             Plus1={pdf(paste("2D_occupancy_Plus1/fivePrime_ends_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/fivePrime_ends_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
    }
  }
  
  if (opt$simplifyPlot == "on") { #Plot simplified plot
    
    if (opt$squeezePlot == "on") {
      layout(matrix(c(1,2), ncol=2), widths=c(4,1))
    } else {
      layout(matrix(c(1,2), ncol=2), widths=c(5,1))
    }
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
    # y axis
    par(tcl= -0.2)
    axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
    par(tcl= -0.5)
    axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    
    abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Fragment length (bp)", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)})}
    box()
    
    # Add scale
    if (opt$squeezePlot == "on") {
      par(mar=c(25,2,2,1))
    } else {
      par(mar=c(15,3,2,1))
    }
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    if (opt$squeezePlot == "off") {
      mtext("Relative coverage (%)", cex=1.4, side=2, padj=-3.5, las=0)
    }
    box()
    
    garbage = dev.off()
    
  } else { #Plot multi-panel figure
    layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
    
    # Empty plot
    plot(0, type='n', axes=FALSE, ann=FALSE)
    
    # Add scale
    par(mar=c(15,6,2,0))
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    mtext("Relative coverage (%)", cex.lab=1.4, side=2, padj=-3.5, las=0)
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
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Avg. 5' end density", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Avg. 5' end density", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Avg. 5' end density", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Avg. 5' end density", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Avg. 5' end density", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Avg. 5' end density", cex.lab=1.4)})}
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
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
}





if ("THREEPRIME_ENDS" %in% plot.type){
  resized_reads = resize(reads, width = 1, fix = "end")
  
  #################
  # Normalization #
  #################
  # Compute rescaling coefficients, such that after rescaling the average dyad density is 1, for each chromosome 
  Occ = coverage(resized_reads)
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
  ## Initialize the 2D occ. matrix
  Occ_matrix = matrix(data = 0, nrow = Lmax-Lmin+1, ncol = 1+beforeRef+afterRef)

  # For each fragment size, compute the corresponding coverage/occupancy
  for(L in Lmin:Lmax) {
    
    # Keep only the reads with the specific length L
    goodReadsInd = (readLength == L)
    goodReadsInd[1:noChr] = TRUE # extra 1bp reads at the ends all chromosomes
    goodReads = resized_reads[goodReadsInd]
    
    if (length(goodReads) > noChr) {
      # Compute average occupancy
      Occ = coverage(goodReads, weight=coverageWeight)
      Occ_matrix[L-Lmin+1,] = colMeans(AlignRegions(Occ, ReferenceGRanges))
    }
  }
  
  # Save the 2D Occupancy matrix
  if (is.null(opt$sites)) {
    switch(selectedReference, 
           TSS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TSS/threePrime_ends_matrix.TSS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           TTS={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_TTS/threePrime_ends_matrix.TTS.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           },
           Plus1={
             save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, file=paste("2D_occupancy_Plus1/threePrime_ends_matrix.Plus1.", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
           }
    )
  } else {
    save(Occ_matrix, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, LengthHist, sample.name, siteLabel, file=paste("2D_occupancy_", opt$siteLabel, "/threePrime_ends_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".RData", sep=""))
  }
  
  
  ################################
  # Plot the 2D Occupancy matrix #
  ################################
  # Plot the figure
  Occ_matrix[Occ_matrix < 0] = 0 # eliminate rounding errors
  if (! is.null(opt$colorScaleMax)) {
    Occ_matrix[Occ_matrix >= opt$colorScaleMax] = opt$colorScaleMax - 1e-10   # set a maximum threshold for the 2D Occ matrix
  }
  
  if (opt$squeezePlot == "on") { #Plot simplified & squeezed plot
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/threePrime_ends_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             TTS={pdf(paste("2D_occupancy_TTS/threePrime_ends_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)},
             Plus1={pdf(paste("2D_occupancy_Plus1/threePrime_ends_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/threePrime_ends_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=5.5, h=8)
    }
  } else {
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={pdf(paste("2D_occupancy_TSS/threePrime_ends_matrix_TSS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             TTS={pdf(paste("2D_occupancy_TTS/threePrime_ends_matrix_TTS.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)},
             Plus1={pdf(paste("2D_occupancy_Plus1/threePrime_ends_matrix_Plus1.", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)})
    } else {  # Specific sites were provided in a separate BED file
      pdf(paste("2D_occupancy_", opt$siteLabel, "/threePrime_ends_matrix.", siteLabel, ".", Lmin, "_", Lmax, ".", sample.name, ".pdf", sep=""), w=7, h=6)
    }
  }

  if (opt$simplifyPlot == "on") { #Plot simplified plot
    
    if (opt$squeezePlot == "on") {
      layout(matrix(c(1,2), ncol=2), widths=c(4,1))
    } else {
      layout(matrix(c(1,2), ncol=2), widths=c(5,1))
    }
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(1, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(1, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
    # y axis
    par(tcl= -0.2)
    axis(2, at=seq(Lmin, Lmax, by=10), labels=F, lwd=1, lwd.ticks=1)
    par(tcl= -0.5)
    axis(2, at=seq(Lmin, Lmax, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    
    abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 1)
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Fragment length (bp)", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Fragment length (bp)", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Fragment length (bp)", cex.lab=1.4)})}
    box()
    
    # Add scale
    if (opt$squeezePlot == "on") {
      par(mar=c(25,2,2,1))
    } else {
      par(mar=c(15,3,2,1))
    }
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    if (opt$squeezePlot == "off") {
      mtext("Relative coverage (%)", cex=1.4, side=2, padj=-3.5, las=0)
    }
    box()
    
    garbage = dev.off()
    
  } else { #Plot multi-panel figure
    layout(matrix(c(1,2,3,4,5,6), ncol=3), widths=c(1,4,2), heights=c(2,4))
    
    # Empty plot
    plot(0, type='n', axes=FALSE, ann=FALSE)
    
    # Add scale
    par(mar=c(15,6,2,0))
    if (is.null(opt$colorScaleMax)) {
      image.scale(col = matlab.like(100), breaks = seq(0, max(Occ_matrix+1e-10), max(Occ_matrix+1e-10)/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(max(Occ_matrix) / 0.005) > 10) {
        if (round(max(Occ_matrix) / 0.01) > 10) {
          axis(2, at=seq(0, max(Occ_matrix), 0.05), labels=100*seq(0, max(Occ_matrix), 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, max(Occ_matrix), 0.01), labels=100*seq(0, max(Occ_matrix), 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, max(Occ_matrix), 0.005), labels=100*seq(0, max(Occ_matrix), 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    } else {
      image.scale(col = matlab.like(100), breaks = seq(0, colorScaleMax, colorScaleMax/100), horiz=FALSE, xlab="", ylab="", yaxt="n")
      if (round(colorScaleMax / 0.005) > 10) {
        if (round(colorScaleMax / 0.01) > 10) {
          axis(2, at=seq(0, colorScaleMax, 0.05), labels=100*seq(0, colorScaleMax, 0.05), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        } else {
          axis(2, at=seq(0, colorScaleMax, 0.01), labels=100*seq(0, colorScaleMax, 0.01), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
        }
      } else {
        axis(2, at=seq(0, colorScaleMax, 0.005), labels=100*seq(0, colorScaleMax, 0.005), cex.axis=1.25, las=2) # las: labels are parallel (=0) or perpendicular(=2) to axis
      }
    }
    mtext("Relative coverage (%)", cex.lab=1.4, side=2, padj=-3.5, las=0)
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
    
    if (is.null(opt$sites)) {
      switch(selectedReference, 
             TSS={title(main=sample.name, font.main = 1, xlab="Position relative to TSS (bp)", ylab="Avg. 3' end density", cex.lab=1.4)},
             TTS={title(main=sample.name, font.main = 1, xlab="Position relative to TTS (bp)", ylab="Avg. 3' end density", cex.lab=1.4)},
             Plus1={title(main=sample.name, font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Avg. 3' end density", cex.lab=1.4)})
    } else {  # Specific sites were provided in a separate BED file
      switch(opt$align, 
             "center"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to center of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Avg. 3' end density", cex.lab=1.4)},
             "fivePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 5' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Avg. 3' end density", cex.lab=1.4)},
             "threePrime"={title(main=sample.name, font.main = 1, 
                             xlab=paste("Position relative to 3' end of ", gsub("_", " ", opt$siteLabel)," (bp)", sep=""), 
                             ylab="Avg. 3' end density", cex.lab=1.4)})}
    
    # 2D Occ
    par(mar=c(5,5,2,2))
    if (is.null(opt$colorScaleMax)) {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, max(Occ_matrix), max(Occ_matrix)/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    } else {
      image(-beforeRef:afterRef, Lmin:Lmax, t(Occ_matrix), col=matlab.like(100), ylim=c(Lmin, Lmax), breaks = seq(0, colorScaleMax, colorScaleMax/100), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
    }
    # x axis
    if (beforeRef >= 500 | afterRef >= 500) {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=100), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=500), lwd=0, lwd.ticks=1, cex.axis=1.25)
    } else {
      par(tcl= -0.2)
      axis(3, at=seq(-1000, 1000, by=10), labels=F, lwd=1, lwd.ticks=1)
      par(tcl= -0.5)
      axis(3, at=seq(-1000, 1000, by=50), lwd=0, lwd.ticks=1, cex.axis=1.25)
    }
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
}