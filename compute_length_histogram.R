#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="Data file names, separated by commas only (BAM format)"),
  make_option(c("-l", "--minLength"), type="integer", default=0, 
              help="The smallest DNA fragment to be considered [default = %default]"),
  make_option(c("-L", "--maxLength"), type="integer", default=500, 
              help="The largest DNA fragment to be considered [default = %default]"),
  make_option(c("-s", "--statistics"), type="character", default="on", 
              help="Include statistics in the plot [options: on, off; default = %default]"),
  make_option(c("-o", "--outputs"), type="character", default="pdf,csv,RData", 
              help="Types of outputs to be generated, separated by commas only [options: pdf, csv, RData; default = %default]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)


if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least the dataset file name must be supplied.", call.=FALSE)
}

minLength = opt$minLength
maxLength = opt$maxLength

##################
# Initialization #
##################
# Load the necessary R packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
  library(Rsamtools)
  library(ggplot2)
  library(gridExtra)
})

# Data files
bam.files = opt$files
bam.files = strsplit(bam.files, ',')[[1]]
noFiles = length(bam.files)

for (f in 1:noFiles){
  #########################################
  # Import the paired-end sequencing data #
  #########################################
  # Data file name
  inputFilename = bam.files[f]
  sample.name = sub(".bam", "", inputFilename)
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
  readLength = width(reads)
  TotalNoReads = length(reads)
  
  #########################################
  # Compute the fragment length histogram #
  #########################################
  output.files = toupper(opt$outputs)
  filesToGenerate = strsplit(output.files, ',')[[1]]
  
  # Compute the histogram
  h = hist(readLength, breaks=seq(from = 0.5, to = 1000.5, by = 1), plot=FALSE)
  
  # Create folder
  dir.create("Length_histograms", showWarnings = FALSE)
  
  if ("PDF" %in% filesToGenerate){
    # Plot the histogram using ggplot2
    
    df = data.frame(x = h$mids, y = 100*h$density)
    myplot = ggplot(df, aes(x = x, y = y)) + geom_line(colour="#56B4E9") +
      scale_x_continuous(expand = c(0, 0), limits = c(minLength, maxLength)) +
      scale_y_continuous(expand = c(0, 0)) + theme_bw() + 
      theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("Fragment length (bp)") + 
      ylab("Percentage (%)") +
      ggtitle(sample.name) +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=18,face="bold"),
            plot.title=element_text(size=18,face="bold"))
    
    if (opt$statistics == 'on') {
      # Add table of quantiles
      mytable = data.frame(Percentile = c("5%", "10%", "25%", "50%", "75%", "90%", "95%"), 
                           Length = quantile(readLength, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
      
      myplot = myplot + annotation_custom(tableGrob(mytable, rows=NULL), xmin=minLength + 0.65*(maxLength-minLength), xmax=minLength + 0.95*(maxLength-minLength), ymin=0.27*max(100*h$density), ymax=0.95*max(100*h$density))
      
      suppressWarnings(
      ggsave(filename=paste("Length_histograms/Length_histogram.", sample.name, ".w_stats.pdf", sep=""), 
             plot=myplot,
             width = 5, height = 4, units = "in"))
    } else {
      suppressWarnings(
      ggsave(filename=paste("Length_histograms/Length_histogram.", sample.name, ".pdf", sep=""), 
             plot=myplot,
             width = 5, height = 4, units = "in"))
    }
    
  }
  
  if ("CSV" %in% filesToGenerate){
    # Save the histogram in a CSV format
    write.csv(data.frame(Length=1:1000, Percentage=100*h$density), 
              file=paste("Length_histograms/Length_histogram.", sample.name, ".csv", sep=""), 
              row.names=FALSE)
  }
  
  if ("RDATA" %in% filesToGenerate){
    fragmentLength = 1:1000
    Percentage = 100*h$density
    save(fragmentLength, Percentage, TotalNoReads, sample.name, file=paste("Length_histograms/Length_histogram.", sample.name, ".RData", sep=""))
  }
}