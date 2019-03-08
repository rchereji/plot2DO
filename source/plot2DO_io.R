suppressPackageStartupMessages({
  library(optparse)
  library(tools)
})

CreateOptionsParser <- function(){
  
  options = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
                help="Name of the file containing aligned sequencing data [options: BAM or BED file]"),
    make_option(c("-t", "--type"), type="character", default="occ", 
                help="Type of distribution to plot [options: occ, dyads, fivePrime_ends, threePrime_ends; default = %default]"),
    make_option(c("-g", "--genome"), type="character", default="sacCer3", 
                help="Genome version\n\t\t[options: sacCer3 (default) (S. cerevisiae); EF2 (S. pombe); dm3, dm6 (D. melanogaster);\n\t\tce10, ce11 (C. elegans); mm9, mm10 (M. musculus);\n\t\thg18, hg19, hg38 (H. sapiens); tair10 (A. thaliana)]"),
    make_option(c("-r", "--reference"), type="character", default="TSS", 
                help="Reference points to be aligned [options: TSS (default), TTS, Plus1]"),
    make_option(c("-s", "--sites"), type="character", default=NULL, 
                help="User-provided sites to be aligned (BED file)"),
    make_option(c("-a", "--align"), type="character", default="center", 
                help="Points of the provided intervals to be aligned? [options: center (default), fivePrime, threePrime]"),
    make_option(c("--siteLabel"), type="character", default="Sites", 
                help="Label for the aligned sites [default = %default]"),
    make_option(c("-l", "--minLength"), type="integer", default=50, 
                help="The smallest DNA fragment to be considered [default = %default]"),
    make_option(c("-L", "--maxLength"), type="integer", default=200, 
                help="The largest DNA fragment to be considered [default = %default]"),
    make_option(c("-u", "--upstream"), type="integer", default=1000, 
                help="Length of the upstream region to be plotted [default = %default]"),
    make_option(c("-d", "--downstream"), type="integer", default=1000, 
                help="Length of the downstream region to be plotted [default = %default]"),
    make_option(c("-m", "--colorScaleMax"), type="double", default=NULL, 
                help="Maximum value on the color scale (e.g. 0.02)"),
    make_option(c("--simplifyPlot"), type="character", default="off", 
                help="Simplify the plot (show only the 2D heat map) [options: on, off (default)]"),
    make_option(c("--squeezePlot"), type="character", default="off", 
                help="Simplify the plot and squeeze the heat map [options: on, off (default)]")
  )
  
  optParser = OptionParser(option_list=options)
  
  return(optParser)
}

LoadArguments <- function(commandLineArgs = NA){

  optParser <- CreateOptionsParser()
  
  if(anyNA(commandLineArgs) | length(commandLineArgs) == 0)  {
    opt <- parse_args(optParser)  
  } else {
    opt <- parse_args(optParser, args =  commandLineArgs)
  }
  
  if (is.null(opt$file)) {
    print_help(optParser)
    stop("At least the dataset file name must be supplied.", call.=FALSE)
  }
  
  if (opt$squeezePlot == "on"){
    opt$simplifyPlot = "on"
  }

  return(opt)
}

InitializeParams <- function(opt) {
  
  # Type of plot
  plotType <- toupper(opt$type)
  plotType <- strsplit(plotType, ',')[[1]]
  
  # Genome
  genome <- opt$genome
  
  inputFilePath <- opt$file
  
  # Data file name
  inputFilename <- basename(inputFilePath)
  
  # Reference label
  if(is.null(opt$sites)) {
      siteLabel <- ""
  } else {
      siteLabel <- opt$siteLabel
  }  
  
  if (! is.null(opt$colorScaleMax)) {
    colorScaleMax <- opt$colorScaleMax
  } else {
    colorScaleMax <- NULL
  }
  
  # Type of alignments, i.e. reference points
  selectedReference <- opt$reference
  
  # Size selection parameters: specify the interval of lengths to be analyzed
  lMin <- opt$minLength   # the smallest DNA fragment to be considered
  lMax <- opt$maxLength  # the largest DNA fragment to be considered
  
  # Window selection parameters
  beforeRef <- opt$upstream  # length of the upstream region that will be plotted
  afterRef <- opt$downstream   # length of the downstream region that will be plotted
  
  inputType = toupper(file_ext(inputFilename))
  
  sampleName <- switch(inputType,
                       BED={
                         sampleName <- sub(".bed", "", inputFilename)
                       },
                       BAM={
                         sampleName <- sub(".bam", "", inputFilename)
                       })
  
  align <- opt$align
  
  referencePointsBed <- opt$sites
  
  if(opt$squeezePlot == "off") {
    squeezePlot <- FALSE
  } else {
    squeezePlot <- TRUE
  }
  
  if(opt$simplifyPlot == "on") {
    simplifyPlot <- TRUE
  } else {
    simplifyPlot <- FALSE
  }
  
  # if (opt$squeezePlot == "on") => opt$simplifyPlot = "on"
  if (squeezePlot){
    simplifyPlot <- TRUE
  }
  
  result <- list(plotType = plotType, genome = genome, align = align,
                 inputFilename = inputFilename, inputFilePath = inputFilePath, 
                 sampleName = sampleName,
                 reference = selectedReference, siteLabel = siteLabel,
                 referencePointsBed = referencePointsBed, 
                 lMin = lMin, lMax = lMax, 
                 beforeRef = beforeRef, afterRef = afterRef,
                 colorScaleMax = colorScaleMax,
                 squeezePlot = squeezePlot, simplifyPlot = simplifyPlot)
  
  return(result)
}

CreateOutputFolders <- function(plotType, sites, selectedReference, siteLabel) {
  # Create folder where figures are saved
  folderPath <- GetOutputBaseFolderPath(plotType, sites, selectedReference, siteLabel)
  dir.create(folderPath, showWarnings = FALSE, recursive = TRUE)
  return(folderPath)
}

# output helper functions:
GetOutputBaseFolderPath <- function(plotType, sites, selectedReference, siteLabel){

  folderNameHead <- paste0("2D_", tolower(plotType), "_")
  
  if (is.null(sites)) {
    folderName <- paste(folderNameHead, selectedReference, sep = "")
  } else {
    folderName <- paste(folderNameHead, siteLabel, sep="")
  }
  
  folderPath <- file.path(outputBasePath, folderName)
  
  return(folderPath)
}

GetOutputPlotFilePath <- function(plotType, sites, selectedReference, siteLabel,
                                      lMin, lMax, sampleName)  {
  extension <- ".pdf"
  result <- GetOutputFilePathBase(plotType, sites, selectedReference, siteLabel, 
                                      lMin, lMax, sampleName, extension)
  return(result)
  
}

GetOutputMatrixFilePath <- function(plotType, sites, selectedReference, siteLabel,
                                        lMin, lMax, sampleName){
  
  extension <- ".RData"
  
  result <- GetOutputFilePathBase(plotType, sites, selectedReference, siteLabel, 
                                        lMin, lMax, sampleName, extension)
  
  return(result)
  
}

GetOutputFilePathBase <- function(plotType, sites, selectedReference, siteLabel, 
                                      lMin, lMax, sampleName, extension){
  
  folder_base_path <- GetOutputBaseFolderPath(plotType, sites, selectedReference, siteLabel)
  
  fileNameHead <- paste0(plotType, "_matrix.")
  fileNameTail <- paste0(".", lMin, "_", lMax, ".", sampleName, extension)
  
  if (is.null(sites)) {
    fileName <- paste0(fileNameHead, selectedReference, fileNameTail)
  } else {
    fileName <- paste0(fileNameHead, siteLabel, fileNameTail)
  }
  
  filePath <- file.path(folder_base_path, fileName)
  
  return(filePath)
  
}
