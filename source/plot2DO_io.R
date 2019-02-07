suppressPackageStartupMessages({
  library("optparse")
})

CreateOptionsParser <- function(){
  
  options = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
                help="Dataset file name [options: BED or BAM format]"),
    make_option(c("-t", "--type"), type="character", default="occ", 
                help="Types of distribution to plot [options: occ, dyads, fivePrime_ends, threePrime_ends; default = %default]"),
    make_option(c("-g", "--genome"), type="character", default="sacCer3", 
                help="Genome version [options: sacCer3 (default) (S. cerevisiae); EF2 (S. pombe); tair10 (A. thaliana); dm3, dm6 (D. melanogaster), ce10, ce11 (C. elegans), mm9, mm10 (M. musculus), hg18, hg19 (H. sapiens)]"),
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

CreateOutputFolders <- function(sites, selectedReference, siteLabel) {
  # Create folder where figures are saved
  folderPath <- GetOutputBaseFolderPath(sites, selectedReference, siteLabel)
  dir.create(folderPath, showWarnings = FALSE, recursive = TRUE)
  return(folderPath)
}

SaveOccupancyMatrix <- function(occMatrix, sites, siteLabel, selectedReference, 
                                  folderPath, 
                                  lMin, lMax, beforeRef, afterRef, 
                                  totalNoReads, lengthHist, sampleName) {
  
  # Save the 2D Occupancy matrix
  if (is.null(sites)) {
    fileNameTail <- selectedReference
  } else  {
    fileNameTail <- siteLabel
  }  
    
  fileName <- paste("/threePrime_ends_matrix.", fileNameTail, ".", lMin, "_", lMax, ".", sampleName, ".RData", sep="")
  filePath <- file.path(folderPath, fileName)

  #save(occ_matrix, lMin, lMax, beforeRef, afterRef, TotalNoReads, LengthHist, sampleName, siteLabel, file=paste("2D_occupancy_", opt$siteLabel, "/threePrime_ends_matrix.", siteLabel, ".", lMin, "_", lMax, ".", sampleName, ".RData", sep=""))
  save(occMatrix, lMin, lMax, beforeRef, afterRef, totalNoReads, lengthHist, sampleName, siteLabel, file=filePath)
  
}


# output helper functions:

GetOutputBaseFolderPath <- function(sites, selectedReference, siteLabel){

  folderNameHead <- "2D_occupancy_"
  
  if (is.null(sites)) {
    folderName <- paste(folderNameHead, selectedReference, sep = "")
  } else {
    folderName <- paste(folderNameHead, siteLabel, sep="")
  }
  
  folderPath <- file.path(testOutputBasePath, folderName)
  
  return(folderPath)
}

GetOutputPlotFilePath <- function(plotType, sites, selectedReference, siteLabel,
                                      lMin, lMax, sampleName)  {
  extension <- ".pdf"
  result <- GetOutputFilePathBase(plotType, sites, selectedReference, siteLabel, 
                                      lMin, lMax, sampleName, extension)
  return(result)
  
}

getOutputMatrixFilePath <- function(plotType, sites, selectedReference, siteLabel,
                                        lMin, lMax, sampleName){
  
  extension <- ".RData"
  
  result <- GetOutputFilePathBase(plotType, sites, selectedReference, siteLabel, 
                                        lMin, lMax, sampleName, extension)
  
  return(result)
  
}

GetOutputFilePathBase <- function(plotType, sites, selectedReference, siteLabel, 
                                      lMin, lMax, sampleName, extension){
  
  folder_base_path <- GetOutputBaseFolderPath(sites, selectedReference, siteLabel)
  
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
