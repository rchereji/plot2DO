io <- file.path(sourceBasePath, "plot2DO_io.R")
data <- file.path(sourceBasePath, "plot2DO_data.R")
core <- file.path(sourceBasePath, "plot2DO_core.R")
plots <- file.path(sourceBasePath, "plot2DO_plots.R")
source(io)
source(data)
source(core)
source(plots)


Main <- function(command_line_args=NA)
{
  opt <- LoadArguments(command_line_args)  
  
  params <- InitializeParams(opt)

  rm(opt) # deleted to be sure using initialized params
  
  annotations <- LoadGenomeAnnotation(params$genome)
  
  rawReads <- LoadReads(params$inputFilePath, params$genome, annotations)
  
  reads <- CleanReads(rawReads, annotations$chrLen, params$lMin, params$lMax)

  # alignment:
  referenceGRanges <- Align(params$referencePointsBed, annotations, params$reference, 
                            params$beforeRef, params$afterRef, params$genome, params$align)
  
  outputFolderPath <- CreateOutputFolders(params$referencePointsBed, params$reference, params$siteLabel)
  
  CalculatePlotData(params, reads, referenceGRanges)

  plot <- PlotFigure(params)
  
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
  siteLabel <- opt$siteLabel
  
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
