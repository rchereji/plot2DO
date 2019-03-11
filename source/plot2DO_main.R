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
  
  # annotations <- LoadGenomeAnnotation(params$genome)
  annotations <- LoadGenomeAnnotation_v2(params$inputFilePath, params$genome)
  
  rawReads <- LoadReads(params$inputFilePath, params$genome, annotations)
  
  reads <- CleanReads(rawReads, annotations$chrLen, params$lMin, params$lMax)

  referenceGRanges <- ConstructReferenceGRanges(params$referencePointsBed, annotations, params$reference, 
                            params$beforeRef, params$afterRef, params$genome, params$align)
  
  outputFolderPath <- CreateOutputFolders(params$plotType, params$referencePointsBed, params$reference, params$siteLabel)
  
  CalculatePlotData(params, reads, referenceGRanges)

  plot <- PlotFigure(params)
    
  # check whether the unwanted file exists and remove it
  somethingToRemove <- file.exists("Rplots.pdf")
  if(somethingToRemove) {
      dummy <- file.remove("Rplots.pdf")  
  }
  
}
