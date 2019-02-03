# libraries:

# biocmanager:
# library(GenomicRanges)
# library(IRanges)
# library(rtracklayer)
# library(AnnotationHub)
# library(biomaRt)
# library(Rsamtools)

# cran:
# library(yaml)
# library(optparse)
# library(ggplot2)
# library(reshape2)
# library(colorRamps)
# library(gridExtra)

# R core:
# library(tools)
# library(grid)


# TODO: download required data (bam files)???


#TODO: create empty folders for annotations and bam files:



biocManagerLibraries <- c("GenomicRanges", "IRanges", "rtracklayer", "AnnotationHub",
                          "biomaRt", "Rsamtools")

cranLibraries <- c("yaml", "optparse", "ggplot2", "reshape2", "colorRamps", "gridExtra")


if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
sapply(biocManagerLibraries, function(x) {
  if (requireNamespace(x, quietly = TRUE)){
    result <- "already installed"
  } else {
    BiocManager::install(x)
    success <- requireNamespace(x, quietly = TRUE)
    if(success) { 
      result <- "successful instalation"
    } else {
      result <- "error"
    }    
  }  
  return(result)
})

sapply(cranLibraries, function(x) {
  if (requireNamespace(x, quietly = TRUE)){
    result <- "already installed"
  } else {
    install.packages(x)
    success <- requireNamespace(x, quietly = TRUE)
    if(success) { 
      result <- "successful instalation"
    } else {
      result <- "error"
    }
  }  
  return(result)
})



