# Libraries to install:
# Bioconductor: GenomicRanges, IRanges, rtracklayer, AnnotationHub, biomaRt, Rsamtools
# CRAN: yaml, optparse, ggplot2, reshape2, colorRamps, gridExtra

# R core libraries (do not need installation): tools, grid


bioconductorLibraries <- c("GenomicRanges", "IRanges", "rtracklayer", "AnnotationHub",
                          "biomaRt", "Rsamtools")

cranLibraries <- c("yaml", "optparse", "ggplot2", "reshape2", "colorRamps", "gridExtra", "pracma")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
sapply(bioconductorLibraries, function(x) {
  if (requireNamespace(x, quietly = TRUE)){
    result <- "already installed"
  } else {
    BiocManager::install(x)
    success <- requireNamespace(x, quietly = TRUE)
    if(success) { 
      result <- "successful installation"
    } else {
      result <- "error"
    }    
  }  
  return(result)
})

# Install CRAN packages
sapply(cranLibraries, function(x) {
  if (requireNamespace(x, quietly = TRUE)){
    result <- "already installed"
  } else {
    install.packages(x, repos="http://cloud.r-project.org")
    success <- requireNamespace(x, quietly = TRUE)
    if(success) { 
      result <- "successful installation"
    } else {
      result <- "error"
    }
  }  
  return(result)
})
