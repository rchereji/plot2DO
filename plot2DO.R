#!/usr/bin/env Rscript

# v1.0.1 - add Arabidopsis thaliana (TAIR10)
#        - change option --organism to --genome
# v1.0.2 - add S. pombe (EF2)
# v1.0.3 - add chr Y and rDNA for dm6 genome
# v1.0.4 - correct RangesList issue (see https://support.bioconductor.org/p/109079/)
# v1.0.5 - refactor

# load app config values:

suppressPackageStartupMessages({
  library(yaml)
})

config <- yaml.load_file("config/config.yaml")

sourceBasePath <- config$application$paths$source
if(is.null(sourceBasePath)) {
  sourceBasePath <- file.path(getwd(), "source")
} 

readsBasePath <- config$application$paths$reads
if(is.null(readsBasePath)) {
  #readsBasePath <- file.path(getwd(), "data")
  readsBasePath <- file.path(getwd(), "Sample_BAM_files")
}

testOutputBasePath <- config$application$paths$outTests
if(is.null(testOutputBasePath)) {
  testOutputBasePath <- file.path(getwd(), "tests")
}

annotationsBasePath <- config$application$paths$annotations
if(is.null(annotationsBasePath)) {
  # annotationsBasePath <- file.path(getwd(), "Annotations")
  annotationsBasePath <- getwd()
}

configBasePath <- config$application$paths$config
if(is.null(configBasePath)) {
  configBasePath <- file.path(getwd(), "config")
}

main <- file.path(sourceBasePath, "plot2DO_main.R")  

source(main)

Main()

