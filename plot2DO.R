#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

# Load default paths from config.yaml
config <- yaml.load_file("config/config.yaml")

sourceBasePath <- config$application$paths$source
if(is.null(sourceBasePath)) {
  sourceBasePath <- file.path(getwd(), "source")
} 

readsBasePath <- config$application$paths$reads
if(is.null(readsBasePath)) {
  readsBasePath <- getwd()
}

outputBasePath <- config$application$paths$outTests
if(is.null(outputBasePath)) {
  outputBasePath <- file.path(getwd(), "output")
}

annotationsBasePath <- config$application$paths$annotations
if(is.null(annotationsBasePath)) {
  annotationsBasePath <- getwd()
}

configBasePath <- config$application$paths$config
if(is.null(configBasePath)) {
  configBasePath <- file.path(getwd(), "config")
}

main <- file.path(sourceBasePath, "plot2DO_main.R")  

source(main)

Main()

