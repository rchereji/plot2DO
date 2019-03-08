#!/usr/bin/env Rscript

# Check the minimum R version: R>=3.5.0 (required for Bioconductor 3.7)
getRversion <- function() {
  rvs <- sub("R version ", "", R.version.string)
  pos <- gregexpr(" \\(", rvs)[[1]][1]
  rVersion <- substr(rvs, 1, pos-1)
  releaseDate <- substr(rvs, pos+2, nchar(rvs)-1)
  result <- list("version" = rVersion, "date" = releaseDate)
  return(result)
}

rvs <- getRversion()
rVersion <- rvs$version
releaseDate <- rvs$date
if (rVersion < "3.5") {
  stop(paste0("Your R version is ", rVersion,", which is a bit too old (release date: ", releaseDate, ").\nTo run plot2DO, please upgrade R to a version >= 3.5.0.\nYou can download the latest version of R from https://www.r-project.org/"))
}

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

