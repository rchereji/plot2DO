GetPlotsLabels <- function(sampleName, type, reference, site, align) {
  
  #defaults:
  referenceLabel <- ""
  typeLabel <- ""
  alignLabel <- ""
  siteLabel <- ""
  
  if(type == "OCC") {
    
    typeLabel <- "occupancy"
  
  } else if(type == "DYADS") {
  
    typeLabel <- "dyad density"
    
  } else if(typeLabel == "fivePrime_ends") {
    
    typeLabel <- "5’ ends density"
    
  } else if(typeLabel == "threePrime_ends") {
    
    typeLabel <- "3’ ends density"
    
  }
  
  if(reference %in% c("TSS", "TTS")) {
    
    referenceLabel <- reference
    
  } else if(reference == "Plus1") {
    
    referenceLabel <- "+1 nuc."
    
  } else if(site != "" & align != "") {
  
    siteLabel <- site
    
    if(align == "threePrime") {
      alignLabel <- "3'"
    } else if( align == "fivePrime") {
      alignLabel <- "5'"
    } else if(align == "center") {
      alignLabel <- "center"
    }
  }
    
  if(reference.Label != "") {
  
    heatmap.xLabel <- "Fragment length (bp)"
    heatmap.yLabel <- paste0("Position relative to ", referenceLabel, " (bp)")
    heatmap.legend <- "Relative coverage (%)"
  
    fragmentLength.xLabel <- "Fragment length (bp)"
    fragmentLength.yLabel <- "Percentage (%)"
  
    average.xLabel <- paste0("Position relative to",  referenceLabel, " (bp)")
    average.yLabel <- paste0("Average ", typeLabel)
    
  } else if(siteLabel != "" & alignLabel != "") {
  
    heatmap.xLabel <- "Fragment length (bp)"
    heatmap.yLabel <- paste0("Position relative to ", alignLabel, " of ", siteLabel, " (bp)")
    heatmap.legend <- "Relative coverage (%)"
    
    fragmentLength.xLabel <- "Fragment length (bp)"
    fragmentLength.yLabel <- "Percentage (%)"
    
    average.xLabel <- paste0("Position relative to ",  alignLabel, " of ", siteLabel, " (bp)")
    average.yLabel <- paste0("Average ", typeLabel)
    
  }
  
  heatmapLabels <- list(xTitle = heatmap.xLabel, yTitle = heatmap.yLabel, legendTitle = heatmap.legend, mainTitle = sampleName)
  averageLabels <- list(xTitle = average.xLabel, yTitle = average.yLabel, legendTitle = "", mainTitle = "")
  fragmentLengthLabels <- list(xTitle = fragmentLength.xLabel, yTitle = fragmentLength.yLabel, legendTitle = "", mainTitle = "")
  
  result <- list(heatmap = heatmapLabels, average = averageLabels, fragmentLength = fragmentLengthLabels)

  return(result)
  
}