suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(colorRamps)
  library(grid)
  library(gridExtra)
  library(yaml)
})

# Construct all plot labels and titles according to the type of plot/alignment
GetPlotsLabels <- function(sampleName, type, reference, site, align) {

  typeLabel <- switch(type, 
                      "OCC"            = "occupancy", 
                      "DYADS"          = "dyad density", 
                      "FIVEPRIME_ENDS" = "density", 
                      "THREEPRIME_ENDS"= "density", 
                      "")
  
  referenceLabel <- switch(reference, 
      "TSS"  = "TSS", 
      "TTS"  = "TTS", 
      "Plus1"= "+1 nuc.", 
      "")

  
  if (site != "" & align != "") {
    siteLabel <- site
    alignLabel <- switch(align, 
                    "fivePrime"  = "5'", 
                    "threePrime" = "3'", 
                    "center"= "center", 
                    "")
  } else {      
      alignLabel <- ""
      siteLabel <- ""
  }

  # Construct axes labels
  heatmapYLabel <- "Fragment length (bp)"
  heatmapLegend <- "Relative coverage (%)"
  
  fragmentLengthXLabel <- "Fragment length (bp)"
  fragmentLengthYLabel <- "Percentage (%)"
  
  if (siteLabel != "" & alignLabel != "") {
    
    heatmapXLabel <- paste0("Position relative to ", alignLabel, " of ", siteLabel, " (bp)")
    averageXLabel <- paste0("Position relative to ", alignLabel, " of ", siteLabel, " (bp)")
    averageYLabel <- paste0("Average ", typeLabel)
    
  } else if (referenceLabel != "") {

    heatmapXLabel <- paste0("Position relative to ", referenceLabel, " (bp)")
    averageXLabel <- paste0("Position relative to ", referenceLabel, " (bp)")
    averageYLabel <- paste0("Average ", typeLabel)
    
  }
      
  heatmapLabels <- list(xTitle = heatmapXLabel, yTitle = heatmapYLabel, legendTitle = heatmapLegend, mainTitle = sampleName)
  averageLabels <- list(xTitle = averageXLabel, yTitle = averageYLabel, legendTitle = "", mainTitle = sampleName)
  fragmentLengthLabels <- list(xTitle = fragmentLengthXLabel, yTitle = fragmentLengthYLabel, legendTitle = "", mainTitle = "")
  
  result <- list(heatmap = heatmapLabels, average = averageLabels, fragmentLength = fragmentLengthLabels)
  
  return(result)
  
}

# Remove some of the labels, i.e. skip noSkippedTicks between consecutive labels
FixTickLabels <- function(labels, noSkippedTicks, paddingLeft) {
  newLabels <- paste0(paddingLeft, labels)
  labelsToSkip <- setdiff(1:length(labels), seq(1, length(labels), noSkippedTicks + 1))
  newLabels[labelsToSkip] <- ""
  return(newLabels)
}

GetHeatmapBreaksAndLabels <- function(occMatrix, colorScaleMax) {

    if (is.null(colorScaleMax)) {
      maxValue <- max(occMatrix)
    } else {
      maxValue <- colorScaleMax      
    }

    if (maxValue > 0.1) {
      step <- 0.05
    } else if (maxValue > 0.05) {
      step <- 0.01
    } else {
      step <- 0.005
    }
        
    breaks <- seq(0, maxValue, step)
    labels <- breaks * 100 # Use percentages instead of fractional numbers  
    limits <- c(0, maxValue)
  
    result <- list(breaks = breaks, labels = labels, limits = limits)
    
    return(result)
}

# debug:
# xTitle <- heatmapTitles$xTitle
# yTitle <- heatmapTitles$yTitle
# mainTitle <- heatmapTitles$mainTitle
# legendTitle <- heatmapTitles$legendTitle 
# customTheme <- graphicalParams$heatmapTheme
# legendTitleTheme <- graphicalParams$legendTitleTheme
# legendLabelTheme <- graphicalParams$legendLabelTheme
# scaleXPosition <- graphicalParams$scaleXPosition
# scaleYPosition <- graphicalParams$scaleYPosition

PlotHeatmap <- function(occMatrix, xTitle, yTitle, mainTitle, legendTitle,
                        beforeRef, afterRef, lMin, lMax, 
                        customTheme, legendTitleTheme, legendLabelTheme, 
                        scaleXPosition, scaleYPosition, colorScaleMax) {
  
  # interpolate = TRUE suaviza los pixels
  occMatrixMelt <- melt(t(occMatrix))
  
  # TODO: pasar a configuracion
  #legend.ticks.count <- 5
  #breaks.size <- round((max(occMatrixMelt$value) - min(occMatrixMelt$value)) / legend.ticks.count, 3)
  #breaks <- c(0:(legend.ticks.count - 1)) * breaks.size
  #labels <- 100 * breaks # relative %
  
  breaksLabels <- GetHeatmapBreaksAndLabels(occMatrix, colorScaleMax)
    
  result <- ggplot(occMatrixMelt, aes(Var1, Var2, fill = value)) + 
    geom_raster(aes(fill = value), interpolate = TRUE) + 
    scale_color_gradientn(colors = matlab.like(100), aesthetics = "fill", 
                          breaks = breaksLabels$breaks, 
                          labels = breaksLabels$labels,
                          limits = breaksLabels$limits)
  
  xbreaks.num <- 20
  step <- (afterRef + beforeRef) / xbreaks.num
  xBreaks <- seq(-beforeRef, afterRef, by=step) 
  matrixIndexes <- xBreaks + rep(beforeRef, length(xBreaks))
  xLabels <- FixTickLabels(as.character(xBreaks), 4, "")  # skip 4 labels
  xBreaks <- matrixIndexes
  xLimits <- c(min(xBreaks), max(xBreaks))
  
  yBreaks <- min(occMatrixMelt$Var2) + seq(min(occMatrixMelt$Var2) - 1, max(occMatrixMelt$Var2) - 1, by=10)
  yLabels <- FixTickLabels(as.character(seq(lMin, lMax, by=10)), 4, "")  # skip 4 labels
  yLimits <- c(min(yBreaks), max(yBreaks))
  
  guideColourbar <- guide_colourbar(title = legendTitle, reverse = FALSE, title.position = "left",
                                    frame.colour = "black", frame.linewidth = 0.5,
                                    ticks.colour = "black", ticks.linewidth = 0.5, 
                                    label.position = "left", 
                                    title.theme = legendTitleTheme,
                                    label.theme = legendLabelTheme,
                                    title.hjust = 0.5,
                                    barheight = 10,                                    
                                    nbin = 1000)  # necessary to shift zero tick to bottom
  
  result <- result + labs(x = xTitle, y = yTitle, title = mainTitle)
  
  scaleY <- scale_y_continuous(breaks = yBreaks, labels = yLabels, limits = yLimits,
                               expand = c(0,0), position = scaleYPosition,
                               sec.axis = dup_axis())
  scaleX <- scale_x_continuous(breaks = xBreaks, labels = xLabels, limits = xLimits,
                               expand = c(0,0), position = scaleXPosition, 
                               sec.axis = dup_axis())

  result <- result + guides(fill = guideColourbar) +
    scaleY + scaleX + geom_vline(xintercept=beforeRef, linetype='longdash', color="white", size=0.4) +
    customTheme

  # result <- result + theme(plot.background = element_rect(fill = "darkblue"))
  
  return(result)

}  

PlotAverageOccupancy <- function(occMatrix, beforeRef, afterRef, xTitle, yTitle, mainTitle, customTheme) {
  
  avgOcc <- colSums(occMatrix)
  avgOcc.df <- as.data.frame(avgOcc)
  avgOcc.df$x <- seq(-beforeRef, afterRef, 1)
  
  deltaYBreaks = max(0.1, round(1.1*max(avgOcc)/10, digits = 1))
  yBreaks <- seq(0, 1.1*max(avgOcc), deltaYBreaks)
  if (length(yBreaks) >= 8){
    yLabels <- FixTickLabels(as.character(yBreaks), 1, " ")  # yBreaks
  } else {
    yLabels <- FixTickLabels(as.character(yBreaks), 0, " ")  # yBreaks
  }
  yLimits <- c(0, 1.1*max(avgOcc))
  xBreaks <- seq(-beforeRef, afterRef, 100)
  xLabels <- FixTickLabels(as.character(xBreaks), 4, "") # xBreaks
  xLimits <- c(-beforeRef, afterRef)
  
  xlineInterceps <- c(-500, 0, 500)
  
  result <- ggplot(data=avgOcc.df) + 
    aes(x = avgOcc.df$x, y = avgOcc.df$avgOcc) + 
    geom_line(color="blue") + 
    scale_y_continuous(breaks = yBreaks, labels = yLabels, limits = yLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    scale_x_continuous(breaks = xBreaks, labels = xLabels, limits = xLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    geom_vline(xintercept = xlineInterceps, linetype = 'longdash', color = "lightgray", size = 0.4) +
    labs(x = xTitle, y = yTitle, title = mainTitle) + customTheme

  return(result)
  
}

PlotFragmentLength <- function(lengthHist, lMin, lMax, xTitle, yTitle, customTheme) {
  
  data <- as.data.frame(lengthHist)
  data$x <- seq(lMin, lMax, 1)
  
  xBreaks <- seq(lMin, lMax, 10)
  xLabels <- FixTickLabels(as.character(xBreaks), 4, "") # skip 4 labels
  xLimits <- c(lMin, lMax)
  
  yBreaks <- round(seq(min(lengthHist), 1.05 * max(lengthHist), 1), 2)
  yLabels <- yBreaks
  yLimits <- c(min(lengthHist), 1.05 * max(lengthHist))
  
  xlineInterceps <- c(100, 150)
  
  result <- ggplot(data=data) + 
    aes(x = data$x, y = data$lengthHist) + 
    geom_line(color="blue") + 
    # plotTheme +
    scale_y_continuous(breaks = yBreaks, labels = yLabels, limits = yLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    scale_x_continuous(breaks = xBreaks, labels = xLabels, limits = xLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    geom_vline(xintercept = xlineInterceps, linetype = 'longdash', color = "lightgray", size = 0.4) +
    labs(x = xTitle, y = yTitle) + 
    customTheme
  
  result <- result + coord_flip()
  
  # result <- result + theme(plot.background = element_rect(fill = "yellow"))
  return(result)
  
}

GetHeatmapLegend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# TODO: consultar con razvan si leer datos guardados o pasar como params
# PlotFigure <- function(opt, occMatrix, lMin, lMax, 
#                        beforeRef, afterRef, sampleName,
#                        plotType, lengthHist)
PlotFigure <- function(params)  {

  suppressWarnings({
    
  dataFilePath <- GetOutputMatrixFilePath(params$plotType, params$referencePointsBed, 
                                          params$reference, params$siteLabel, 
                                          params$lMin, params$lMax, params$sampleName)

  # Load data:
  load(dataFilePath)

  occMatrix[occMatrix < 0] = 0 # eliminate rounding errors
  if (! is.null(params$colorScaleMax)) {
    occMatrix[occMatrix > params$colorScaleMax] = params$colorScaleMax # override the default colorbar and set everything above the threshold with the threshold value
  }
  
  plotsConfig <- GetPlotsLabels(sampleName, params$plotType, params$reference, siteLabel, params$align)
  
  heatmapTitles <- plotsConfig$heatmap
  
  graphicalParams <- GetGraphicalParams(params$simplifyPlot, params$squeezePlot)
  
  heatmap <- PlotHeatmap(occMatrix, heatmapTitles$xTitle, heatmapTitles$yTitle, 
                         heatmapTitles$mainTitle, heatmapTitles$legendTitle, 
                         beforeRef, afterRef, lMin, lMax, 
                         graphicalParams$heatmapTheme, graphicalParams$legendTitleTheme,
                         graphicalParams$legendLabelTheme,
                         graphicalParams$scaleXPosition, graphicalParams$scaleYPosition,
                         params$colorScaleMax)

  if(params$simplifyPlot) { 
    avgOccupancy <- NA
    fragmentLength <- NA
    result.grob <- arrangeGrob(heatmap, 
                               ncol = ncol(graphicalParams$layout), nrow = nrow(graphicalParams$layout),
                               layout_matrix = graphicalParams$layout, 
                               widths = graphicalParams$gridWidths, 
                               heights = graphicalParams$gridHeights)
    # grid.draw(result.grob)
  } else {
    
    averageTitles <- plotsConfig$average
    avgOccupancy <- PlotAverageOccupancy(occMatrix, beforeRef, afterRef, 
                                         averageTitles$xTitle, averageTitles$yTitle, averageTitles$mainTitle,
                                         graphicalParams$avgOccupancyTheme) 
    
    fragmentLengthTitles <- plotsConfig$fragmentLength
    fragmentLength <- PlotFragmentLength(lengthHist, lMin, lMax, 
                                         fragmentLengthTitles$xTitle, fragmentLengthTitles$yTitle,
                                         graphicalParams$fragmentLengthTheme)
    
    # Separate legend from heatmap:
    legend <- GetHeatmapLegend(heatmap)
    heatmapWithoutLegend <- heatmap + theme(legend.position="none")
  
    result.grob <- arrangeGrob(avgOccupancy, heatmapWithoutLegend, 
                               fragmentLength, legend, 
                               ncol = ncol(graphicalParams$layout), nrow = nrow(graphicalParams$layout), 
                               layout_matrix = graphicalParams$layout,
                               widths = graphicalParams$gridWidths, 
                               heights = graphicalParams$gridHeights)
    # grid.draw(result.grob)
  }

  plotFilePath <- GetOutputPlotFilePath(params$plotType, params$referencePointsBed, 
                                        params$reference, siteLabel, 
                                        lMin, lMax, sampleName)
  ggsave(plotFilePath, result.grob, width = graphicalParams$plotWidth, height = graphicalParams$plotHeight, 
         units = "in", dpi = 300, scale = 1)   
  
  })
  
  return(result.grob)
    
}

GetGraphicalParams <- function(simplifyPlot, squeezePlot) {

  baseTheme <- theme_bw() + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5,  vjust = 0.5, colour="black"),
    axis.text.y = element_text(size = 12, angle = 0, hjust = 0.5,  vjust = 0.5, colour="black"),
    axis.title = element_text(size = 14, hjust = 0.5, vjust = 0.5),
    axis.ticks.length = unit(.15, "cm")
  )
  
  legendTitleTheme <- element_text(size = 14, angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(t = 0, r = 7, b = 0, l = 0))
  legendLabelTheme <- element_text(size = 12, angle = 0, hjust = 0.5)
  
  # TODO: move to config ?
  if(simplifyPlot) {
    if(squeezePlot) { 
      layout <- cbind(c(1, 1, 1), c(1, 1, 1))
      gridWidths <- c(1, 1)
      gridHeights <- c(1, 1, 1)
      plotWidth <- 5.5
      plotHeight <- 8
    } else {
      layout <- cbind(c(1))
      gridWidths <- c(1)
      gridHeights <- c(1)
      plotWidth <- 7
      plotHeight <- 6
    }
    
    scaleXPosition <- "bottom"
    scaleYPosition <- "left"
    
    heatmapTheme <- baseTheme + theme(legend.position="right",
                                      axis.ticks.x.top = element_blank(),                                      
                                      axis.text.x.top = element_blank(),
                                      axis.title.x.top = element_blank(),
                                      axis.title.x.bottom = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                                      axis.ticks.y.right = element_blank(),                                      
                                      axis.text.y.right = element_blank(),
                                      axis.title.y.right = element_blank(),
                                      axis.title.y.left = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
    avgOccupancyTheme <- NA
    fragmentLengthTheme <- NA
    
  } else {
    
    layout <- cbind(c(NA, 4, 4), c(1,2,2), c(NA,3,3))  
    gridWidths <- c(0.5, 2, 1)
    gridHeights <- c(1.2, 1, 1)
    plotWidth <- 10
    plotHeight <- 7
    
    scaleXPosition <- "top"
    scaleYPosition <- "right"

    # i have to duplicate axis and add element text to align all the plots.
    # they should use that space but not be shown
    heatmapTheme <- baseTheme + theme(plot.title = element_blank(),                                                                           
                                      axis.ticks.x.bottom = element_blank(),                                      
                                      axis.ticks.y.left = element_blank(),
                                      axis.text.y.left = element_text(color="white"),                                      
                                      axis.text.x.bottom = element_text(color="white"),
                                      axis.title.x.bottom = element_text(color="white"),
                                      axis.title.x.top = element_blank(),                                      
                                      axis.title.y.left = element_text(color="white", margin = margin(t = 0, r = 10, b = 0, l = 0)),
                                      axis.title.y.right = element_blank())
    
    avgOccupancyTheme <- baseTheme + theme(axis.ticks.x.top = element_blank(),
                                           axis.text.x.top = element_blank(),
                                           axis.title.x.top = element_blank(),
                                           axis.title.x.bottom = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)),
                                           axis.ticks.y.right = element_blank(),
                                           axis.text.y.right = element_text(color="white"),
                                           axis.title.y.right = element_blank(),
                                           axis.title.y.left = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
    
    fragmentLengthTheme <- baseTheme + theme(axis.ticks.x.top = element_blank(),
                                             axis.text.x.top = element_text(color="white"),
                                             axis.title.x.top = element_blank(),                                             
                                             axis.ticks.y.right = element_blank(),
                                             axis.text.y.right = element_blank(),
                                             axis.title.y.right = element_blank(),
                                             axis.title.y.left = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
  }  

     
  result = list(layout = layout, 
                gridWidths = gridWidths, gridHeights = gridHeights, 
                plotWidth = plotWidth, plotHeight = plotHeight, 
                heatmapTheme = heatmapTheme,
                legendTitleTheme = legendTitleTheme,
                legendLabelTheme = legendLabelTheme, 
                avgOccupancyTheme = avgOccupancyTheme,
                fragmentLengthTheme = fragmentLengthTheme,
                scaleXPosition = scaleXPosition,
                scaleYPosition = scaleYPosition)
  
  return(result)
  
}
