suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(colorRamps)
  library(grid)
  library(gridExtra)
  library(yaml)
})

#debug parameters
# source("/home/paulati/Documents/ingebi/2018/razvan/plot2DO-develop/config.R")
# io <- file.path(sourceBasePath, "plot2DO_io.R")
# source(io)
# opt <- LoadArguments(arguments_1)  
# fin debug parameters

GetPlotsConfig <- function(plotType, selectedReference, siteLabel, align) {
  
  configFilePath <- file.path(configBasePath, "plot_config.yaml")
  config <- yaml.load_file(configFilePath)
  
  plotConfigType <- config$plot[[plotType]] 
  if(is.null(plotConfigType)) {
    print("error config")  
  } else {
    
    plotConfigSite <- plotConfigType[[siteLabel]] 
    if(is.null(plotConfigSite)) {
      
      plotConfigReference <- plotConfigType[[selectedReference]]
      if(is.null(plotConfigReference)) {
        print ("config error")
      } else {
        heatmapConfig <- plotConfigReference$heatmap
        averageConfig <- plotConfigReference$average
        fragmentLengthConfig <- plotConfigReference$fragmentLength
      }
      
    } else {
      plotConfigAlign <- plotConfigSite[[align]]
      if(is.null(plotConfigAlign)) {
        print ("config error")
      } else  {
        heatmapConfig <- plotConfigAlign$heatmap
        averageConfig <- plotConfigAlign$average
        fragmentLengthConfig <- plotConfigAlign$fragmentLength
      }
    }
  }
  result <- list(heatmap = heatmapConfig, average = averageConfig,  fragmentLength = fragmentLengthConfig)
  return(result)
  
}

GetPlotTitles <- function(plotConfig, sampleName) {
  
  mainTitle <- sampleName
  if(is.null(plotConfig)) { #set defaults
    xTitle <- ""
    yTitle <- ""
    legendTitle <- ""
  } else { 
    xTitle <- plotConfig$xLabel
    if(is.null(xTitle)) { xTitle <- "" }
    yTitle <- plotConfig$yLabel
    if(is.null(yTitle)) { yTitle <- ""  }
    legendTitle <- plotConfig$legend
    if(is.null(legendTitle)) { legendTitle <- ""  }
  }  
  result <- list(xTitle = xTitle, yTitle = yTitle, mainTitle = mainTitle, legendTitle = legendTitle)
  return(result)  
}

GetMinorTicksAxisLabels <- function(labels, gapLabelsLength) {
  result <- c()
  i <- 1    
  while(i <= length(labels)) {
    label <- labels[i]
    toProcessLength <- length(labels) - length(result)
    if(toProcessLength > gapLabelsLength) {
      effectiveGapLength <- gapLabelsLength  
    } else {
      effectiveGapLength <- toProcessLength - 1 # 1 is for label
    }
    gap <- rep("", effectiveGapLength)
    result <- c(result, label, gap)
    i <- i + effectiveGapLength + 1
  }

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
                        scaleXPosition, scaleYPosition) {
  
  # interpolate = TRUE suaviza los pixels
  occMatrixMelt <- melt(t(occMatrix))
  
  # TODO: pasar a configuracion
  legend.ticks.count <- 5
  breaks.size <- round((max(occMatrixMelt$value) - min(occMatrixMelt$value)) / legend.ticks.count, 3)
  breaks <- c(0:(legend.ticks.count - 1)) * breaks.size
  labels <- 100 * breaks # relative %
  
  result <- ggplot(occMatrixMelt, aes(Var1, Var2, fill = value)) + 
    geom_raster(aes(fill = value), interpolate = TRUE) + 
    scale_color_gradientn(colors = matlab.like(100), aesthetics = "fill", 
                          breaks = breaks, labels = labels) # + plotTheme 
  
  xbreaks.num <- 20
  step <- (afterRef + beforeRef) / xbreaks.num
  xBreaks <- seq(-beforeRef, afterRef, by=step) 
  matrixIndexes <- xBreaks + rep(beforeRef, length(xBreaks))
  xLabels <- GetMinorTicksAxisLabels(as.character(xBreaks), 4)  # 4 sin etiquetas
  xBreaks <- matrixIndexes
  xLimits <- c(min(xBreaks), max(xBreaks))
  
  yBreaks <- min(occMatrixMelt$Var2) + seq(min(occMatrixMelt$Var2) - 1, max(occMatrixMelt$Var2) - 1, by=10)
  yLabels <- GetMinorTicksAxisLabels(as.character(seq(lMin, lMax, by=10)), 4)  # 4 sin etiquetas
  yLimits <- c(min(yBreaks), max(yBreaks))
  
  guideColourbar <- guide_colourbar(title = legendTitle, reverse = FALSE, title.position = "left",
                                    frame.colour = "black", frame.linewidth = 0.5,
                                    ticks.colour = "black", ticks.linewidth = 0.5, 
                                    label.position = "left", 
                                    title.theme = legendTitleTheme,
                                    label.theme = legendLabelTheme)  
  
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
  
  yBreaks <- seq(0, 1.1*max(avgOcc), 0.2)
  yLabels <- GetMinorTicksAxisLabels(as.character(yBreaks), 1)  # yBreaks
  yLimits <- c(0, 1.1*max(avgOcc))
  xBreaks <- seq(-beforeRef, afterRef, 100)
  xLabels <- GetMinorTicksAxisLabels(as.character(xBreaks), 4) # xBreaks
  xLimits <- c(-beforeRef, afterRef)
  
  xlineInterceps <- c(-500, 0, 500)
  
  result <- ggplot(data=avgOcc.df) + 
    aes(x = avgOcc.df$x, y = avgOcc.df$avgOcc) + 
    geom_line(color="blue") + 
    # plotTheme + 
    scale_y_continuous(breaks = yBreaks, labels = yLabels, limits = yLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    scale_x_continuous(breaks = xBreaks, labels = xLabels, limits = xLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    geom_vline(xintercept = xlineInterceps, linetype = 'longdash', color = "lightgray", size = 0.4) +
    labs(x = xTitle, y = yTitle, title = mainTitle) + customTheme
  
  # result <- result + theme(plot.background = element_rect(fill = "red"))
  
  return(result)
  
}

PlotFragmentLength <- function(lengthHist, lMin, lMax, xTitle, yTitle, customTheme) {
  
  data <- as.data.frame(lengthHist)
  data$x <- seq(lMin, lMax, 1)
  
  xBreaks <- seq(lMin, lMax, 10)
  xLabels <- GetMinorTicksAxisLabels(as.character(xBreaks), 4) # xBreaks
  xLimits <- c(lMin, lMax)
  
  yBreaks <- round(seq(min(lengthHist), max(lengthHist), 1), 2)
  yLabels <- yBreaks
  # yLimits <- c(min(lengthHist), 1.05 * max(lengthHist))
  yLimits <- c(min(lengthHist),  max(lengthHist))
  
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

GetHeatmapLegend<-function(myggplot){
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
  
  dataFilePath <- getOutputMatrixFilePath(params$plotType, params$referencePointsBed, 
                                          params$reference, params$siteLabel, 
                                          params$lMin, params$lMax, params$sampleName)
  # load data:
  load(dataFilePath)
  # load occMatrix + histogram + ...
  
  # paula TODO: ver por que hay valores negativos:
  occMatrix[occMatrix < 0] = 0 # eliminate rounding errors
  
  # paula TODO: ???
  if (! is.null(params$colorScaleMax)) {
    occMatrix[occMatrix >= params$colorScaleMax] = params$colorScaleMax - 1e-10   # set a maximum threshold for the 2D Occ matrix
  }
  
  plotsConfig <- GetPlotsConfig(params$plotType, params$reference, siteLabel, params$align)

  heatmapTitles <- GetPlotTitles(plotsConfig$heatmap, sampleName)
  
  graphicalParams <- GetGraphicalParams(params$simplifyPlot, params$squeezePlot)
  
  heatmap <- PlotHeatmap(occMatrix, heatmapTitles$xTitle, heatmapTitles$yTitle, 
                         heatmapTitles$mainTitle, heatmapTitles$legendTitle, 
                         beforeRef, afterRef, lMin, lMax, 
                         graphicalParams$heatmapTheme, graphicalParams$legendTitleTheme,
                         graphicalParams$legendLabelTheme,
                         graphicalParams$scaleXPosition, graphicalParams$scaleYPosition)

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
    
    averageTitles <- GetPlotTitles(plotsConfig$average, sampleName)
    avgOccupancy <- PlotAverageOccupancy(occMatrix, beforeRef, afterRef, 
                                         averageTitles$xTitle, averageTitles$yTitle, averageTitles$mainTitle,
                                         graphicalParams$avgOccupancyTheme) 
    
    fragmentLengthTitles <- GetPlotTitles(plotsConfig$fragmentLength, sampleName)
    fragmentLength <- PlotFragmentLength(lengthHist, lMin, lMax, 
                                         fragmentLengthTitles$xTitle, fragmentLengthTitles$yTitle,
                                         graphicalParams$fragmentLengthTheme)
    
    #separate legend from heatmap:
    legend <- GetHeatmapLegend(heatmap)
    heatmapWithoutLegend <- heatmap + theme(legend.position="none")
  
    # showGrob(legend)

    result.grob <- arrangeGrob( #avgOccupancyTable, heatmapTable, 
                               avgOccupancy, heatmapWithoutLegend, 
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
         units = "cm", dpi = 300, scale = 1)   
  
}

GetGraphicalParams <- function(simplifyPlot, squeezePlot) {

  baseTheme <- theme_bw() + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5,  vjust = 0.5),
    axis.title = element_text(size = 8, hjust = 0.5, vjust = 0.5),
    #axis.title.x = element_text(size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 8, angle = 0, hjust = 0.5,  vjust = 0.5),
    #axis.title.y = element_text(size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.ticks = element_line(size = 0.2)
  )
  
  # TODO: move to config ?
  if(simplifyPlot) {
    if(squeezePlot) { 
      layout <- cbind(c(1, 1, 1), c(1, 1, 1))
      gridWidths <- c(1, 1)
      gridHeights <- c(1, 1, 1)
      plotWidth <- 14 # 5.5 
      plotHeight <- 20
    } else {
      layout <- cbind(c(1))
      gridWidths <- c(1)
      gridHeights <- c(1)
      plotWidth <- 18
      plotHeight <- 15
    }
    
    scaleXPosition <- "bottom"
    scaleYPosition <- "left"
    
    heatmapTheme <- baseTheme + theme(legend.position="right",
                                      axis.line.x.top = element_line(color="white"),
                                      axis.ticks.x.top = element_line(color="white"),
                                      axis.text.x.top = element_text(color="white"),
                                      axis.title.x.top = element_text(color="white"),
                                      axis.line.y.right = element_line(color="white"),
                                      axis.ticks.y.right = element_line(color="white"),
                                      axis.text.y.right = element_text(color="white"),
                                      axis.title.y.right = element_text(color="white"))
    avgOccupancyTheme <- NA
    fragmentLengthTheme <- NA
    
  } else {
    
    layout <- cbind(c(NA, 4, 4), c(1,2,2), c(NA,3,3))  
    gridWidths <- c(0.25, 2, 1)
    gridHeights <- c(1, 1, 1)
    plotWidth <- 21 
    plotHeight <- 18
    
    scaleXPosition <- "top"
    scaleYPosition <- "right"
    
    # i have to duplicate axis and add element text to align all the plots.
    # they should use tha space but not be shown
    heatmapTheme <- baseTheme + theme(legend.position="left",
                                      plot.title = element_blank(), 
                                      axis.title.x = element_text(),
                                      axis.title.y = element_text(),
                                      axis.line.x.bottom = element_line(color="white"),
                                      axis.ticks.x.bottom = element_line(color="white"), 
                                      axis.text.x.top = element_text(color="white"),
                                      axis.title.x.top = element_text(color="white"),
                                      axis.title.x.bottom = element_text(color="white"),
                                      axis.line.y.left = element_line(color="white"),
                                      axis.ticks.y.left = element_line(color="white"),
                                      axis.text.y.left = element_text(color="white"),
                                      axis.title.y.left = element_text(color="white"),
                                      axis.title.y.right = element_text(color="white"))
    
    avgOccupancyTheme <- baseTheme + theme( axis.line.x.top = element_line(color="white"),
                                            axis.ticks.x.top = element_line(color="white"), 
                                            axis.text.x.top = element_text(color="white"),
                                            axis.title.x.top = element_text(colour = "white"),
                                            axis.line.y.right = element_line(color="white"),
                                            axis.ticks.y.right = element_line(color="white"),
                                            axis.text.y.right = element_text(color="white"),
                                            axis.title.y.right = element_text(colour = "white"))
    
    fragmentLengthTheme <- baseTheme + theme( axis.line.x.top = element_line(color="white"),
                                             axis.ticks.x.top = element_line(color="white"), 
                                             axis.text.x.top = element_text(color="white"),
                                             axis.title.x.top = element_text(colour = "white"),
                                             axis.line.y.right = element_line(color="white"),
                                             axis.ticks.y.right = element_line(color="white"),
                                             axis.text.y.right = element_text(color="white"),
                                             axis.title.y.right = element_text(colour = "white"))
  }  

  legendTitleTheme <- element_text(size = 8, angle = 90, hjust = 0.5)
  legendLabelTheme <- element_text(size = 8, angle = 0)
  
      
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
