#' @title rbioFS_PCA
#'
#' @description A simple to use wrapper for PCA (Principal Component Analysis) and visualization
#' @param objTitle File name prefix for output files and the PCA object. Default is \code{"data"}.
#' @param input Input data, data frame.
#' @param sampleIDVar Sample variable name. It's a character string.
#' @param groupIDVar Group variable name. It's a character string.
#' @param scaleData If to scale the data when performing PCA. Default is \code{TRUE}.
#' @param boxplotTitle The boxplot title. Default is \code{NULL}.
#' @param boxplotWidth The boxplot width. Default is \code{170}.
#' @param boxplotHeight	The boxplot height. Default is \code{150}.
#' @param biplotPC The two PCs to plot for the biplot. Default is \code{c("PC1", "PC2")}.
#' @param biplotTitle The biplot title. Default is \code{NULL}.
#' @param biplotSymbolSize The symbol size for the scatter plot portion of the biplot. Default is \code{2}.
#' @param ellipse	If to draw ellipses. Default is \code{FALSE}.
#' @param ellipse_conf The confidence value for the ellipses. Default is \code{0.93}.
#' @param loadingPlot	If to superimpose loading plot. Default is \code{TRUE}.
#' @param loadingSize	The font size of the loading plot labels. Default is \code{3}.
#' @param biplotWidth	The biplot width. Default is \code{170}.
#' @param biplotHeight The biplot height. Default is \code{150}.
#' @param fontType Font for the figure texts.
#' @param xLabel_biplot	X-axis label for the biplot. Default is NULL.
#' @param yLabel_biplot	Y-axis label for the biplot. Default is NULL.
#' @return Outputs a PCA object, a boxplot (proportion of variance) and a biplot from PCA analysis. The format is \code{pdf}.
#' @details Make sure to arrange input data with first two columns for smaple ID and conditions, and the rest for features (e.g., genes).
#' @import ggplot2
#' @examples
#' \dontrun{
#' rbioFS_PCA(input = pcaDfm, idx = pcaDfm$Conditions, ellipse = TRUE, loadingPlot = TRUE, biplotWidth = 200, biplotHeight = 170)
#' }
#' @export
rbioFS_PCA <- function(objTitle = "data", input = NULL, sampleIDVar = NULL, groupIDVar = NULL, scaleData = TRUE,
                               boxplotTitle = NULL,
                               boxplotWidth = 170, boxplotHeight = 150,
                               biplotPC = c("PC1", "PC2"),
                               biplotTitle = NULL,
                               biplotSymbolSize = 2,
                               ellipse = FALSE, ellipse_conf = 0.93,
                               loadingPlot = TRUE, loadingSize = 3,
                               biplotWidth = 170, biplotHeight = 150,
                               fontType = "sans", xTickLblSize = 10, yTickLblSize = 10){


  ## PCA
  PCA <- prcomp(input[, ! names(input) %in% c(sampleIDVar, groupIDVar)], scale. = scaleData)

  ## plotting
  grid.newpage()

  # boxplot
  varpp <- 100 * summary(PCA)$importance[2, ] # extract and calcualte the proportion of variance
  boxdfm <- data.frame(PC = as.numeric(gsub("PC", "", names(varpp))), varpp = varpp)

  boxplt <- ggplot(data = boxdfm, aes(x = PC, y = varpp, group = 1)) +
    geom_bar(position = "dodge", stat = "identity", color = "black", fill = "gray66") +
    scale_x_continuous() +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, with(boxdfm, ceiling(max(varpp))))) +
    ggtitle(boxplotTitle) +
    xlab("PC") +
    ylab("Proportion of variance (%)") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          plot.title = element_text(face = "bold", family = fontType, hjust = 0.5),
          axis.title = element_text(face = "bold", family = fontType),
          legend.position = "bottom",legend.title = element_blank(),legend.key = element_blank(),
          axis.text.x = element_text(size = xTickLblSize, family = fontType),
          axis.text.y = element_text(size = yTickLblSize, family = fontType, hjust = 0.5))

  ggsave(filename = paste(objTitle,".PCA.boxplot.pdf", sep = ""), plot = boxplt,
         width = boxplotWidth, height = boxplotHeight, units = "mm",dpi = 600)

  grid.draw(boxplt) # preview


  # biplot. Inspired by ggord package.
  # argument check
  if (length(biplotPC) != 2){

    stop("Please properly set biplotPC argument so that two and only two PCs are used for biplot.")

  } else {
    # prepare for scatter plot values (i.e. sample score)
    sampleScore <- data.frame(PCA$x[, biplotPC], check.names = FALSE) # extract rotated sample scores
    sampleScore$Group <- factor(input[, groupIDVar], levels = unique(input[, groupIDVar]))
    names(sampleScore)[1:2] <- c("axis1", "axis2")

    # prepare for loading plot values (i.e. loading value for variables)
    loadingValue <- data.frame(PCA$rotation[, biplotPC], check.names = FALSE) # extract loading/rotation/eigenvectors for variables
    names(loadingValue) <- c("axis1", "axis2") # give a generic variable name for the ratation dataframe
    loadingScale <- max(max(abs(sampleScore$axis1)) / max(abs(loadingValue$axis1)),
                        max(abs(sampleScore$axis2)) / max(abs(loadingValue$axis2))) * 0.85 # determine scaling for loading values
    loadingValuePlot <- loadingValue * loadingScale
    loadingValuePlot$lbl <- rownames(loadingValuePlot)

    # perpare for the axis labels
    varpp_biplot <- varpp[biplotPC] # extract the proportion of variance for the selected PCs
    pc_axis_lbl <- paste(biplotPC, " (", round(varpp_biplot, digits = 2), "%)", sep = "")


    biplt <- ggplot(sampleScore, aes(x = axis1, y = axis2)) +
      geom_point(aes(shape = Group, colour = Group), size = biplotSymbolSize) + # plot the sample score scatter plot
      ggtitle(biplotTitle) +
      xlab(pc_axis_lbl[1]) +
      ylab(pc_axis_lbl[2]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = fontType, hjust = 0.5),
            axis.title = element_text(face = "bold", family = fontType),
            legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
            axis.text.x = element_text(size = xTickLblSize, family = fontType),
            axis.text.y = element_text(size = yTickLblSize, family = fontType, hjust = 0.5))

    if (ellipse){ # circles
      biplt <- biplt +
        stat_ellipse(aes(colour = Group, group = Group), type = "norm", level = ellipse_conf)
    }

    if (loadingPlot){ # superimpose loading plot
      biplt <- biplt +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_text(data = loadingValuePlot, aes(x = axis1, y = axis2), label = loadingValuePlot$lbl, colour = "gray30", size = loadingSize)

    }


    ggsave(filename = paste(objTitle,".PCA.biplot.pdf", sep = ""), plot = biplt,
           width = biplotWidth, height = biplotHeight, units = "mm",dpi = 600)

    grid.draw(biplt) # preview

  }

  return(assign(paste(objTitle, "_PCA", sep = ""), PCA, envir = .GlobalEnv)) # output PCA object

}


#' @title rbioFS_PCA.file
#'
#' @description A simple to use wrapper for PCA (Principal Component Analysis) and visualization. This is the version that loads \code{csv} file directly.
#' @param file Input file name. Format is \code{csv}.
#' @param ... Augument for \code{\link{rbioFS_PCA}}.
#' @return Outputs a PCA object, a boxplot (proportion of variance) and a biplot from PCA analysis. The format is \code{pdf}.
#' @details Make sure to arrange input data with first two columns for smaple ID and conditions, and the rest for features (e.g., genes).
#' @import ggplot2
#' @examples
#' \dontrun{
#' rbioFS_PCA.file(file = "data.csv", ellipse = TRUE, loadingPlot = TRUE, biplotWidth = 200, biplotHeight = 170)
#' }
#' @export
rbioFS_PCA.file <- function(file, ...){

  ## import file
  raw <- read.csv(file = file, header = TRUE, na.strings = c("NA", ""), stringsAsFactors = FALSE, check.names = FALSE)

  ## PCA
  rbioFS_PCA(input = raw, ...)
}
