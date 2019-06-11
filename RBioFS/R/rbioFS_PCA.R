#' @title rbioFS_PCA
#'
#' @description A simple to use wrapper for PCA (Principal Component Analysis) and visualization
#' @param input Input data, data frame.
#' @param sampleIDVar Sample variable name. It's a character string.
#' @param groupIDVar Group variable name. It's a character string.
#' @param scaleData If to scale the data when performing PCA. Default is \code{TRUE}.
#' @param centerData If to center the data when performing PCA. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{prcomp} function.
#' @param boxplot If to plot a boxplot for variance distribution. Default is \code{TRUE}.
#' @param boxplot.Title The boxplot title. Default is \code{NULL}.
#' @param boxplot.Width The boxplot width. Default is \code{170}.
#' @param boxplot.Height The boxplot height. Default is \code{150}.
#' @param biplot If to plot a scoreplot, with the option of superimpose a loading plot, i.e. biplot. Default is \code{TRUE}.
#' @param biplot.comps Integer or vector of integers. Index number(s) for principal component(s) to plot. Default is \code{c(1, 2)}.
#' @param biplot.Title The biplot title. Default is \code{NULL}.
#' @param biplot.sampleLabel.type  If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param biplot.sampleLabelSize Only set when \code{biplot.sampleLabel.type} is not \code{"none"}, The size of the sample label. Default is \code{2}.
#' @param biplot.sampleLabel.padding Set only when \code{biplot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param biplot.SymbolSize The symbol size for the scatter plot portion of the biplot. Default is \code{2}.
#' @param biplot.ellipse If to draw ellipses. Default is \code{FALSE}.
#' @param biplot.ellipse_conf The confidence value for the ellipses. Default is \code{0.93}.
#' @param biplot.xAngle The rotation angle (degrees) of the biplot x axis marks. Default is \code{0} - horizontal.
#' @param biplot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param biplot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param biplot.loadingplot If to superimpose loading plot. Default is \code{FALSE}.
#' @param biplot.loadingplot.textsize The font size of the loading plot labels. Default is \code{3}.
#' @param biplot.mtx.densityplot If to display a density plot on the diagonal for the correlation scoreplot matrix. Default is \code{FALSE}.
#' @param biplot.mtx.stripLblSize The label font size for the correlation scoreplot matrix strips. Default is \code{10}.
#' @param biplot.Width The biplot width. Default is \code{170}.
#' @param biplot.Height The biplot height. Default is \code{150}.
#' @param rightsideY If to show the right side y-axis for both boxplot and biplot. For biplot, only applicble when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}.
#' @param fontType Font for the figure texts. Default is \code{"sans"}.
#' @param xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Outputs a PCA object, a boxplot (proportion of variance) and a biplot from PCA analysis. The format is \code{pdf}.
#' @details Make sure to arrange input data with first two columns for smaple ID and conditions, and the rest for features (e.g., genes).
#' @import ggplot2
#' @import ggrepel
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @examples
#' \dontrun{
#' rbioFS_PCA(input = pcaDfm, sampleIDVar = "conditions", groupIDVar = "samples",
#' biplot.ellipse = TRUE, biplot.loadingplot = TRUE, biplot.Width = 200, biplot.Height = 170)
#' }
#' @export
rbioFS_PCA <- function(input = NULL, sampleIDVar = NULL, groupIDVar = NULL, scaleData = TRUE, centerData = TRUE, ...,
                       boxplot = TRUE,
                       boxplot.Title = NULL,
                       boxplot.Width = 170, boxplot.Height = 150,
                       biplot = TRUE, biplot.comps = c(1:2),
                       biplot.Title = NULL,
                       biplot.sampleLabel.type = "none", biplot.sampleLabelSize = 2,
                       biplot.sampleLabel.padding = 0.5,
                       biplot.SymbolSize = 2,
                       biplot.ellipse = FALSE, biplot.ellipse_conf = 0.95,
                       biplot.xAngle = 0, biplot.xhAlign = 0.5, biplot.xvAlign = 0.5,
                       biplot.loadingplot  = FALSE, biplot.loadingplot.textsize = 3,
                       biplot.mtx.densityplot = FALSE, biplot.mtx.stripLblSize = 10,
                       biplot.Width = 170, biplot.Height = 150,
                       rightsideY = FALSE,
                       fontType = "sans", xTickLblSize = 10, yTickLblSize = 10,
                       verbose = TRUE){
  ## check the argument
  if (!all(c(sampleIDVar, groupIDVar) %in% names(input))) stop("sampleIDvar and/or groupIDvar not found in the input dataframe.")
  biplot.sampleLabel.type <- match.arg(tolower(biplot.sampleLabel.type), c("none", "direct", "indirect"))

  ## set up input
  x <- input[, !names(input) %in% c(sampleIDVar, groupIDVar)]

  ## argument check
  if (biplot & length(biplot.comps) > ncol(x))stop("biplot.comps length exceeded the maximum PC length.")
  if (biplot & !all(biplot.comps %in% seq(ncol(x))))stop("biplot.comps contain non-existant PC.")

  ## PCA
  PCA <- prcomp(x, scale. = scaleData, center = centerData, ...)
  varpp_x <- 100 * summary(PCA)$importance[2, ] # extract and calcualte the proportion of variance
  boxdfm_x <- data.frame(PC = as.numeric(gsub("PC", "", names(varpp_x))), varpp = varpp_x)

  ## Boxplot
  if (boxplot){
    if (verbose) cat(paste("Boxplot being saved to file: ", deparse(substitute(input)), ".pca.boxplot.pdf...", sep = ""))  # initial message
    # grid.newpage()
    boxplt <- ggplot(data = boxdfm_x, aes(x = PC, y = varpp_x, group = 1)) +
      geom_bar(position = "dodge", stat = "identity", color = "black", fill = "gray66") +
      scale_x_continuous() +
      scale_y_continuous(expand = c(0, 0),
                         limits = c(0, with(boxdfm_x, ceiling(max(varpp_x))) * 1.1)) +
      ggtitle(boxplot.Title) +
      xlab("PC") +
      ylab("Proportion of variance (%)") +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = fontType, hjust = 0.5),
            axis.title = element_text(face = "bold", family = fontType),
            legend.position = "bottom",legend.title = element_blank(),legend.key = element_blank(),
            axis.text.x = element_text(size = xTickLblSize, family = fontType),
            axis.text.y = element_text(size = yTickLblSize, family = fontType, hjust = 0.5))

    if (rightsideY){ # add the right-side y axis
      boxplt <- rightside_y(boxplt)  # port back to boxplt object
    }

    ggsave(filename = paste(deparse(substitute(input)),".pca.boxplot.pdf", sep = ""), plot = boxplt,
           width = boxplot.Width, height = boxplot.Height, units = "mm",dpi = 600)
    if (verbose) cat("Done!\n")
    grid.draw(boxplt) # preview
  }

  ## biplot = scoreplot + loadingplot
  # prepare for scatter plot values (i.e. sample score)
  score_x <- data.frame(PCA$x[, biplot.comps, drop = FALSE], check.names = FALSE) # extract rotated sample scores
  score_x$group <- factor(input[, groupIDVar], levels = unique(input[, groupIDVar]))
  score_x$sample.label <- input[, sampleIDVar]
  var_percentage_x <- varpp_x[paste0("PC", biplot.comps)] # extract the proportion of variance for the selected PCs
  pc_axis_lbl <- paste("PC ", biplot.comps, " (", round(var_percentage_x, digits = 2), "%)", sep = "")

  # plotting
  if (biplot){
    if (length(biplot.comps) == 1){
      if (biplot.loadingplot){
        cat("loading plot is not applicable when length(biplot.comps) == 1. Proceed without one.\n")
      }
      if (biplot.ellipse){
        cat("ellipse is not applicable when length(biplot.comps) == 1. Proceed without one.\n")
      }

      score_x$sample <- as.numeric(rownames(score_x))
      names(score_x)[1] <- "axis1"

      if (verbose) cat(paste("Single PC biplot being saved to file: ", deparse(substitute(input)), ".pca.biplot.pdf...", sep = ""))  # initial message
      biplt <- ggplot(score_x, aes(x = sample, y = axis1)) +
        geom_line(aes(colour = group, linetype = group))

      # sample labels
      if (biplot.sampleLabel.type == "direct"){
        biplt <- biplt +
          geom_text(aes(colour = group, label = sample.label), size = biplot.SymbolSize)
      } else if (biplot.sampleLabel.type == "indirect"){
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = biplot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(biplot.sampleLabel.padding, "lines"),
                          show.legend = FALSE, size = biplot.sampleLabelSize)
      } else {
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = biplot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
      }

      biplt <- biplt +
        ggtitle(biplot.Title) +
        ylab(pc_axis_lbl[1]) +
        theme_bw() +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = xTickLblSize, family = fontType, angle = biplot.xAngle, hjust = biplot.xhAlign, vjust = biplot.xhAlign),
              axis.text.y = element_text(size = yTickLblSize, family = fontType, hjust = 0.5))

      # grid.newpage()
      if (rightsideY){ # add the right-side y axis
        biplt <- RBioplot::rightside_y(biplt)
      }
    } else if (length(biplot.comps) == 2){
      names(score_x)[1:2] <- c("axis1", "axis2")

      if (verbose) cat(paste("Biplot being saved to file: ", deparse(substitute(input)), ".pca.biplot.pdf...", sep = ""))  # initial message
      biplt <- ggplot(score_x, aes(x = axis1, y = axis2))

      # sample labels
      if (biplot.sampleLabel.type == "direct"){
        biplt <- biplt +
          geom_text(aes(colour = group, label = sample.label), size = biplot.SymbolSize)
      } else if (biplot.sampleLabel.type == "indirect"){
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = biplot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
        geom_text_repel(aes(label = sample.label), point.padding = unit(biplot.sampleLabel.padding, "lines"),
                        show.legend = FALSE, size = biplot.sampleLabelSize)
      } else {
        biplt <- biplt +
          geom_point(aes(shape = group, colour = group), size = biplot.SymbolSize) +
          scale_shape_manual(values=1:nlevels(score_x$group))
      }

      biplt <- biplt +
        ggtitle(biplot.Title) +
        xlab(pc_axis_lbl[1]) +
        ylab(pc_axis_lbl[2]) +
        theme_bw() +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = xTickLblSize, family = fontType, angle = biplot.xAngle, hjust = biplot.xhAlign, vjust = biplot.xvAlign),
              axis.text.y = element_text(size = yTickLblSize, family = fontType, hjust = 0.5))

      if (biplot.ellipse){ # circles
        biplt <- biplt +
          stat_ellipse(aes(colour = group, group = group), type = "norm", level = biplot.ellipse_conf)
      }

      if (biplot.loadingplot){ # superimpose loading plot
        # prepare for loading plot values (i.e. loading value for variables)
        loadingValue <- data.frame(PCA$rotation[, biplot.comps], check.names = FALSE) # extract loading/rotation/eigenvectors for variables
        names(loadingValue) <- c("axis1", "axis2") # give a generic variable name for the ratation dataframe
        loadingScale <- max(max(abs(score_x$axis1)) / max(abs(loadingValue$axis1)),
                            max(abs(score_x$axis2)) / max(abs(loadingValue$axis2))) * 0.85 # determine scaling for loading values
        loadingValuePlot <- loadingValue * loadingScale
        loadingValuePlot$lbl <- rownames(loadingValuePlot)

        # plot
        biplt <- biplt +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_text(data = loadingValuePlot, aes(x = axis1, y = axis2), label = loadingValuePlot$lbl, colour = "gray30", size = biplot.loadingplot.textsize)
      }

      # grid.newpage()
      if (rightsideY){ # add the right-side y axis
        biplt <- RBioplot::rightside_y(biplt)
      }

    } else if (length(biplot.comps) > 2) {
      if (rightsideY){ # add the right-side y axis
        cat("Rightside y-axis is not applicable when length(biplot.comps) > 2. Proceed without one.\n")
      }

      if (biplot.loadingplot) {
        cat("For now, loadingplot is unavailable when using more then two PCs. \n")
      }

      # custom functions for the paired scoreplot
      ellipsefunc <- function(data = score_x, mapping, label.method = biplot.sampleLabel.type,
                              ellipse = biplot.ellipse, ellipse_conf = biplot.ellipse_conf,
                              ...){
        g <- ggplot(data = data, mapping = mapping)

        if (label.method == "direct"){
          g <- g + geom_text(aes(colour = group, label = sample.label), ...)
        } else if (label.method == "indirect"){
          g <- g +  geom_point(...) +
            scale_shape_manual(values = 1:nlevels(data$group)) +
            geom_text_repel(aes(label = sample.label), point.padding = unit(biplot.sampleLabel.padding, "lines"),
                            show.legend = FALSE, size = biplot.sampleLabelSize)
        } else {
          g <- g + geom_point(...) +
            scale_shape_manual(values = 1:nlevels(data$group))
        }
        if (ellipse){
          g <- g +
            stat_ellipse(aes(colour = group, group = group), type = "norm", level = ellipse_conf)
        }
        return(g)
      }

      densityfunc <- function(data = score_x, mapping, alpha = 0.1, densityplot = biplot.mtx.densityplot, ...){
        g <- ggplot(data = data, mapping = mapping)
        if (densityplot){
          g <- g + geom_density(alpha = alpha, aes(colour = group, linetype = group, ...))
        }
        return(g)
      }

      # matrix scoreplot
      if (verbose) cat(paste("Biplot matrix being saved to file: ", deparse(substitute(input)),".pca.biplot.pdf...", sep = ""))  # initial message
      biplt <- ggpairs(score_x, columns = biplot.comps, aes(colour = group, shape = group),
                       axisLabels = "show", columnLabels = pc_axis_lbl,
                       showStrips = NULL,
                       lower = list(continuous = ellipsefunc),
                       upper = list(continuous = ellipsefunc),
                       diag = list(continuous = densityfunc),
                       legend = 2)
      biplt <- biplt +
        ggtitle(biplot.Title) +
        theme(plot.title = element_text(face = "bold", family = fontType, hjust = 0.5),
              axis.title = element_text(face = "bold", family = fontType),
              strip.background = element_blank(),  # no strip background colour
              strip.text = element_text(face = "bold", size = biplot.mtx.stripLblSize),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
              axis.text.x = element_text(size = xTickLblSize, family = fontType, angle = biplot.xAngle, hjust = biplot.xhAlign, vjust = biplot.xvAlign),
              axis.text.y = element_text(size = yTickLblSize, family = fontType))

      # grid.newpage()
    }
    ggsave(filename = paste(deparse(substitute(input)),".pca.biplot.pdf", sep = ""), plot = biplt,
           width = biplot.Width, height = biplot.Height, units = "mm",dpi = 600)
    if (verbose) cat("Done!\n") # final message
    grid.draw(biplt)
  }

  assign(paste(deparse(substitute(input)), "_pca", sep = ""), PCA, envir = .GlobalEnv)
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
