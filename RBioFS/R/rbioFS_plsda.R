#' @title dummy
#'
#' @description dummification of catagorical variables
#' @param x Input vector. Make sure it is a factor.
#' @param scale Logical, whether to scale the data or not. Default is \code{TRUE}.
#' @return Outputs a matrix with dummified variables.
#' @details This function is needed to pre-process data when conducting plsda analysis. And the function can also be used for other purposes when needed.
#' @examples
#' \dontrun{
#' y <- dummy(y)
#' }
#' @export
dummy <- function (x, drop2nd = FALSE){  # integrate into the main function eventually
  if (!is.factor(x)){  # use this one for the arugments
    stop("'x' should be a factor")
  }
  y <- model.matrix(~x - 1)
  # below: remove unecessary array attributes
  colnames(y) <- gsub("^x", "", colnames(y))
  attributes(y)$assign <- NULL
  attributes(y)$contrasts <- NULL

  if (length(levels(x)) == 2 & drop2nd) { # if to only keep the factor binary
    y <- y[, 1]
  }
  return(y)
}



#' @title rbioFS_plsda_jackknife
#'
#' @description Jack-Knife procedure for the \code{PLS} models, e.g. \code{PLS-DA} or \code{PLS-R}.
#' @param object A \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ncomp Defaults is all the components the \code{mvr} object has.
#' @param use.mean Defaults is \code{FALSE}.
#' @param sig.p Alpha value for the jack-knife coffecient p values. Defaults is \code{0.05}.
#' @param plot.title Whether to display plot title on top of the plot. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param outlineCol The outline colour for the bar gars. Default is \code{"black"}.
#' @param greyScale To set the graph in grey scale. Default is \code{TRUE}.
#' @param errorbar Set the type of errorbar. Options are standard error of the mean (\code{"SEM"}, \code{"standard error"}, \code{"standard error of the mean"}), or standard deviation (\code{"SD"}, \code{"standard deviation"}), case insensitive. Default is \code{"SEM"}.
#' @param errorbarWidth Set the width for errorbar. Default is \code{0.2}.
#' @param errorbarLblSize Set the label size for the errorbar. Default is \code{6}.
#' @param fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param xLabel x axis label. Type with quotation marks. Default is \code{NULL}.
#' @param xLabelSize x axis label size. Default is \code{10}.
#' @param xTickLblSize Font size of x axis ticks. Default is \code{10}.
#' @param xTickItalic Set x axis tick font to italic. Default is \code{FALSE}.
#' @param xTickBold Set x axis tick font to bold. Default is \code{FALSE}.
#' @param xAngle The rotation angle (degrees) of the x axis marks. Default is \code{0} - horizontal.
#' @param xAlign The alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param rightsideY If to display the right side y-axis. Default is \code{TRUE}.
#' @param yLabel y axis label. Type with quotation marks. Default is \code{NULL}.
#' @param yLabelSize y axis label size. Default is \code{10}.
#' @param yTickLblSize Font size of y axis ticks. Default is \code{10}.
#' @param yTickItalic Set y axis tick font to italic. Default is \code{FALSE}.
#' @param yTickBold Set y axis tick font to bold. Default is \code{FALSE}.
#' @param legendSize Legend size. Default is \code{9}.
#' @param legendTtl Hide/Display legend title. If \code{TRUE} or \code{T}, the name of the first column of the raw data file will display as the legend title. Default is \code{FALSE}.
#' @param legendTtlSize Set when \code{legendTtl = TRUE}, font size of the legend title. Default is \code{9}.
#' @param plotWidth The width of the plot (unit: mm). Default is 170. Default will fit most of the cases.
#' @param plotHeight The height of the plot (unit: mm). Default is 150. Default will fit most of the cases.
#' @return Outputs two list objects to the environment, one for raw Jack-Knife results matrices, one for plot dataframe. Also the function also generates the pdf figure files to the working directory.
#' @details \code{use.mean = FALSE} is more main stream. Make sure to use cross validated and optimized component number for \code{ncomp}.
#' @importFrom reshape2 melt
#' @importFrom multcompView multcompLetters
#' @importFrom multcomp glht mcp
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @importFrom scales rescale_none
#' @import pls
#' @import ggplot2
#' @examples
#' \dontrun{
#' rbioFS_plsda_jackknife(object = new_model_optm, use.mean = FALSE,
#'                        sig.p = 0.05, plot = TRUE, plot.title = TRUE, plot.titleSize = 10,
#'                        outlineCol = "black", errorbar = "SEM", errorbarWidth = 0.2,
#'                        errorbarLblSize = 6, fontType = "sans", xLabel = "Features",
#'                        xLabelSize = 10, xTickLblSize = 10, xTickItalic = FALSE,
#'                        xTickBold = FALSE, xAngle = 90, xAlign = 1, rightsideY = TRUE,
#'                        yLabel = "Coefficients", yLabelSize = 10, yTickLblSize = 10,
#'                        yTickItalic = FALSE, yTickBold = FALSE, legendSize = 9,
#'                        legendTtl = FALSE, legendTtlSize = 9, plotWidth = 170,
#'                        plotHeight = 150)
#' }
#' @export
rbioFS_plsda_jackknife <- function(object, ncomp = object$ncomp, use.mean = FALSE, sig.p = 0.05,
                                   plot = TRUE, plot.title = FALSE, plot.titleSize = 10,
                                   outlineCol = "black",
                                   errorbar = "SEM", errorbarWidth = 0.2, errorbarLblSize = 6,
                                   fontType = "sans",
                                   xLabel = NULL, xLabelSize = 10, xTickLblSize = 10, xTickItalic = FALSE, xTickBold = FALSE, xAngle = 0, xAlign = 0.5,
                                   rightsideY = TRUE,
                                   yLabel = NULL, yLabelSize = 10, yTickLblSize = 10, yTickItalic = FALSE, yTickBold = FALSE,
                                   legendSize = 9, legendTtl = FALSE, legendTtlSize = 9,
                                   plotWidth = 170, plotHeight = 150){
  # compute sd, df, t-value, p-value for jackknife
  nresp <- dim(object$coefficients)[2]
  sdjack <- sqrt(var.jack(object, ncomp = ncomp, covariance = FALSE,
                          use.mean = use.mean))  # var.jack from pls pacakge
  sem <- sdjack / sqrt(length(object$validation$segments))
  B <- coef(object, ncomp = ncomp)
  df <- length(object$validation$segments) - 1
  tvals <- B/sdjack
  pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)

  # output
  B <- as.matrix(B[,,])
  tvals <- as.matrix(tvals[,,])
  pvals <- as.matrix(pvals[,,])
  sdjack <- as.matrix(sdjack[,,])
  sem <- as.matrix(sem[,,])
  out <- list(coefficients = B, sd = sdjack, sem = sem , df = df, tvalues = tvals, pvalues = pvals, ncomp = ncomp)

  if (plot){
    loclEnv <- environment()
    plot_list <- foreach(i = 1:dim(object$coefficients)[2]) %do% {
      tmp <- data.frame(features = dimnames(object$coefficients)[[1]], coefficients = out$coefficients[, i],
                        sd = out$sd[, i], sem = out$sem[, i], tvalues = out$tvalues[, i],
                        pvalues = out$pvalues[, i],
                        row.names = NULL, stringsAsFactors = FALSE)
      tmp$sig <- ifelse(tmp$pvalues < sig.p, "*", "")
      return(tmp)
    }
    names(plot_list) <- dimnames(object$coefficients)[[2]]

    for (i in 1:length(plot_list)){
      DfPlt <- plot_list[[i]]

      if (tolower(errorbar) %in% c("sem", "standard error", "standard error of the mean")){  # error bar
        err <- DfPlt$sem
      } else if (tolower(errorbar) %in% c("sd", "standard deviation")){
        err <- DfPlt$sd
      }

      ymax <- (max(DfPlt$coefficients[DfPlt$coefficients > 0] + err[DfPlt$coefficients > 0]) * 1.05) * 1.2
      ymin <- (min(DfPlt$coefficients[DfPlt$coefficients < 0] - err[DfPlt$coefficients < 0]) * 1.15) * 1.2
      y_axis_Mx <- max(abs(ymax), abs(ymin))
      y_axis_Mn <- max(abs(ymax), abs(ymin)) * sign(ymin)  # make sure the y-axes have the same abs value

      baseplt <- ggplot(data = DfPlt, aes(x = features, y = coefficients)) +
        geom_bar(position = "dodge", stat = "identity", color = outlineCol) +
        ggtitle(names(plot_list)[i]) +
        geom_errorbar(aes(ymin = coefficients - err, ymax = coefficients + err),
                      position = position_dodge(0.9), color = "black", width = errorbarWidth) +
        geom_text(aes(y = ifelse(sign(coefficients) > 0, (coefficients + err) * 1.05, (coefficients - err) * 1.15), label = sig),
                  position = position_dodge(width = 0.9), color = "black", size = 10) +
        scale_y_continuous(expand = c(0, 0), limits = c(y_axis_Mn, y_axis_Mx),
                           oob = rescale_none) +
        xlab(xLabel) +
        ylab(yLabel) +
        geom_hline(yintercept = 0) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", size = plot.titleSize, family = fontType),
              axis.title.x = element_text(face = "bold", size = xLabelSize, family = fontType),
              axis.title.y = element_text(face = "bold", size = xLabelSize, family = fontType),
              legend.position = "bottom",
              legend.text = element_text(size = legendSize),
              axis.text.x = element_text(size = xTickLblSize, family = fontType, angle = xAngle, hjust = xAlign),
              axis.text.y = element_text(size = yTickLblSize, family = fontType, hjust = 0.5))

      if (plot.title){
        baseplt <- baseplt + ggtitle(names(plot_list)[i])
      } else (
        baseplt <- baseplt + ggtitle(NULL)
      )

      if (xTickItalic & xTickBold){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "bold.italic"))
      } else if (xTickItalic & !xTickBold){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "italic"))
      } else if (xTickBold & !xTickItalic){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "bold"))
      }

      if (yTickItalic & yTickBold){
        baseplt <- baseplt +
          theme(axis.text.y  = element_text(face = "bold.italic"))
      } else if (yTickItalic & !yTickBold){
        baseplt <- baseplt +
          theme(axis.text.y = element_text(face = "italic"))
      } else if (yTickBold & !yTickItalic){
        baseplt <- baseplt +
          theme(axis.text.y = element_text(face = "bold"))
      }

      if (legendTtl == FALSE){
        baseplt <- baseplt + theme(legend.title = element_blank())
      } else {
        baseplt <- baseplt + theme(legend.title = element_text(size = legendTtlSize))
      }

      plt <- baseplt
      ## finalize the plot
      grid.newpage()
      if (rightsideY){ # add the right-side y axis
        # extract gtable
        pltgtb <- ggplot_gtable(ggplot_build(plt))
        # add the right side y axis
        Aa <- which(pltgtb$layout$name == "axis-l")
        pltgtb_a <- pltgtb$grobs[[Aa]]
        axs <- pltgtb_a$children[[2]]
        axs$widths <- rev(axs$widths)
        axs$grobs <- rev(axs$grobs)
        axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
        Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
        pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
        pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)
      } else { # no right side y-axis
        pltgtb <- plt
      }

      ## export the file and draw a preview
      cat(paste("Plot saved to file: ", deparse(substitute(object)), ".", names(plot_list)[i], ".", ".jackknife.pdf...", sep = "")) # initial message
      ggsave(filename = paste(deparse(substitute(object)), ".", names(plot_list)[i], ".", ".jackknife.pdf", sep = ""), plot = pltgtb,
             width = plotWidth, height = plotHeight, units = "mm",dpi = 600)
      cat("Done!\n") # final message
      grid.draw(pltgtb) # preview
    }
  }

  assign(paste(deparse(substitute(object)), "_jackknife_mtx_list", sep = ""), out, envir = .GlobalEnv)
  assign(paste(deparse(substitute(object)), "_jackknife_plot_list", sep = ""), plot_list, envir = .GlobalEnv)
}
