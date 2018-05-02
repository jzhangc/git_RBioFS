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


#' @title rbioFS_plsda
#'
#' @description PLS-DA modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param ncomp Number of components to be used for modelling.
#' @param method PLS-DA modelling method. Four PLSR algorithms are available: the kernel algorithm ("kernelpls"), the wide kernel algorithm ("widekernelpls"), SIMPLS ("simpls") and the classical orthogonal scores algorithm ("oscorespls"). Default is the popular \code{"simpls}.
#' @param scale Logical, whether to scale the data or not. Default is \code{TRUE}.
#' @param validation Cross validation methods. Options are "none", "CV" (fold), "LOO" (leave-one-out). Default is \code{"CV"}.
#' @param segments Set only when \code{validation = "CV}, the number of segement to be set. Default is \code{10}.
#' @param segments.type Method to set up the segments. Options are \code{"random", "consecutive", "interleaved"}.Default is \code{"random"}.
#' @param jackknife If to use jack-knife procedure. Default is \code{TRUE}.
#' @param ... Additional arguments for \code{mvr} function from \code{pls} pacakge.
#' @return Returns a PLS-DA model object, with class "mvr".
#' @details For sequencing data, the input x needs to be either tranformed with the function like \code{clr_ilr_transfo()} from \code{RBioArray} package, or normalized using methods like "TMM" or "RLE" implemented in \code{edgeR} pacakge.
#' @importFrom pls plsr
#' @examples
#' \dontrun{
#' y <- dummy(y)
#' }
#' @export
rbioFS_plsda <- function(x, y, ncomp, method = "simpls", scale = TRUE, validation = c("none", "CV", "LOO"),
                         segments = 10, segments.type = "random",
                         jackknife = TRUE, ...){
  ## check arguments
  if (!class(x) %in% c("matrix", "data.frame")){
    stop(cat("x has to be either a matrix or data frame."))
  }
  if (class(x) == "data.frame"){
    message(cat("data frame x converted to a matrix object."))
    x <- as.matrix(x)
  }
  if (!is.factor(y))stop("y has to be a factor object.")
  if (is.null(ncomp))stop("please set the ncomp number.")

  ## data processing
  # X
  cat(paste0("data centered with the scale option ", ifelse(scale, "\"ON\" ", "\"OFF\" "), "prior to modelling..."))
  centered_X <- center_scale(x, scale = scale)  # center data with the option of scaling
  cat("DONE!\n")
  X <- centered_X$centerX
  # Y
  Y <- dummy(y)
  # constract dataframe for modelling
  model_dfm <- data.frame(n = paste("row", 1:nrow(Y), sep = ""))
  model_dfm$Y <- Y  # y is the y matrix as a whole, not the content of y
  model_dfm$X <- X  # x is the x matrix as a whole, not the content of x

  ## modelling
  out_model <- plsr(Y ~ X, data = model_dfm, ncomp = ncomp,
                    method = method, scale = FALSE,
                    validation = validation, segments = segments, segments.type = segments.type, jackknife = TRUE, ...)
  out_model$centerX <- centered_X
  out_model$inputX <- x
  out_model$inputY <- y

  class(out_model) <- c("rbiomvr", "mvr")
  return(out_model)
}


#' @title rbioFS_plsda_tuplot
#'
#' @description T-U plot function for PLS-DA models.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param comps Integer vector. Components to plot. The index of the components are intergers. The vector length should be between 1 and the total number of components, inclusive. Can be Default is \code{c(1, 2)}.
#' @param multi_tuplot.ncol Set only when \code{length(comps) > 1}, number of columns in one figure page. Default is \code{length(comps)}.
#' @param multi_tuplot.nrow Set only when \code{length(comps) > 1}, number of rows in one figure page. Default is \code{1}.
#' @param rightsideY If to show the right side y-axis. Only applicble when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}. Note: the right side Y is ignored when \code{length(comps) > 1}
#' @param sampleLabel.type If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param sampleLabel.vector Set only when \code{sampleLabel.type} is not set to \code{"none"}, a character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param sampleLabel.padding Set only when \code{sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param tuplot.Title Scoreplot title. Default is \code{NULL}.
#' @param tuplot.SymbolSize Symbol size. Default is \code{2}.
#' @param tuplot..fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param tuplot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param tuplot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param tuplot.Width Scoreplot width. Default is \code{170}.
#' @param tuplot.Height Scoreplot height. Default is \code{150}.
#' @return Returns a pdf file for scoreplot.
#' @details The T-U plot shows the correlation betweem the decomposed x and y matrices for PLS-DA anlaysis. Such plot is useful to inspecting X-Y decomposition and outlier detection. Since PLS-DA is a classification modelling method, the well-correlated T-U plot will likely show a logistic regression-like plot, as opposed to a linear regressoin plot. The \code{sampleLabel} series arguments make it possible to show exact sample, something useful for outlier detection. The function supports plotting multiple components at the same time, i.e. multiple plots on one page. The right side y-axis is not applicable when plotting multiple components.
#' @import ggplot2
#' @import ggrepel
#' @import foreach
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioFS_plsda_tuplot(new_model, comps = c(1, 2, 3), scoreplot.ellipse = TRUE)
#' }
#' @export
rbioFS_plsda_tuplot <- function(object, comps = 1, multi_tuplot.ncol = length(comps), multi_tuplot.nrow = 1,
                                rightsideY = TRUE,
                                sampleLabel.type = "none", sampleLabel.vector = NULL, sampleLabel.padding = 0.5,
                                tuplot.SymbolSize = 5, tuplot.Title = NULL,
                                tuplot.fontType = "sans", tuplot.xTickLblSize = 10, tuplot.yTickLblSize = 10,
                                tuplot.Width = 170, tuplot.Height = 150){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.\n")
  if (length(comps) > object$ncomp) stop("comps length exceeded the maximum comp length.\n")
  if (!all(comps %in% seq(object$ncomp))) stop("comps contain non-existant comp.\n")
  if (!tolower(sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\". \n")
  if (rightsideY){
    if (length(comps) > 1){
      cat("right side y-axis ignored for multi-plot figure.\n")
    }
  }

  ## extract and construt t-u score plot dataframe
  tu_dfm <- data.frame(y = object$inputY)

  varpp_x <- 100 * object$Xvar / object$Xtotvar
  var_percentage_x <- varpp_x[paste0("Comp ", comps)] # extract the proportion of variance for the selected comps

  ## plot
  cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.tuplot.pdf...", sep = ""))  # initial message
  plt_list <- vector(mode = "list", length = length(comps))  # list
  plt_list[] <- foreach(i = comps) %do% {
    tu_dfm <- data.frame(t = object$scores[, i], u = object$Yscores[, i], tu_dfm)
    var_percentage_x <- varpp_x[paste0("Comp ", i)] # extract the proportion of variance for the selected comps
    comp_axis_lbl <- paste("(comp ", i, ", ", round(var_percentage_x, digits = 2), "%)", sep = "")
    lbl <- paste(c("t ", "u "), comp_axis_lbl, sep = "")
    plt <- ggplot(data = tu_dfm)  # base plot

    if (sampleLabel.type != "none"){
      if (is.null(sampleLabel.vector)){
        cat("sampleLabel.vector not provided. proceed without sampole labels.\n")
        plt <- plt + geom_point(size = tuplot.SymbolSize, aes(x = t, y = u, colour = y, shape = y))
      } else if (length(sampleLabel.vector) != nrow(tu_dfm)){
        cat("sampleLabel.vector not the same length as the number of samples. proceed without sampole labels.\n")
        plt <- plt + geom_point(size = tuplot.SymbolSize, aes(x = t, y = u, colour = y, shape = y))
      } else {
        tu_dfm$samplelabel <- as.character(sampleLabel.vector)
        if (tolower(sampleLabel.type) == "direct"){
          plt <- plt + geom_text(data = tu_dfm, aes(x = t, y = u, colour = y, label = samplelabel), size = tuplot.SymbolSize)
        } else if (tolower(sampleLabel.type) == "indirect") {
          plt <- plt + geom_point(size = tuplot.SymbolSize, aes(x = t, y = u, colour = y, shape = y)) +
            geom_text_repel(data = tu_dfm, aes(x = t, y = u, label = samplelabel), point.padding = unit(sampleLabel.padding, "lines"))
        }
      }
    } else {
      plt <- plt + geom_point(size = tuplot.SymbolSize, aes(x = t, y = u, colour = y, shape = y))
    }

    plt <- plt +
      ggtitle(tuplot.Title) +
      xlab(lbl[1]) +
      ylab(lbl[2]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = tuplot.fontType, hjust = 0.5),
            axis.title = element_text(face = "bold", family = tuplot.fontType),
            legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
            axis.text.x = element_text(size = tuplot.xTickLblSize, family = tuplot.fontType),
            axis.text.y = element_text(size = tuplot.yTickLblSize, family = tuplot.fontType, hjust = 0.5))

    if (rightsideY & length(comps) == 1){
      plt <- rightside_y(plt)
    }
    plt
  }
  names(plt_list) <- paste0("g", comps)

  if (length(comps) > 1) plt <- multi_plot_shared_legend(plt_list, ncol = multi_tuplot.ncol, nrow = multi_tuplot.nrow, position = "bottom")

  ## save
  grid.newpage()
  ggsave(filename = paste(deparse(substitute(object)),".plsda.tuplot.pdf", sep = ""), plot = plt,
         width = tuplot.Width, height = tuplot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  cat("Done!\n")

  ## return invinsible object
  invisible(plt)
}


#' @title rbioFS_plsda_scoreplot
#'
#' @description scoreplot function for PLS-DA models.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param y Set when object class is \code{mvr}. Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param comps Integer vector. Components to plot. The index of the components are intergers. The vector length should be between 1 and the total number of components, inclusive. Can be Default is \code{c(1, 2)}.
#' @param rightsideY If to show the right side y-axis. Only applicble when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}.
#' @param scoreplot.Title Scoreplot title. Default is \code{NULL}.
#' @param scoreplot.SymbolSize Symbol size. Default is \code{2}.
#' @param scoreplot.ellipse If to draw ellipses. Default is \code{FALSE}.
#' @param scoreplot.ellipse_conf The confidence value for the ellipses. Default is \code{0.95}.
#' @param scoreplot.mtx.densityplot If to display a density plot on the diagonal for the correlation scoreplot matrix. Default is \code{FALSE}.
#' @param scoreplot.mtx.stripLblSize The label font size for the correlation scoreplot matrix strips. Default is \code{10}.
#' @param scoreplot.xAngle The rotation angle (degrees) of the x axis marks. Default is \code{0} - horizontal.
#' @param scoreplot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param scoreplot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param scoreplot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param scoreplot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param scoreplot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param scoreplot.Width Scoreplot width. Default is \code{170}.
#' @param scoreplot.Height Scoreplot height. Default is \code{150}.
#' @return Returns a pdf file for scoreplot.
#' @details When \code{length(comps) == 1}, the function generates a scatter plot plotting sample vs score for the comp of interest. When \code{length(comps) == 2}, the function generates a scatter plot plotting the two comps of interest against each other. When \code{length(comps) > 2}, the function generates a multi-panel correlation scoreplot matrix for the comps of interest - might be slow if the there are many comps.
#' @import ggplot2
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @examples
#' \dontrun{
#' rbioFS_plsda_scoreplot(new_model, comps = c(1, 2, 3), scoreplot.ellipse = TRUE)
#' }
#' @export
rbioFS_plsda_scoreplot <- function(object, y = NULL, comps = c(1, 2),
                                   rightsideY = FALSE,
                                   scoreplot.Title = NULL,
                                   scoreplot.SymbolSize = 2,
                                   scoreplot.ellipse = FALSE, scoreplot.ellipse_conf = 0.95,
                                   scoreplot.mtx.densityplot = FALSE, scoreplot.mtx.stripLblSize = 10,
                                   scoreplot.xAngle = 0, scoreplot.xhAlign = 0.5, scoreplot.xvAlign = 0.5,
                                   scoreplot.fontType = "sans", scoreplot.xTickLblSize = 10, scoreplot.yTickLblSize = 10,
                                   scoreplot.Width = 170, scoreplot.Height = 150){
  ## check variables
  if (class(object) == "rbiomvr"){
    y <- object$inputY
  } else if (class(object) == "mvr"){
    if (is.null(y))stop("for mvr class objects, please provide response factor vector y.\n")
    if (!is.factor)stop("for mvr class objects, y has to be a factor vector.\n")
    y <- response
  } else stop("object needs to be \"rbiomvr\" or \"mvr\" class, e.g. created from RBioFS_plsda() function.\n")
  if (length(comps) > object$ncomp)stop("comps length exceeded the maximum comp length.\n")
  if (!all(comps %in% seq(object$ncomp)))stop("comps contain non-existant comp.\n")

  ## extract information
  score_x <- data.frame(object$scores[, comps, drop = FALSE], check.names = FALSE)
  score_x$group <- y
  varpp_x <- 100 * object$Xvar / object$Xtotvar
  boxdfm_x <- data.frame(comp_x = as.numeric(gsub("Comp", "", names(varpp_x))), varpp_x = varpp_x)
  var_percentage_x <- varpp_x[paste0("Comp ", comps)] # extract the proportion of variance for the selected PCs
  comp_axis_lbl <- paste("comp ", comps, " (", round(var_percentage_x, digits = 2), "%)", sep = "")

  ## scoreplot plotting
  if (length(comps) == 1){  # one component plot
    score_x$sample <- as.numeric(rownames(score_x))
    names(score_x)[1] <- "axis1"

    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggplot(score_x, aes(x = sample, y = axis1)) +
      geom_line(aes(colour = group, linetype = group)) +
      geom_point(aes(shape = group, colour = group), size = scoreplot.SymbolSize) + # plot the sample score scatter plot
      ggtitle(scoreplot.Title) +
      ylab(comp_axis_lbl[1]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = scoreplot.fontType, hjust = 0.5),
            axis.title = element_text(face = "bold", family = scoreplot.fontType),
            legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
            axis.text.x = element_text(size = scoreplot.xTickLblSize, family = scoreplot.fontType, angle = scoreplot.xAngle, hjust = scoreplot.xhAlign, vjust = scoreplot.xvAlign),
            axis.text.y = element_text(size = scoreplot.yTickLblSize, family = scoreplot.fontType, hjust = 0.5))

    grid.newpage()
    if (rightsideY){ # add the right-side y axis
      # extract gtable
      pltgtb <- rightside_y(scoreplt)
    } else { # no right side y-axis
      pltgtb <- scoreplt
    }

  } else if (length(comps) == 2){  # two components plot
    names(score_x)[1:2] <- c("axis1", "axis2")

    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggplot(score_x, aes(x = axis1, y = axis2)) +
      geom_point(aes(shape = group, colour = group), size = scoreplot.SymbolSize) + # plot the sample score scatter plot
      ggtitle(scoreplot.Title) +
      xlab(comp_axis_lbl[1]) +
      ylab(comp_axis_lbl[2]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = scoreplot.fontType, hjust = 0.5),
            axis.title = element_text(face = "bold", family = scoreplot.fontType),
            legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
            axis.text.x = element_text(size = scoreplot.xTickLblSize, family = scoreplot.fontType, angle = scoreplot.xAngle, hjust = scoreplot.xhAlign, vjust = scoreplot.xvAlign),
            axis.text.y = element_text(size = scoreplot.yTickLblSize, family = scoreplot.fontType, hjust = 0.5))
    if (scoreplot.ellipse){ # circles
      scoreplt <- scoreplt +
        stat_ellipse(aes(colour = group, group = group), type = "norm", level = scoreplot.ellipse_conf)
    }

    grid.newpage()
    if (rightsideY){ # add the right-side y axis
      pltgtb <- rightside_y(scoreplt)
    } else { # no right side y-axis
      pltgtb <- scoreplt
    }

  } else if (length(comps) > 2){  # over two components plot matrix
    if (rightsideY){
      cat("Right side y-axis ignored for comps more than 2...\n")
    }

    # custom functions for the paired scoreplot
    if (scoreplot.ellipse){  # ellipse
      ellipsefunc <- function(data = score_x, mapping, ellipse_conf = scoreplot.ellipse_conf, ...){
        ggplot(data = data, mapping = mapping) +
          geom_point(...) +
          stat_ellipse(aes(colour = group, group = group), type = "norm", level = ellipse_conf)
      }
    } else {
      ellipsefunc <- function(data = score_x, mapping, ellipse_conf = scoreplot.ellipse_conf, ...){
        ggplot(data = data, mapping = mapping) +
          geom_point(...)
      }
    }
    if (scoreplot.mtx.densityplot){  # diag densityplot
      densityfunc <- function(data = score_x, mapping, alpha = 0.1, ...){
        ggplot(data = data, mapping = mapping) +
          geom_density(alpha = alpha, aes(colour = group, linetype = group, ...))
      }
    } else {
      densityfunc <- function(data = score_x, mapping, alpha = 0.1){
        ggplot(data = data, mapping = mapping)
      }
    }

    # matrx scoreplot
    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggpairs(score_x, columns = comps, aes(colour = group, shape = group),
                        axisLabels = "show", columnLabels = comp_axis_lbl,
                        showStrips = NULL,
                        lower = list(continuous = ellipsefunc),
                        upper = list(continuous = ellipsefunc),
                        diag = list(continuous = densityfunc),
                        legend = 2)
    scoreplt <- scoreplt +
      ggtitle(scoreplot.Title) +
      theme(plot.title = element_text(face = "bold", family = scoreplot.fontType, hjust = 0.5),
            axis.title = element_text(face = "bold", family = scoreplot.fontType),
            strip.background = element_blank(),  # no strip background colour
            strip.text = element_text(face = "bold", size = scoreplot.mtx.stripLblSize),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            legend.position = "bottom", legend.title = element_blank(), legend.key = element_blank(),
            axis.text.x = element_text(size = scoreplot.xTickLblSize, family = scoreplot.fontType, angle = scoreplot.xAngle, hjust = scoreplot.xhAlign, vjust = scoreplot.xvAlign),
            axis.text.y = element_text(size = scoreplot.yTickLblSize, family = scoreplot.fontType))

    grid.newpage()
    pltgtb <- scoreplt
  }
  ggsave(filename = paste(deparse(substitute(object)),".plsda.scoreplot.pdf", sep = ""), plot = pltgtb,
         width = scoreplot.Width, height = scoreplot.Height, units = "mm",dpi = 600)
  cat("Done!\n") # final message
  grid.draw(pltgtb)

  ## return invinsible object
  invisible(pltgtb)
}


#' @title rbioFS_plsda_jackknife
#'
#' @description Jack-Knife procedure for the \code{PLS} models, e.g. \code{PLS-DA} or \code{PLS-R}.
#' @param object A \code{mvr} or \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
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
#' @param xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
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
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
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
#'                        xTickBold = FALSE, xAngle = 90, xhAlign = 1, xvAligh = 0.2, rightsideY = TRUE,
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
                                   xLabel = NULL, xLabelSize = 10, xTickLblSize = 10, xTickItalic = FALSE, xTickBold = FALSE, xAngle = 0,
                                   xhAlign = 0.5, xvAlign = 0.5,
                                   rightsideY = TRUE,
                                   yLabel = NULL, yLabelSize = 10, yTickLblSize = 10, yTickItalic = FALSE, yTickBold = FALSE,
                                   legendSize = 9, legendTtl = FALSE, legendTtlSize = 9,
                                   plotWidth = 170, plotHeight = 150){
  # check arguments
  if (!any(class(object) %in% c("mvr", "rbiomvr")))stop("object has to be a mvr or rbiomvr class.\n")

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
              axis.text.x = element_text(size = xTickLblSize, family = fontType, angle = xAngle,
                                         hjust = xhAlign, vjust = xvAlign),
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
        pltgtb <- rightside_y(plt)
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
