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
#' @param ncomp Number of components to be used for modelling. Default is \code{length(unique(y)) - 1}.
#' @param method PLS-DA modelling method. Four PLSR algorithms are available: the kernel algorithm ("kernelpls"), the wide kernel algorithm ("widekernelpls"), SIMPLS ("simpls") and the classical orthogonal scores algorithm (also known as NIPALS) ("oscorespls"). Default is the popular \code{"simpls}.
#' @param scale Logical, whether to scale the data or not. Default is \code{TRUE}.
#' @param validation Cross validation methods. Options are "none", "CV" (fold), "LOO" (leave-one-out). Default is \code{"CV"}.
#' @param segments Set only when \code{validation = "CV}, the number of segement to be set. Default is \code{10}.
#' @param segments.type Method to set up the segments. Options are \code{"random", "consecutive", "interleaved"}.Default is \code{"random"}.
#' @param jackknife If to use jack-knife procedure. Default is \code{TRUE}.
#' @param ... Additional arguments for \code{mvr} function from \code{pls} pacakge.
#' @return Returns a PLS-DA model object, with class "mvr".
#' @details For \code{ncomp} value, the default (full model) compares feature number with \code{class number - 1}, instead of observation number. For sequencing data, the input x needs to be either tranformed with the function like \code{clr_ilr_transfo()} from \code{RBioArray} package, or normalized using methods like "TMM" or "RLE" implemented in \code{edgeR} pacakge.
#' @importFrom pls plsr
#' @examples
#' \dontrun{
#' rbioFS_plsda(x, y, ncomp = 20)
#' }
#' @export
rbioFS_plsda <- function(x, y, ncomp = length(unique(y)) - 1, method = "simpls", scale = TRUE, validation = c("none", "CV", "LOO"),
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
  out_model$dummy_y <- Y
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
#' @param multi_plot.ncol Set only when \code{length(comps) > 1}, number of columns on one figure page. Default is \code{length(comps)}.
#' @param multi_plot.nrow Set only when \code{length(comps) > 1}, number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.rightsideY If to show the right side y-axis. Only applicble when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}. Note: the right side Y is ignored when \code{length(comps) > 1}
#' @param plot.sampleLabel.type If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param plot.sampleLabel.vector Set only when \code{plot.sampleLabel.type} is not set to \code{"none"}, a character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param plot.sampleLabelSize Only set when \code{plot.sampleLabel.type} is not \code{"none"}. The size of the sample label. Default is \code{2}.
#' @param plot.sampleLabel.padding Set only when \code{plot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.Title Scoreplot title. Default is \code{NULL}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
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
#' rbioFS_plsda_tuplot(new_model, comps = c(1, 2, 3))
#' }
#' @export
rbioFS_plsda_tuplot <- function(object, comps = 1, multi_plot.ncol = length(comps), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                                plot.rightsideY = TRUE,
                                plot.sampleLabel.type = "none", plot.sampleLabel.vector = NULL,
                                plot.sampleLabelSize = 2, plot.sampleLabel.padding = 0.5,
                                plot.SymbolSize = 5, plot.Title = NULL,
                                plot.fontType = "sans",
                                plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                plot.legendSize = 9,
                                plot.Width = 170, plot.Height = 150){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.\n")
  if (length(comps) > object$ncomp) stop("comps length exceeded the maximum comp length.\n")
  if (!all(comps %in% seq(object$ncomp))) stop("comps contain non-existant comp.\n")
  if (!tolower(plot.sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\". \n")
  if (plot.rightsideY){
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

    if (plot.sampleLabel.type != "none"){
      if (is.null(plot.sampleLabel.vector)){
        cat("sampleLabel.vector not provided. proceed without sampole labels.\n")
        plt <- plt + geom_point(size = plot.SymbolSize, aes(x = t, y = u, colour = y, shape = y))
      } else if (length(plot.sampleLabel.vector) != nrow(tu_dfm)){
        cat("sampleLabel.vector not the same length as the number of samples. proceed without sampole labels.\n")
        plt <- plt + geom_point(size = plot.SymbolSize, aes(x = t, y = u, colour = y, shape = y))
      } else {
        tu_dfm$samplelabel <- as.character(plot.sampleLabel.vector)
        if (tolower(plot.sampleLabel.type) == "direct"){
          plt <- plt + geom_text(data = tu_dfm, aes(x = t, y = u, colour = y, label = samplelabel), size = plot.SymbolSize)
        } else if (tolower(plot.sampleLabel.type) == "indirect") {
          plt <- plt + geom_point(size = plot.SymbolSize, aes(x = t, y = u, colour = y, shape = y)) +
            geom_text_repel(data = tu_dfm, aes(x = t, y = u, label = samplelabel),
                            point.padding = unit(plot.sampleLabel.padding, "lines"), size = plot.sampleLabelSize,
                            show.legend = FALSE)
        }
      }
    } else {
      plt <- plt + geom_point(size = plot.SymbolSize, aes(x = t, y = u, colour = y, shape = y))
    }

    plt <- plt +
      ggtitle(plot.Title) +
      xlab(lbl[1]) +
      ylab(lbl[2]) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
            axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
            axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
            legend.position = "bottom", legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
            legend.key = element_blank(),
            axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
            axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

    if (plot.rightsideY & length(comps) == 1){
      plt <- rightside_y(plt)
    }
    plt
  }
  names(plt_list) <- paste0("g", comps)

  if (length(comps) > 1) {
    if (multi_plot.ncol * multi_plot.nrow < length(comps)){
      stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. make sure they match the number of components used.")
    } else {
      plt <- RBioplot::multi_plot_shared_legend(plt_list, ncol = multi_plot.ncol, nrow = multi_plot.nrow, position = multi_plot.legend.pos)
    }
  }

  ## save
  grid.newpage()
  ggsave(filename = paste(deparse(substitute(object)),".plsda.tuplot.pdf", sep = ""), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  cat("Done!\n")

  ## return invinsible object
  invisible(plt)
}


#' @title rbioFS_plsda_q2r2()
#'
#' @description q2-r2 (i.e. Q^2 and R^2 scores) caluclation and plot for plsda models
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param intercept Wether to include intercept term, i.e. comps = 0. Default is \code{TRUE}.
#' @param q2r2plot If to generate a q2r2 plot. Default is \code{TRUE}.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is the number of responding classes, i.e. unique y classes.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}. Note: doesn't seem to be necessasry as PLS-DA always has at least two y classes.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Components"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"R2 & Q2"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @return Prints the selected number of components for each y class. Returns RMSEP values for each y class to the environment, as well as a pdf file for the RMSEP plot if \code{rmsepplot = TRUE}.
#' @details A vertical line indicating the number of component with minimum q2-r2 distance.
#' @import ggplot2
#' @import foreach
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @importFrom pls mvrValstats R2
#' @examples
#' \dontrun{
#' rbioFS_plsda_q2r2(object = new_model, multi_plot.ncol = 2, multi_plot.nrow = 2, intercept = TRUE)
#' }
#' @export
rbioFS_plsda_q2r2 <- function(object, intercept = TRUE, q2r2plot = TRUE,
                              plot.display.Title = TRUE,
                              multi_plot.ncol = length(dimnames(object$coefficients)[[2]]), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                              plot.rightsideY = TRUE, plot.fontType = "sans",
                              plot.SymbolSize = 2,
                              plot.xLabel = "Components", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                              plot.yLabel = "R2 & Q2", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                              plot.legendSize = 9,
                              plot.Width = 170, plot.Height = 150){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.\n")
  if (is.null(object$validation) || is.null(object$validation$coefficients)) stop("'object' was not fit with jackknifing enabled")  # from pls pacakge

  ## construct q2r2 dataframes
  # calculate q2 and r2
  sst <- pls::mvrValstats(object, estimate = "train")$SST
  press_mtx <- object$validation$PRESS
  press_mtx <- cbind(object$validation$PRESS0, object$validation$PRESS)
  colnames(press_mtx)[1] <- "0 comps"
  q2 <- 1 - press_mtx / c(sst)  # calculate q2. q2 = 1 - PRESS/SST
  r2 <- pls::R2(object, intercept = TRUE, estimate = "train")$val  # R2() calulates R2. FYI: r2 = 1 - SSR/SST. Use estimate = "train" to use SSR
  dimnames(r2)[[3]][1] <- "0 comps"

  # construct and output a multi-dfm list
  q2r2_dfm_list <- vector(mode = "list", length = length(rownames(press_mtx)))
  q2r2_dfm_list[] <- foreach(i = 1:length(rownames(press_mtx))) %do% {
    q2_dfm <- data.frame(components = as.integer(gsub(" comps", "", colnames(q2))), Q2 = q2[i, ],
                         stringsAsFactors = FALSE, row.names = NULL)
    r2_dfm <- data.frame(components = as.integer(gsub(" comps", "", dimnames(r2)[[3]])), R2 = r2[, i, ],
                         row.names = NULL)
    q2r2_dfm <- merge(q2_dfm, r2_dfm, by = "components")
    if (!intercept){
      q2r2_dfm <- q2r2_dfm[-1, ]
    }
    q2r2_dfm
  }
  names(q2r2_dfm_list) <- rownames(press_mtx)
  assign(paste(deparse(substitute(object)), "_q2r2_list", sep = ""), q2r2_dfm_list, envir = .GlobalEnv)

  ## plot
  if (q2r2plot){
    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.q2r2plot.pdf...", sep = ""))  # initial message
    plt_list <- vector(mode = "list", length = length(q2r2_dfm_list))  # list
    plt_list[] <- foreach(i = 1:length(q2r2_dfm_list)) %do% {
      q2r2_plot_dfm <- reshape2::melt(q2r2_dfm_list[[i]], id = "components")
      plt <- ggplot(data = q2r2_plot_dfm, aes(x = components, y = value, color = variable)) +
        geom_line(aes(linetype = variable)) +
        geom_point(size = plot.SymbolSize) +
        scale_x_continuous(breaks = c(ifelse(intercept, 0, 1), q2r2_dfm_list[[i]][which.min(abs(q2r2_dfm_list[[i]][, 2] - q2r2_dfm_list[[i]][, 3])), 1],
                                      floor(median(seq(ifelse(intercept, 0, 1), nrow(q2r2_dfm_list[[i]])))),
                                      q2r2_dfm_list[[i]][nrow(q2r2_dfm_list[[i]]), 1])) +
        scale_y_continuous() +
        ggtitle(ifelse(plot.display.Title, names(q2r2_dfm_list)[i], NULL)) +
        geom_vline(xintercept = q2r2_dfm_list[[i]][which.min(abs(q2r2_dfm_list[[i]][, 2] - q2r2_dfm_list[[i]][, 3])), 1], linetype = "dashed", colour = "red") +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              legend.position = "bottom", legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
              legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

      if (plot.rightsideY & length(names(q2r2_dfm_list)) == 1){
        plt <- RBioplot::rightside_y(plt)
      }
      plt
    }
    names(plt_list) <- names(q2r2_dfm_list)

    if (length(names(q2r2_dfm_list)) > 1) {
      if (multi_plot.ncol * multi_plot.nrow < length(q2r2_dfm_list)){
        stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. make sure they match the number of response groups.")
      } else {
        plt <- RBioplot::multi_plot_shared_legend(plt_list, ncol = multi_plot.ncol, nrow = multi_plot.nrow, position = multi_plot.legend.pos)
      }
    }

    # save
    grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.q2r2plot.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    cat("Done!\n")
  }
}


#' @title randomiz.test
#'
#' @description Randomization function from \code{pls} pacakge. For \code{\link{rbioFS_plsda_ncomp_select}} function.
randomiz.test <- function(residualsNew, residualsReference, nperm){
  d <- residualsNew^2 - residualsReference^2
  md <- mean(d)
  N <- length(d)

  signs <- round(matrix(runif(N * nperm), N, nperm)) * 2 - 1
  dsigns <- d * signs
  mdsigns <- colMeans(dsigns)

  count <- sum(mdsigns >= md) ## equality will never occur, right?
  (count + 0.5) / (nperm + 1)
}


#' @title rbioFS_plsda_ncomp_select()
#'
#' @description Optimal number of components selection for PLS-DA model, with RMSEP plot funcitonality. Selection methods are modified based on \code{selectNcomp()} from \code{pls} pacakge.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ... Additional argument for \code{RMSEP} function from \code{pls} package.
#' @param ncomp.selection.method Optimal numbers of components selection method. Options are \code{"min"}, \code{"1sd"}, and \code{"randomization"}. Default is \code{"1sd"}.
#' @param randomization.nperm Set only when \code{ncomp.selection.method = "randomization"}, number of permutations. Default is \code{999}.
#' @param randomization.alpha Set only when \code{ncomp.selection.method = "randomization"}, alpha for the p values used during "randomization" selection. Default is \code{0.05}.
#' @param rmsepplot If to generate a RMSEP plot. Default is \code{TRUE}.
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}. Note: doesn't seem to be necessasry as PLS-DA always has at least two y classes.
#' @param plot.optm.ncomp.line If to display the vertical line indicting the optimal number of components. Default is \code{TRUE}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is the number of responding classes, i.e. y.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Components"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"RMSEP"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @return Prints the selected number of components for each y class. Returns RMSEP values for each y class to the environment, as well as a pdf file for the RMSEP plot if \code{rmsepplot = TRUE}.
#' @details The RMSEP figure shows both CV estimates and adjusted CV estimates, which is CV estimiates corrected for bias. Three methods are used for components number selection: \code{"min"} simply chooses the number of components to reach te minimum RMSEP; \code{"1sd"} chooses the number of components when its RMSEP first reaches minimum as well as within one standard deviation; For "randomization", see the help file for \code{selectNcomp()} function from  \code{pls} pacakge.
#' @import ggplot2
#' @import foreach
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @importFrom pls RMSEP
#' @examples
#' \dontrun{
#' rbioFS_plsda_ncomp_select(new_model,  multi_plot.ncol = 2, multi_plot.nrow = 2, plot.optm.ncomp.line = T,
#'          ncomp.selection.method = "randomization", randomization.nperm = 999, randomization.alpha = 0.05)
#' }
#' @export
rbioFS_plsda_ncomp_select <- function(object, ...,
                                      ncomp.selection.method = "1sd", randomization.nperm = 999, randomization.alpha = 0.05,
                                      rmsepplot = TRUE,
                                      plot.rightsideY = TRUE,
                                      plot.optm.ncomp.line = TRUE,
                                      multi_plot.ncol = length(dimnames(object$coefficients)[[2]]), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                                      plot.display.Title = TRUE,
                                      plot.SymbolSize = 2,
                                      plot.fontType = "sans",
                                      plot.xLabel = "Components", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                      plot.yLabel = "RMSEP", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                      plot.legendSize = 9,
                                      plot.Width = 170, plot.Height = 150){
  ## check arguments
  if (!any(class(object) %in% c("mvr", "rbiomvr"))) stop("object has to be a mvr or rbiomvr class.\n")
  if (!"validation" %in% names(object)) stop("PLS-DA model has to include Cross-Validation.\n")
  if (!tolower(ncomp.selection.method) %in% c("min", "1sd", "randomization")) stop("ncomp.selection.method needs to be \"min\", \"1sd\", or \"randomization\" exactly.")

  ## calcuate RMSEP
  rmsep <- RMSEP(object, ...)

  ## ncomp selection and plot
  # prepare dataframe for ncomp selection and plotting
  rmsep_dfm_list <- vector(mode = "list", length = dim(rmsep$val)[2])
  plt_list <- vector(mode = "list", length = dim(rmsep$val)[2])
  rmsep_dfm_list[] <- foreach(i = 1:dim(rmsep$val)[2]) %do% {
    dfm <- as.data.frame(t(rmsep$val[,i,]))
    dfm <- data.frame(comps = seq(0, dim(rmsep$val)[3] - 1), dfm, row.names = NULL)
    dfm <- melt(dfm, id = "comps")
  }
  names(rmsep_dfm_list) <- dimnames(rmsep$val)[[2]]  # dimnames(tstrmsep$val)[2] is a list object. use [] to get the characters in string class.

  # ncomp selection
  ncompsel_mtx <- foreach(i = 1:length(rmsep_dfm_list), .combine = "rbind") %do% {
    dfm_cv <- rmsep_dfm_list[[i]][rmsep_dfm_list[[i]]$variable == "CV", ]
    dfm_adjcv <- rmsep_dfm_list[[i]][rmsep_dfm_list[[i]]$variable == "adjCV", ]

    if (ncomp.selection.method == "min"){
      cv_x <- dfm_cv$comps[which.min(dfm_cv$value)]
      adjcv_x <- dfm_adjcv$comps[which.min(dfm_adjcv$value)]
    } else {  # modified from pls::selectNcomp()
      origResponse <- object$dummy_y[, i]
      allresids <- cbind(origResponse - (sum(origResponse) - origResponse) / (length(origResponse) - 1),
                         object$validation$pred[,i,] - origResponse)  # CV prediction for all the comps used
      residsds <- apply(allresids, 2, sd) / sqrt(nrow(allresids))

      if (ncomp.selection.method == "1sd"){
        dfm_cv$SD <- residsds
        min_rmsep_cv_idx <- which.min(dfm_cv$value)  # index for the minimum mean oob feature group
        sd_min_cv <- dfm_cv$SD[min_rmsep_cv_idx]  # oob SD for the feature group above
        cv_x <- min(with(dfm_cv, which(value <= (value[min_rmsep_cv_idx] + sd_min_cv))))  # 1sd minimum selection

        dfm_adjcv$SD <- residsds
        min_rmsep_adjcv_idx <- which.min(dfm_adjcv$value)  # index for the minimum mean oob feature group
        sd_min_adjcv <- dfm_adjcv$SD[min_rmsep_adjcv_idx]  # oob SD for the feature group above
        adjcv_x <- min(with(dfm_adjcv, which(value <= (value[min_rmsep_adjcv_idx] + sd_min_adjcv))))  # 1sd minimum selection
      } else if (ncomp.selection.method == "randomization"){
        absBest_cv <- which.min(dfm_cv$value)
        pvals_cv <- sapply(seq_len(absBest_cv - 1), function(m) randomiz.test(allresids[, m], allresids[, absBest_cv], nperm = randomization.nperm))
        idx_cv <- which(pvals_cv > randomization.alpha)
        cv_x <- min(c(idx_cv, absBest_cv)) - 1

        absBest_adjcv <- which.min(dfm_adjcv$value)
        pvals_adjcv <- sapply(seq_len(absBest_adjcv - 1), function(m) randomiz.test(allresids[, m], allresids[, absBest_adjcv], nperm = randomization.nperm))
        idx_adjcv <- which(pvals_cv > randomization.alpha)
        adjcv_x <- min(c(idx_adjcv, absBest_adjcv)) - 1
      }
    }
    c(cv_x, adjcv_x)
  }
  rownames(ncompsel_mtx) <- names(rmsep_dfm_list)
  colnames(ncompsel_mtx) <- c("CV", "adjCV")
  cat(paste0("Optimial ncomp (selection method: ", ncomp.selection.method, ") : \n", sep = ""))
  print(ncompsel_mtx)

  # plotting
  if (rmsepplot){
    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.rmsepplot.pdf...", sep = ""))  # initial message
    plt_list[] <- foreach(i = 1:length(rmsep_dfm_list)) %do% {
      plt <- ggplot(data = rmsep_dfm_list[[i]], aes(x = comps, y = value, colour = variable)) +
        geom_line(aes(linetype = variable)) +
        geom_point(aes(shape = variable), size = plot.SymbolSize) +
        ggtitle(ifelse(plot.display.Title, names(rmsep_dfm_list)[i], NULL)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              legend.position = "bottom", legend.key = element_blank(),
              legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

      if (plot.optm.ncomp.line){
        plt <- plt +
          geom_vline(xintercept = ncompsel_mtx[i, ], linetype = "dashed", colour = "red")
      }

      if (plot.rightsideY & length(rmsep_dfm_list) == 1){
        plt <- rightside_y(plt)
      }
      plt
    }
    names(plt_list) <- paste0("g_", names(rmsep_dfm_list))

    if (length(rmsep_dfm_list) > 1) {
      if (multi_plot.ncol * multi_plot.nrow < length(rmsep_dfm_list)){
        stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. make sure they match the number of response groups.")
      } else {
        plt <- RBioplot::multi_plot_shared_legend(plt_list, ncol = multi_plot.ncol, nrow = multi_plot.nrow, position = multi_plot.legend.pos)
      }
    }

    ## save
    grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.rmsepplot.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    cat("Done!\n")
  }

  ## return RMSEP values
  assign(paste(deparse(substitute(object)), "_plsda_rmsep_list", sep = ""), rmsep_dfm_list, envir = .GlobalEnv)
}


#' @title rbioFS_plsda_scoreplot
#'
#' @description scoreplot function for PLS-DA models.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param y Set when object class is \code{mvr}. Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param comps Integer vector. Components to plot. The index of the components are intergers. The vector length should be between 1 and the total number of components, inclusive. Can be Default is \code{c(1, 2)}.
#' @param plot.rightsideY If to show the right side y-axis. Only applicble when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}.
#' @param plot.Title Scoreplot title. Default is \code{NULL}.
#' @param plot.sampleLabel.type If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param plot.sampleLabel.vector Set only when \code{plot.sampleLabel.type} is not set to \code{"none"}, a character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param plot.sampleLabelSize Only set when \code{plot.sampleLabel.type} is not \code{"none"}. The size of the sample label. Default is \code{2}.
#' @param plot.sampleLabel.padding Set only when \code{plot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.ellipse If to draw ellipses. Default is \code{FALSE}.
#' @param plot.ellipse_conf The confidence value for the ellipses. Default is \code{0.95}.
#' @param plot.mtx.densityplot If to display a density plot on the diagonal for the correlation scoreplot matrix. Default is \code{FALSE}.
#' @param plot.mtx.stripLblSize The label font size for the correlation scoreplot matrix strips. Default is \code{10}.
#' @param plot.xAngle The rotation angle (degrees) of the x-axis marks. Default is \code{0} - horizontal.
#' @param plot.xhAlign The horizontal alignment type of the x-axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.xvAlign The vertical alignment type of the x-axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @return Returns a pdf file for scoreplot.
#' @details When \code{length(comps) == 1}, the function generates a scatter plot plotting sample vs score for the comp of interest. When \code{length(comps) == 2}, the function generates a scatter plot plotting the two comps of interest against each other. When \code{length(comps) > 2}, the function generates a multi-panel correlation scoreplot matrix for the comps of interest - might be slow if the there are many comps.
#' @import ggplot2
#' @import ggrepel
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @examples
#' \dontrun{
#' rbioFS_plsda_scoreplot(new_model, comps = c(1, 2, 3), plot.ellipse = TRUE)
#' }
#' @export
rbioFS_plsda_scoreplot <- function(object, y = NULL, comps = c(1, 2),
                                   plot.rightsideY = FALSE,
                                   plot.Title = NULL,
                                   plot.sampleLabel.type = "none", plot.sampleLabel.vector = NULL,
                                   plot.sampleLabelSize = 2, plot.sampleLabel.padding = 0.5,
                                   plot.SymbolSize = 2,
                                   plot.ellipse = FALSE, plot.ellipse_conf = 0.95,
                                   plot.mtx.densityplot = FALSE, plot.mtx.stripLblSize = 10,
                                   plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                                   plot.fontType = "sans",
                                   plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                   plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                   plot.legendSize = 9,
                                   plot.Width = 170, plot.Height = 150){
  ## check variables
  if (any(class(object) == "rbiomvr")){
    y <- object$inputY
  } else if (!any(class(object) == "rbiomvr") & any(class(object) == "mvr")){
    if (is.null(y)) stop("for mvr class objects, please provide response factor vector y.\n")
    if (!is.factor(y)) stop("for mvr class objects, y has to be a factor vector.\n")
    y <- y
  } else stop("object needs to be \"rbiomvr\" or \"mvr\" class, e.g. created from RBioFS_plsda() function.\n")
  if (length(comps) > object$ncomp)stop("comps length exceeded the maximum comp length.\n")
  if (!all(comps %in% seq(object$ncomp)))stop("comps contain non-existant comp.\n")
  if (!tolower(plot.sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\". \n")
  if (tolower(plot.sampleLabel.type != "none")){
    if (!is.null(plot.sampleLabel.vector) & length(plot.sampleLabel.vector) != nrow(object$inputX)) {
      stop("plot.sampleLabel.vector has to be the same length as number of samples in the training set from object, i.e. nrow(object$inputX). \n")}
  }

  ## extract information
  score_x <- data.frame(object$scores[, comps, drop = FALSE], check.names = FALSE)
  score_x$group <- y
  if (is.null(plot.sampleLabel.vector)){  # sample label
    score_x$sample.label <- 1:nrow(object$inputX)
  } else {
    score_x$sample.label <- plot.sampleLabel.vector
  }
  varpp_x <- 100 * object$Xvar / object$Xtotvar
  boxdfm_x <- data.frame(comp_x = as.numeric(gsub("Comp", "", names(varpp_x))), varpp_x = varpp_x)
  var_percentage_x <- varpp_x[paste0("Comp ", comps)] # extract the proportion of variance for the selected PCs
  comp_axis_lbl <- paste("comp ", comps, " (", round(var_percentage_x, digits = 2), "%)", sep = "")

  ## scoreplot plotting
  if (plot.sampleLabel.type != "none" & is.null(plot.sampleLabel.vector)){
    cat("plot.sampleLabel.vector not provided. Proceed with row numbers as sampole labels.\n")
  }
  if (length(comps) == 1){  # one component plot
    score_x$sample <- as.numeric(rownames(score_x))
    names(score_x)[1] <- "axis1"

    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggplot(score_x, aes(x = sample, y = axis1)) +
      geom_line(aes(colour = group, linetype = group))

    if (plot.sampleLabel.type != "none"){  # labels
      if (tolower(plot.sampleLabel.type) == "direct"){
        scoreplt <- scoreplt + geom_text(aes(label = sample.label, colour = group), size = plot.SymbolSize)
      } else if (tolower(plot.sampleLabel.type) == "indirect") {
        scoreplt <- scoreplt + geom_point(alpha = 0.6, size = plot.SymbolSize, aes(colour = group, shape = group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          size = plot.sampleLabelSize, show.legend = FALSE)
      }
    } else {
      scoreplt <- scoreplt + geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) # plot the sample score scatter plot
    }

    scoreplt <- scoreplt +
      ggtitle(plot.Title) +
      ylab(comp_axis_lbl[1]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
            axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
            legend.position = "bottom", legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
            legend.key = element_blank(),
            axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xvAlign),
            axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

    grid.newpage()
    if (plot.rightsideY){ # add the right-side y axis
      # extract gtable
      pltgtb <- rightside_y(scoreplt)
    } else { # no right side y-axis
      pltgtb <- scoreplt
    }

  } else if (length(comps) == 2){  # two components plot
    names(score_x)[1:2] <- c("axis1", "axis2")

    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggplot(score_x, aes(x = axis1, y = axis2))

    if (plot.sampleLabel.type != "none"){  # labels
      if (tolower(plot.sampleLabel.type) == "direct"){
        scoreplt <- scoreplt + geom_text(aes(label = sample.label, colour = group), size = plot.SymbolSize)
      } else if (tolower(plot.sampleLabel.type) == "indirect") {
        scoreplt <- scoreplt + geom_point(alpha = 0.6, size = plot.SymbolSize, aes(colour = group, shape = group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          size = plot.sampleLabelSize, show.legend = FALSE)
      }
    } else {
      scoreplt <- scoreplt + geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) # plot the sample score scatter plot
    }

    scoreplt <- scoreplt +
      ggtitle(plot.Title) +
      xlab(comp_axis_lbl[1]) +
      ylab(comp_axis_lbl[2]) +
      theme_bw() +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
            axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
            axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
            legend.position = "bottom", legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
            legend.key = element_blank(),
            axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xvAlign),
            axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
    if (plot.ellipse){ # circles
      scoreplt <- scoreplt +
        stat_ellipse(aes(colour = group, group = group), type = "norm", level = plot.ellipse_conf)
    }

    grid.newpage()
    if (plot.rightsideY){ # add the right-side y axis
      pltgtb <- rightside_y(scoreplt)
    } else { # no right side y-axis
      pltgtb <- scoreplt
    }

  } else if (length(comps) > 2){  # over two components plot matrix
    if (plot.rightsideY){
      cat("Right side y-axis ignored for comps more than 2...\n")
    }

    # custom functions for the paired scoreplot
    ellipsefunc <- function(data = score_x, mapping, label.method = plot.sampleLabel.type,
                            ellipse = plot.ellipse, ellipse_conf = plot.ellipse_conf,
                            ...){
      g <- ggplot(data = data, mapping = mapping)

      if (label.method == "direct"){
        g <- g + geom_text(aes(colour = group, label = sample.label), ...)
      } else if (label.method == "indirect"){
        g <- g +  geom_point(...) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          size = plot.sampleLabelSize, show.legend = FALSE)
      } else {
        g <- g + geom_point(...)
      }
      if (ellipse){
        g <- g +
          stat_ellipse(aes(colour = group, group = group), type = "norm", level = ellipse_conf)
      }
      return(g)
    }

    densityfunc <- function(data = score_x, mapping, alpha = 0.1, densityplot = plot.mtx.densityplot, ...){
      g <- ggplot(data = data, mapping = mapping)
      if (densityplot){
        g <- g + geom_density(alpha = alpha, aes(colour = group, linetype = group, ...))
      }
      return(g)
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
      ggtitle(plot.Title) +
      theme(plot.title = element_text(face = "bold", family = plot.fontType, hjust = 0.5),
            axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
            axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
            strip.background = element_blank(),  # no strip background colour
            strip.text = element_text(face = "bold", size = plot.mtx.stripLblSize),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            legend.position = "bottom", legend.text = element_text(size = plot.legendSize), legend.title = element_blank(),
            legend.key = element_blank(),
            axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle, hjust = plot.xhAlign, vjust = plot.xvAlign),
            axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType))

    grid.newpage()
    pltgtb <- scoreplt
  }
  ggsave(filename = paste(deparse(substitute(object)),".plsda.scoreplot.pdf", sep = ""), plot = pltgtb,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
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
#' @param plot.outlineCol The outline colour for the bar gars. Default is \code{"black"}.
#' @param plot.errorbar Set the type of errorbar. Options are standard error of the mean (\code{"SEM"}, \code{"standard error"}, \code{"standard error of the mean"}), or standard deviation (\code{"SD"}, \code{"standard deviation"}), case insensitive. Default is \code{"SEM"}.
#' @param plot.errorbarWidth Set the width for errorbar. Default is \code{0.2}.
#' @param plot.errorbarLblSize Set the label size for the errorbar. Default is \code{6}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Features"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize Font size of x-axis ticks. Default is \code{10}.
#' @param plot.xTickItalic Set x-axis tick font to italic. Default is \code{FALSE}.
#' @param plot.xTickBold Set x-axis tick font to bold. Default is \code{FALSE}.
#' @param plot.xAngle The rotation angle (degrees) of the x axis marks. Default is \code{0} - horizontal.
#' @param plot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.rightsideY If to display the right side y-axis. Default is \code{TRUE}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"Coefficients"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Font size of y axis ticks. Default is \code{10}.
#' @param plot.yTickItalic Set y-axis tick font to italic. Default is \code{FALSE}.
#' @param plot.yTickBold Set y-axis tick font to bold. Default is \code{FALSE}.
#' @param plot.Width The width of the plot (unit: mm). Default is 170. Default will fit most of the cases.
#' @param plot.Height The height of the plot (unit: mm). Default is 150. Default will fit most of the cases.
#' @return Outputs a jacknife summary list object to the environment. The function also generates the pdf figure files to the working directory.
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
#'                        plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
#'                        plot.errorbarLblSize = 6, plot.fontType = "sans", plot.xLabel = "Features",
#'                        plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
#'                        plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 1, plot.xvAligh = 0.2, plot.rightsideY = TRUE,
#'                        plot.yLabel = "Coefficients", plot.yLabelSize = 10, plot.yTickLblSize = 10,
#'                        plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
#'                        plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
#'                        plot.Height = 150)
#' }
#' @export
rbioFS_plsda_jackknife <- function(object, ncomp = object$ncomp, use.mean = FALSE, sig.p = 0.05,
                                   plot = TRUE,
                                   plot.title = FALSE, plot.titleSize = 10,
                                   plot.outlineCol = "black",
                                   plot.errorbar = "SEM", plot.errorbarWidth = 0.2, plot.errorbarLblSize = 6,
                                   plot.fontType = "sans",
                                   plot.xLabel = NULL, plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
                                   plot.xTickBold = FALSE, plot.xAngle = 0,
                                   plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                                   plot.rightsideY = TRUE,
                                   plot.yLabel = NULL, plot.yLabelSize = 10, plot.yTickLblSize = 10, plot.yTickItalic = FALSE, plot.yTickBold = FALSE,
                                   plot.Width = 170, plot.Height = 150){
  # check arguments
  if (!any(class(object) %in% c("mvr", "rbiomvr"))) stop("object has to be a mvr or rbiomvr class.\n")

  # compute sd, df, t-value, p-value for jackknife
  nresp <- dim(object$coefficients)[2]
  sdjack <- sqrt(pls::var.jack(object, ncomp = ncomp, covariance = FALSE,
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

  plot_list <- vector(mode = "list", length = dim(object$coefficients)[2])
  plot_list[] <- foreach(i = 1:dim(object$coefficients)[2]) %do% {
    tmp <- data.frame(features = dimnames(object$coefficients)[[1]], coefficients = out$coefficients[, i],
                        sd = out$sd[, i], sem = out$sem[, i], tvalues = out$tvalues[, i],
                        pvalues = out$pvalues[, i],
                        row.names = NULL, stringsAsFactors = FALSE)
    tmp$sig <- ifelse(tmp$pvalues < sig.p, "*", "")
    return(tmp)
  }
  names(plot_list) <- dimnames(object$coefficients)[[2]]

  if (plot){
    loclEnv <- environment()
    if (plot.xTickLblSize == 0) cat("Due to plot.xTickLblSize = 0, x-axis ticks are hidden.\n")
    for (i in 1:length(plot_list)){
      cat(paste("Plot saved to file: ", deparse(substitute(object)), ".", names(plot_list)[i], ".jackknife.pdf...", sep = "")) # initial message
      DfPlt <- plot_list[[i]]

      if (tolower(plot.errorbar) %in% c("sem", "standard error", "standard error of the mean")){  # error bar
        err <- DfPlt$sem
      } else if (tolower(plot.errorbar) %in% c("sd", "standard deviation")){
        err <- DfPlt$sd
      }

      # use tryCatch in the case of non of the coef are > or < 0.
      ymax <- tryCatch((max(DfPlt$coefficients[DfPlt$coefficients > 0] + err[DfPlt$coefficients > 0]) * 1.05),
                       warning = function(err) return(0))
      ymin <- tryCatch((min(DfPlt$coefficients[DfPlt$coefficients < 0] - err[DfPlt$coefficients < 0]) * 1.15),
                       warning = function(err) return(0))
#      y_axis_Mx <- ifelse(ymax == 0, 0, max(abs(ymax), abs(ymin)))
#      y_axis_Mn <- ifelse(ymin == 0, 0, max(abs(ymax), abs(ymin)) * sign(ymin)) # make sure the y-axes have the same abs value
      y_axis_Mx <- ifelse(ymax == 0, 0, abs(ymax))
      y_axis_Mn <- ifelse(ymin == 0, 0, abs(ymin) * sign(ymin))

      baseplt <- ggplot(data = DfPlt, aes(x = features, y = coefficients)) +
        geom_bar(position = "dodge", stat = "identity", color = plot.outlineCol) +
        ggtitle(names(plot_list)[i]) +
        geom_errorbar(aes(ymin = coefficients - err, ymax = coefficients + err),
                      position = position_dodge(0.9), color = "black", width = plot.errorbarWidth) +
        geom_text(aes(y = ifelse(sign(coefficients) > 0, (coefficients + err) * 1.05, (coefficients - err) * 1.15), label = sig),
                  position = position_dodge(width = 0.9), color = "black", size = plot.errorbarLblSize) +
        scale_x_discrete(expand = c(0.01, 0.01)) +
        scale_y_continuous(expand = c(0.01, 0.01), limits = c(y_axis_Mn, y_axis_Mx),
                           oob = rescale_none) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        geom_hline(yintercept = 0) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle,
                                         hjust = plot.xhAlign, vjust = plot.xvAlign),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5),
              axis.ticks.x = if(plot.xTickLblSize == 0) element_blank())

      if (plot.title){
        baseplt <- baseplt + ggtitle(names(plot_list)[i])
      } else (
        baseplt <- baseplt + ggtitle(NULL)
      )

      if (plot.xTickItalic & plot.xTickBold){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "bold.italic"))
      } else if (plot.xTickItalic & !plot.xTickBold){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "italic"))
      } else if (plot.xTickBold & !plot.xTickItalic){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "bold"))
      }

      if (plot.yTickItalic & plot.yTickBold){
        baseplt <- baseplt +
          theme(axis.text.y  = element_text(face = "bold.italic"))
      } else if (plot.yTickItalic & !plot.yTickBold){
        baseplt <- baseplt +
          theme(axis.text.y = element_text(face = "italic"))
      } else if (plot.yTickBold & !plot.yTickItalic){
        baseplt <- baseplt +
          theme(axis.text.y = element_text(face = "bold"))
      }

      plt <- baseplt
      ## finalize the plot
      grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        pltgtb <- rightside_y(plt)
      } else { # no right side y-axis
        pltgtb <- plt
      }

      ## export the file and draw a preview
      ggsave(filename = paste(deparse(substitute(object)), ".", names(plot_list)[i], ".jackknife.pdf", sep = ""), plot = pltgtb,
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      cat("Done!\n") # final message
      grid.draw(pltgtb) # preview
    }
  }
  cat(paste0("jackknife is built on ", ncomp, " components.\n"))
  assign(paste(deparse(substitute(object)), "_plsda_jackknife_summary_list", sep = ""), plot_list, envir = .GlobalEnv)
}


#' @title rbioFS_plsda_VIP
#'
#' @description VIP, or variable importance in projection, calcualtion and plotting for plsda models. This is another FS method, and can be used independently.
#' @param object A \code{mvr} or \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
#' @param vip.alpha Alpha value (threshold) for VIP values. Any VIP above this is considered important. Defaults is \code{1}.
#' @param plot.title Whether to display plot title on top of the plot. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.sig.line Wether to display a horizontal line indicating the VIP threshold. Default is \code{TRUE}.
#' @param plot.outlineCol The outline colour for the bar gars. Default is \code{"black"}.
#' @param plot.errorbar Set the type of errorbar. Only applicable if the object model is built using more than one component. Options are standard error of the mean (\code{"SEM"}, \code{"standard error"}, \code{"standard error of the mean"}), or standard deviation (\code{"SD"}, \code{"standard deviation"}), case insensitive. Default is \code{"SEM"}.
#' @param plot.errorbarWidth Set the width for errorbar. Default is \code{0.2}.
#' @param plot.errorbarLblSize Set the label size for the errorbar. Default is \code{6}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Features"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize Font size of x axis ticks. Default is \code{10}.
#' @param plot.xTickItalic Set X-axis tick font to italic. Default is \code{FALSE}.
#' @param plot.xTickBold Set X-axis tick font to bold. Default is \code{FALSE}.
#' @param plot.xAngle The rotation angle (degrees) of the x axis marks. Default is \code{0} - horizontal.
#' @param plot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.rightsideY If to display the right side y-axis. Default is \code{TRUE}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"VIP"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Font size of y axis ticks. Default is \code{10}.
#' @param plot.yTickItalic Set Y-axis tick font to italic. Default is \code{FALSE}.
#' @param plot.yTickBold Set Y-axis tick font to bold. Default is \code{FALSE}.
#' @param plot.Width The width of the plot (unit: mm). Default is 170. Default will fit most of the cases.
#' @param plot.Height The height of the plot (unit: mm). Default is 150. Default will fit most of the cases.
#' @return Outputs a list objects to the environment with the VIP raw values, VIP summary and the ncomp value. Also the function also generates the pdf figure files to the working directory.
#' @details Only works when the plsda model is fitted with the orthorgonal score algorithm, or NIPALS. Such model can be built using \code{\link{rbioFS_plsda}} with \code{method = "oscorespls"}. For each feature, the boxplot is the mean of the VIP values from all the components, hence with errorbars. However, if the model is fitted using only one component, the function will automatically adjust. The VIP threshold of 1 is the most commonly accepted value. However it is also acceptable to set according to the data and objectives of the study.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @importFrom scales rescale_none
#' @examples
#' \dontrun{
#' rbioFS_plsda_VIP(object = new_model_optm,
#'                        vip.alpha = 0.05,
#'                        plot = TRUE, plot.title = TRUE, plot.titleSize = 10,
#'                        plot.sig.line = TRUE,
#'                        plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
#'                        plot.errorbarLblSize = 6, plot.fontType = "sans", plot.xLabel = "Features",
#'                        plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
#'                        plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 1, plot.xvAligh = 0.2, plot.rightsideY = TRUE,
#'                        plot.yLabel = "Coefficients", plot.yLabelSize = 10, plot.yTickLblSize = 10,
#'                        plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
#'                        plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
#'                        plot.Height = 150)
#' }
#' @export
rbioFS_plsda_VIP <- function(object, vip.alpha = 1,
                    plot = TRUE, plot.title = TRUE, plot.titleSize = 10,
                    plot.sig.line = TRUE,
                    plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
                    plot.errorbarLblSize = 6, plot.fontType = "sans",
                    plot.xLabel = "Features", plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
                    plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 0.95, plot.xvAlign = 0.5,
                    plot.rightsideY = TRUE,
                    plot.yLabel = "VIP", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                    plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
                    plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
                    plot.Height = 150){
  ## argument check
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.\n")
  if (object$method != "oscorespls") stop("Object needs to fit using oscorespls (i.e. NIPALS) algorithm. Please re-fit using RBioFS_plsda with method = \"oscorespls\".") # only oscorespls algorithm is applicable to VIP

  ## VIP calculation
  score <- object$scores
  lodw <- object$loading.weights
  Wnorm2 <- colSums(lodw^2)

  vip_raw_list <- vector(mode = "list", length = dim(object$Yloadings)[1])
  vip_raw_list[] <- foreach(i = 1:dim(object$Yloadings)[1]) %do% {
    ylod <- object$Yloadings[i, ] # one y variable at a time
    ss <- c(ylod)^2 * colSums(score^2)
    ssw <- sweep(lodw^2, 2, ss / Wnorm2, "*")
    vip <- sqrt(nrow(ssw) * apply(ssw, 1, cumsum) / cumsum(ss))
    return(vip)
  }
  names(vip_raw_list) <- dimnames(object$Yloadings)[[1]]

  vip_list <- vector(mode = "list", length = dim(object$Yloadings)[1])
  vip_list[] <- foreach(i = 1:dim(object$Yloadings)[1]) %do% {
    if (object$ncomp > 1){
      vip_mean <- foreach(m = vip_raw_list[[i]], .combine = "c") %do% mean(m)
      vip_sd <- foreach(n = vip_raw_list[[i]], .combine = "c") %do% sd(n)
      vip_sem <- foreach(o = vip_sd, .combine = "c") %do% sqrt(o / nrow(vip_raw_list[[i]]))
      vip_plot_dfm <- data.frame(features = colnames(vip_raw_list[[i]]), VIP = vip_mean, sd = vip_sd, sem = vip_sem)
    } else { # no sem or sd is needed in ncomp == 1
      vip_plot_dfm <- data.frame(features = names(vip), VIP = vip_raw_list[[i]])
    }
    vip_plot_dfm <- vip_plot_dfm[order(vip_plot_dfm$VIP, decreasing = TRUE), ]
    vip_plot_dfm$features <- factor(vip_plot_dfm$features, levels = unique(vip_plot_dfm$features))
    return(vip_plot_dfm)
  }
  names(vip_list) <- dimnames(object$Yloadings)[[1]]

  ## plot
  if (plot){
   if (plot.xTickLblSize == 0) cat("Due to plot.xTickLblSize = 0, x-axis ticks are hidden.\n")

    for (j in 1:dim(object$Yloadings)[1]){
      cat(paste("Plot saved to file: ", deparse(substitute(object)), ".", names(vip_list)[j], ".vip.pdf...", sep = "")) # initial message

      if (object$ncomp > 1){  #  detect if to use errorbar
        if (tolower(plot.errorbar) %in% c("sem", "standard error", "standard error of the mean")){  # error bar
          err <- vip_list[[j]]$sem
        } else if (tolower(plot.errorbar) %in% c("sd", "standard deviation")){
          err <- vip_list[[j]]$sd
        }
        y_axis_Mx <- max((vip_list[[j]]$VIP + err) * 1.05) * 1.2

      } else {
        y_axis_Mx <- max((vip_list[[j]]$VIP) * 1.05) * 1.2
      }
      y_axis_Mn <- 0

      baseplt <- ggplot(data = vip_list[[j]], aes(x = features, y = VIP)) +
        geom_bar(position = "dodge", stat = "identity", color = plot.outlineCol) +
        ggtitle(names(vip_list)[j]) +
        scale_x_discrete(expand = c(0.01, 0.01)) +
        scale_y_continuous(expand = c(0, 0), limits = c(y_axis_Mn, y_axis_Mx),
                           oob = rescale_none) +
        xlab(plot.xLabel) +
        geom_hline(yintercept = 0) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              legend.position = "bottom",
              legend.text = element_text(size = plot.legendSize),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle,
                                         hjust = plot.xhAlign, vjust = plot.xvAlign),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5),
              axis.ticks.x = if(plot.xTickLblSize == 0) element_blank())

      if (plot.title){
        baseplt <- baseplt + ggtitle(names(vip_list)[j])
      } else (
        baseplt <- baseplt + ggtitle(NULL)
      )

      if (plot.sig.line){
        baseplt <- baseplt +
          geom_hline(yintercept = vip.alpha, linetype = "dashed", colour = "red")
      }

      if (object$ncomp > 1){
        baseplt <- baseplt +
          geom_errorbar(aes(ymin = VIP - err, ymax = VIP + err),
                      position = position_dodge(0.9), color = "black", width = plot.errorbarWidth) +
          ylab(paste0(plot.yLabel, " (", object$ncomp, " comps)"))
      } else {
        baseplt <- baseplt +
          ylab(paste0(plot.yLabel, " (", object$ncomp, " comp)"))
      }

      if (plot.xTickItalic & plot.xTickBold){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "bold.italic"))
      } else if (plot.xTickItalic & !plot.xTickBold){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "italic"))
      } else if (plot.xTickBold & !plot.xTickItalic){
        baseplt <- baseplt +
          theme(axis.text.x = element_text(face = "bold"))
      }

      if (plot.yTickItalic & plot.yTickBold){
        baseplt <- baseplt +
          theme(axis.text.y  = element_text(face = "bold.italic"))
      } else if (plot.yTickItalic & !plot.yTickBold){
        baseplt <- baseplt +
          theme(axis.text.y = element_text(face = "italic"))
      } else if (plot.yTickBold & !plot.yTickItalic){
        baseplt <- baseplt +
          theme(axis.text.y = element_text(face = "bold"))
      }

      plt <- baseplt
      ## finalize the plot
      grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        pltgtb <- RBioplot::rightside_y(plt)
      } else { # no right side y-axis
        pltgtb <- plt
      }

      ## export the file and draw a preview
      ggsave(filename = paste(deparse(substitute(object)), ".", names(vip_list)[j], ".vip.pdf", sep = ""), plot = pltgtb,
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      cat("Done!\n") # final message
      grid.draw(pltgtb) # preview
    }
  }

  ## output
  # important features lsit
  ipf_list <- vector(mode = "list", length = length(vip_list))
  ipf_list[] <- foreach(m = 1:length(vip_list)) %do% {
    as.character(vip_list[[m]][which(vip_list[[m]]$VIP > vip.alpha), 1])
  }
  names(ipf_list) <- names(vip_list)

  # output
  out <- list(vip_summary = vip_list, features_above_alpha = ipf_list, vip_raw = vip_raw_list, ncomp = object$ncomp)
  assign(paste(deparse(substitute(object)), "_plsda_vip_summary_list", sep = ""), out, envir = .GlobalEnv)
}


#' @title rbioFS_plsda_roc_auc()
#'
#' @description ROC-AUC analysis and ploting for plsda model
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.smooth If to smooth the curves. Uses binormal method to smooth the curves. Default is \code{FALSE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is \code{length(plot.comps)}.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}. Note: doesn't seem to be necessasry as PLS-DA always has at least two y classes.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"1 - specificity"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"sensitivity"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @return Prints AUC values in the console. And a pdf file for ROC plot
#' @details Uses pROC module to calculate ROC.
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioFS_plsda_roc_auc(object = model_binary, rocplot = TRUE, plot.comps = 1:2)
#' }
#' @export
rbioFS_plsda_roc_auc <- function(object, rocplot = TRUE,
                                 plot.comps = 1:object$ncomp,
                                 plot.smooth = FALSE,
                                 multi_plot.ncol = length(plot.comps), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                                 plot.rightsideY = TRUE,
                                 plot.SymbolSize = 2, plot.display.Title = TRUE, plot.titleSize = 10,
                                 plot.fontType = "sans",
                                 plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                 plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                 plot.legendSize = 9,
                                 plot.Width = 170, plot.Height = 150){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.\n")
  if(plot.smooth) cat("ROC smooth: ON.\n") else cat("ROC smooth: OFF.\n")

  ## calcuate ROC-AUC
  pred_raw <- predict(object, newdata = object$inputX)
  outcome <- object$inputY

  roc_dfm_list <- vector(mode = "list", length = length(plot.comps))
  roc_dfm_list[] <- foreach(i = 1:length(plot.comps)) %do% {
    tmp_pred_raw <- pred_raw[, , i]
    out <- foreach(j = 1:length(levels(outcome)), .combine = "rbind") %do% {
      response <- outcome
      levels(response)[-j] <- "others"
      predictor <- as.matrix(tmp_pred_raw[, j], ncol = 1)
#      pred <- ROCR::prediction(predictions = predictor, labels = response)
#      perf <- ROCR::performance(prediction.obj = pred, "tpr", "fpr")
      splt <- split(predictor, response)  # split function splist array according to a factor
      controls <- splt$others
      cases <- splt[[levels(outcome)[j]]]
      perf <- tryCatch(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth),
                       error = function(err){
                         cat("Curve not smoothable. Proceed without smooth.\n")
                         pROC::roc(controls = controls, cases = cases, smooth = FALSE)
                       })
      if (length(levels(outcome)) == 2){
        cat(paste0("comp ", i, " AUC - ", levels(outcome)[j], ": ", perf$auc, "\n"))
      } else {
        cat(paste0("comp ", i, " AUC - ", levels(outcome)[j], " (vs Others): ", perf$auc, "\n"))
      }

#      fpr <- as.numeric(unlist(perf@x.values))
#      tpr <- as.numeric(unlist(perf@y.values))
      fpr <- 1 - perf$specificities
      tpr <- perf$sensitivities
      mtx <- cbind(fpr, tpr)
      if (length(levels(outcome)) == 2){
        df <- data.frame(mtx, group = rep(levels(outcome)[j], times = nrow(mtx)), row.names = NULL)
      } else {
        df <- data.frame(mtx, group = rep(paste0(levels(outcome)[j], " (vs Others)"), times = nrow(mtx)), row.names = NULL)
      }
      return(df)
    }
  }
  names(roc_dfm_list) <- paste0("comp ", 1:length(plot.comps))

  ## plot
  if (rocplot){
    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.roc.pdf...", sep = ""))  # initial message

    # extract comp variance distribution
    varpp_x <- 100 * object$Xvar / object$Xtotvar
    boxdfm_x <- data.frame(comp_x = as.numeric(gsub("Comp", "", names(varpp_x))), varpp_x = varpp_x)
    var_percentage_x <- varpp_x[paste0("Comp ", plot.comps)] # extract the proportion of variance for the selected PCs
    comp_axis_lbl <- paste("comp ", plot.comps, " (", round(var_percentage_x, digits = 2), "%)", sep = "")

    plt_list <- vector(mode = "list", length = length(plot.comps))
    plt_list[] <- foreach(k = 1:length(plot.comps)) %do% {
      plt <- ggplot(data = roc_dfm_list[[k]], aes(x = fpr, y = tpr, group = group, colour = group)) +
        geom_line(aes(linetype = group)) +
        geom_point(aes(shape = group), size = plot.SymbolSize) +
        geom_abline(intercept = 0) +
        ggtitle(ifelse(plot.display.Title, comp_axis_lbl[k], NULL)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
              legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

      if (plot.rightsideY & length(plot.comps) == 1){
        plt <- RBioplot::rightside_y(plt)
      }
      plt
    }
    names(plt_list) <- names(roc_dfm_list)[plot.comps]
    if (length(plot.comps) > 1) {
      if (multi_plot.ncol * multi_plot.nrow < length(plot.comps)){
        stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. Make sure they match the length of plot.comps.")
      } else {
        plt <- RBioplot::multi_plot_shared_legend(plt_list, ncol = multi_plot.ncol, nrow = multi_plot.nrow, position = multi_plot.legend.pos)
      }
    }

    # save
    grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.roc.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    cat("Done!\n")
  }
  assign(paste(deparse(substitute(object)), "_plsda_roc_list", sep = ""), roc_dfm_list, envir = .GlobalEnv)
}


#' @title rbioFS_plsda_predict
#'
#' @description Prediction function for PLS-DA analysis. The function calculates the predicted value for unknown sample data using the input PLS-DA model.
#' @param object A \code{rbiomvr} or \code{mvr} object.
#' @param comps  Number of PLS-DA components used in the model. Default is \code{object$ncomp}.
#' @param newdata Input data to be classified. Make sure it is a \code{matrix} class and has the same variables as the model, i.e. same number of columns as the training data.
#' @param prob.method Method to calculate classification probability. Options are \code{"softmax"} and \code{"Bayes"}. See details for more information. Default is \code{"Bayes"}.
#' @param threshold  Classification threshold. Should be a number between \code{0} and \code{1}. Default is \code{0.2}.
#' @param predplot If to generate a prediction value plot. Default is \code{TRUE}.
#' @param plot.sampleLabel.type If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param plot.sampleLabel.vector Set only when \code{plot.sampleLabel.type} is not set to \code{"none"}, a character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param plot.sampleLabelSize Only set when \code{plot.sampleLabel.type} is not \code{"none"}. The size of the sample label. Default is \code{2}.
#' @param plot.sampleLabel.padding Set only when \code{plot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is \code{length(plot.comps)}.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}. Note: doesn't seem to be necessasry as PLS-DA always has at least two y classes.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.unclassifiedColour Colour for the unclassified samples. Default is \code{"gray"}.
#' @param plot.classifiedColour Colour for the classified samples. Default is \code{"red"}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Samples"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.xAngle The rotation angle (degrees) of the x axis marks. Default is \code{0} - horizontal.
#' @param plot.xhAlign The horizontal alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.xvAlign The vertical alignment type of the x axis marks. Options are \code{0}, \code{0.5} and \code{1}. The default value at \code{0} is especially useful when \code{xAngle = 90}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"Predicted values"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @return  A \code{prediction} obejct, as well as pdf figure file for predicted values if \code{predplot = TRUE}.
#' @details Regarding \code{threshold}, the value should between \code{0} and \code{1}. It's the flank region around the dummified classification values \code{0} (i.e. "control") and \code{1} (i.e. "case").
#' This is only for information sake.
#'
#' The "Bayes" method uses the klaR package implementation of naive Bayes algorithm, based on the Bayes theorem: \code{P(A|B) = P(B|A)P(A)/P(B)}.
#' Specifically, the equation models classification against predicted values from pls-da model and training data: \code{P(class|predicted values) = P(predicted values|)P(class)/P(predicted values)}.
#' The resulted Bayesian classification model is then applied to newdata, thereby prosterior probabilites are calcuated.
#'
#' The "softmax" method uses the equation: \code{probability = exp(predicted.value) / sum(exp(predicted.values))}. This method doesn't relay on training data.
#'
#' The conventional wisdom is to use "Bayes" method for unbalanced classification, and {"softmax"} for balanced situation.
#'
#' For classification plot, the output \code{prediction} object should be used with function \code{\link{rbioFS_plsda_classification()}}.
#' @import ggplot2
#' @import pls
#' @import ggrepel
#' @importFrom klaR NaiveBayes
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioFS_plsda_predict(object = new_model_optm, newdata = newdata, prob.method = "Bayes",
#'                      plot.sampleLabel.type = "none", plot.sampleLabel.vector = NULL, plot.sampleLabel.padding = 0.5,
#'                      multi_plot.ncol = length(levels(object$inputY)), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
#'                      plot.SymbolSize = 2, plot.display.Title = TRUE, plot.titleSize = 10,
#'                      plot.fontType = "sans", plot.unclassifiedColour = "gray", plot.classifiedColour = "red",
#'                      plot.xLabel = "Samples", plot.xLabelSize = 10, plot.xTickLblSize = 6, plot.xTickItalic = FALSE,
#'                      plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
#'                      plot.yLabel = "Predicted values", plot.yLabelSize = 10, plot.yTickLblSize = 10,
#'                      plot.legendSize = 9,
#'                      plot.Width = 170, plot.Height = 150)
#' }
#' @export
rbioFS_plsda_predict <- function(object, comps = object$ncomp, newdata, prob.method = "Bayes",
                                 threshold = 0.2,
                                 predplot = TRUE,
                                 plot.sampleLabel.type = "none", plot.sampleLabel.vector = NULL,
                                 plot.sampleLabelSize = 2, plot.sampleLabel.padding = 0.5,
                                 multi_plot.ncol = length(levels(object$inputY)), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                                 plot.rightsideY = TRUE,
                                 plot.SymbolSize = 2, plot.display.Title = TRUE, plot.titleSize = 10,
                                 plot.fontType = "sans", plot.unclassifiedColour = "gray", plot.classifiedColour = "red",
                                 plot.xLabel = "Samples", plot.xLabelSize = 10, plot.xTickLblSize = 6, plot.xTickItalic = FALSE,
                                 plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                                 plot.yLabel = "Predicted values", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                 plot.legendSize = 9,
                                 plot.Width = 170, plot.Height = 150){
  ## argument check
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.\n")
  if (!class(newdata) %in% "matrix") stop("newdata has to be a matrix object. \n")
  if (ncol(newdata) != ncol(object$inputX)) stop("newdata needs to have the same number of variables, i.e. columns, as the object. \n")
  if (!prob.method %in% c("softmax", "Bayes")) stop("Probability method should be either \"softmax\" or \"Bayes\".\n")
  if (!tolower(plot.sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\". \n")
  if (tolower(plot.sampleLabel.type != "none")){
    if (!is.null(plot.sampleLabel.vector) & length(plot.sampleLabel.vector) != nrow(newdata)) {
      stop("plot.sampleLabel.vector has to be the same length as nrow(newdata). \n")}
  }

  ## prediction and class assign
  rownames(newdata) <- 1:nrow(newdata)
  pred <- predict(object = object, ncomp = comps, newdata = newdata, type = "response")
  if (is.null(plot.sampleLabel.vector)){
    sample.label <- seq(nrow(newdata))
  } else {
    sample.label <- plot.sampleLabel.vector
  }

  ## classification and probability calculation
  pred_mtx <- pred[, ,1]
  if (is.null(dim(pred_mtx))) {  # if only one sample
    pred_mtx <- t(as.matrix(pred_mtx))
  }
  rownames(pred_mtx) <- plot.sampleLabel.vector
  if (prob.method == "Bayes"){ # Naive Bayes probability calculation and prediction for sample data
    training_mtx <- as.matrix(object$centerX$centerX)
    trainingpred <- predict(object = object, ncomp = comps, newdata = training_mtx, type = "response")
    bayes.prob <- klaR::NaiveBayes(x = trainingpred,
                                   grouping = object$inputY, usekernel = TRUE)
    bayes.prob$train.posterior <- predict(bayes.prob)$posterior  # calcuate posterior probability for
    bayes.prob$x <- NULL

    bayespred <- predict(object = bayes.prob, newdata = pred_mtx)
    prob <- bayespred$posterior

    prob <- foreach(i = 1:nrow(pred_mtx), .combine = "rbind") %do% {
      prob_dfm <- data.frame(Sample = rep(rownames(prob)[i], times = ncol(prob)),
                             Class = colnames(prob), Probability = prob[i, ], stringsAsFactors = FALSE)
      prob_dfm$Sample <- factor(prob_dfm$Sample, unique(prob_dfm$Sample))
      prob_dfm$repel.label.pos <- rev(cumsum(rev(prob_dfm$Probability)) - rev(prob_dfm$Probability) / 2)  # calculate the repel lable position, seemingly from bottom up
      prob_dfm$precent.label <- paste0(signif(prob_dfm$Probability, 4) * 100, "%")
      prob_dfm$Class <- factor(prob_dfm$Class, unique(prob_dfm$Class))
      return(prob_dfm)
    }
  } else {
    group <- colnames(pred_mtx)
    # calcuate probability
    prob_mtx <- apply(pred_mtx, 1, FUN = function(x) exp(x) / sum(exp(x)))
    rownames(prob_mtx) <- group

    # constuct plot dataframe list
    prob <- foreach(i = 1:ncol(prob_mtx), .combine = "rbind") %do% {
      prob_dfm <- data.frame(Sample = rep(rownames(pred_mtx)[i], times = length(group)),
                             Class = group, Probability = prob_mtx[, i], stringsAsFactors = FALSE)
      prob_dfm$Sample <- factor(prob_dfm$Sample, unique(prob_dfm$Sample))
      prob_dfm$repel.label.pos <- rev(cumsum(rev(prob_dfm$Probability)) - rev(prob_dfm$Probability) / 2)  # calculate the repel lable position, seemingly from bottom up
      prob_dfm$precent.label <- paste0(signif(prob_dfm$Probability, 4) * 100, "%")
      prob_dfm$Class <- factor(prob_dfm$Class, unique(prob_dfm$Class))
      return(prob_dfm)
    }
  }

  ## prediction plot
  if (predplot){
    # prepare prediction value list
    predlist <- vector(mode = "list", length = length(levels(object$inputY)))
    predlist[] <- foreach(i = 1:length(levels(object$inputY))) %do% {
      preddfm <- data.frame(sample = as.integer(rownames(newdata)), sample.label = sample.label, predicted.value = pred[, i,])
      preddfm$classification <- sapply(preddfm$predicted.value, FUN = function(x)ifelse(x > 1 - threshold & x < (1 + threshold), levels(object$inputY)[i], ifelse(x > - threshold & x < threshold, "rest", "undermined")))
      preddfm$`Within threshold` <- ifelse(preddfm$classification == "undermined", "N", "Y")
      preddfm$`Within threshold` <- factor(preddfm$`Within threshold`, levels = c("Y", "N"))
      return(preddfm)
    }
    names(predlist) <- levels(object$inputY)

    # plot
    if (plot.sampleLabel.type != "none" & is.null(plot.sampleLabel.vector)){  # message when no label vector is provided.
      cat("plot.sampleLabel.vector not provided. Proceed with row numbers as sampole labels.\n")
    }
    cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.predict.pdf...", sep = ""))  # initial message
    plt_list <- vector(mode = "list", length = length(levels(object$inputY)))
    plt_list[] <- foreach(j = 1:length(levels(object$inputY))) %do% {
      pltdfm <- predlist[[j]]
      plt <- ggplot(data = pltdfm)  # base plot

      if (plot.sampleLabel.type != "none"){  # labels
        if (tolower(plot.sampleLabel.type) == "direct"){
          plt <- plt + geom_text(aes(x = sample, y = predicted.value, label = sample.label, colour = `Within threshold`), size = plot.SymbolSize)
        } else if (tolower(plot.sampleLabel.type) == "indirect") {
          plt <- plt + geom_point(alpha = 0.6, size = plot.SymbolSize, aes(x = sample, y = predicted.value, colour = `Within threshold`)) +
            geom_text_repel(data = pltdfm[!pltdfm$classification %in% c("rest", "undetermined"), ], aes(x = sample, y = predicted.value, label = sample.label),
                            point.padding = unit(plot.sampleLabel.padding, "lines"), size = plot.sampleLabelSize, show.legend = FALSE)
        }
      } else {
        plt <- plt + geom_point(alpha = 0.4, size = plot.SymbolSize, aes(x = sample, y = predicted.value, colour = `Within threshold`))
      }

      if (length(unique(pltdfm$`Within threshold`)) == 1){
        if (unique(pltdfm$`Within threshold`) == "N"){
          plt <- plt +
            scale_color_manual(values = plot.unclassifiedColour)
        } else if (unique(pltdfm$`Within threshold`) == "Y"){
          plt <- plt +
            scale_color_manual(values = plot.classifiedColour)
        }
      } else {
        plt <- plt +
          scale_color_manual(values = c(plot.classifiedColour, plot.unclassifiedColour))
      }

      plt <- plt +
        ggtitle(ifelse(plot.display.Title, levels(object$inputY)[j], NULL)) +
        #       scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(breaks = c(-threshold, 0, threshold, 1 - threshold, 1, 1 + threshold)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        geom_hline(yintercept = c(0, 1)) +
        geom_hline(yintercept = c(-threshold, threshold, 1 - threshold, 1 + threshold), linetype = "dashed") +
        #       geom_ribbon(aes(x = sample, ymax = threshold, ymin = -threshold), fill = "pink", alpha = 0.4) +  # colour between lines: lower
        #       geom_ribbon(aes(x = sample, ymax = 1 + threshold, ymin = 1 - threshold), fill = "pink", alpha = 0.4) +  # colour between lines: higher
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              legend.position = "bottom",
              #              legend.title = element_blank(),
              legend.text = element_text(size = plot.legendSize),
              legend.key = element_blank(),
              axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle,
                                         hjust = plot.xhAlign, vjust = plot.xvAlign),
              axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))

      if (plot.rightsideY & length(plt_list) == 1){
        plt <- RBioplot::rightside_y(plt)
      }

      return(plt)
    }
    names(plt_list) <- levels(object$inputY)

    if (length(levels(object$inputY)) > 1) {
      if (multi_plot.ncol * multi_plot.nrow < length(levels(object$inputY))){
        stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. Make sure they match the number of y groups.")
      } else {
        plt <- RBioplot::multi_plot_shared_legend(plt_list, ncol = multi_plot.ncol, nrow = multi_plot.nrow, position = multi_plot.legend.pos)
      }
    }

    # save
    grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.predict.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    cat("Done!\n")
  }

  ## export
  out <- list(predicted.value = pred_mtx, probability.summary = prob, probability.method = prob.method)
  class(out) <- "prediction"
  assign(paste(deparse(substitute(object)), "_plsda_predict", sep = ""), out, envir = .GlobalEnv)
}


#' @title rbioFS_plsda_classplot
#'
#' @description Classification plot function for PLS-DA analysis. The function uses \code{prediction} object generated from \code{\link{rbioFS_plsda_predict}} to generate classification probablity pie charts.
#' @param pred.object A \code{prediction} object, which can be obtained from funciton \code{\link{rbioFS_plsda_predict}}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is \code{nrow(pred.obj)}.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.Title Plot title. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.probLabelSize The size of the sample label. Default is \code{2}.
#' @param plot.probLabel.padding Set only when \code{plot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.probLabel.outside If to put the probability label outside of the pies. Default is \code{"TRUE"}.
#' @param plot.probLabel.outside.nudge Set only when \code{plot.probLabel.outside = TRUE}, adjustment to nudge the starting position of each label. Default is \code{"1.5"}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @return  A \code{classification} obejct with classification probability summary for each sample, as well as pdf figure file fif \code{classplot = TRUE}.
#' @details The function operates in conjunction with the prediction function \code{\link{rbioFS_plsda_predict}}, to which the sample(s) of intestested is provided.
#' @import ggplot2
#' @import ggrepel
#' @importFrom grid grid.newpage grid.draw
#' @examples
#' \dontrun{
#' rbioFS_plsda_classplot(pred.obj = new_model_optm_plsda_predict, multi_plot.ncol = 4, multi_plot.nrow = 4, plot.probLabelSize = 2)
#' }
#' @export
rbioFS_plsda_classplot <- function(pred.obj,
                                   multi_plot.ncol = nrow(pred.obj), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                                   multi_plot.stripLblSize = 10,
                                   plot.Title = NULL, plot.titleSize = 10,
                                   plot.probLabelSize = 5, plot.probLabel.padding = 0,
                                   plot.probLabel.outside = FALSE, plot.probLabel.outside.nudge = 1.5,
                                   plot.fontType = "sans",
                                   plot.legendSize = 9,
                                   plot.Width = 170, plot.Height = 150){
  ## check arguments
  if (!any(class(pred.obj) %in% "prediction")) stop("pred.obj needs to be a  \"prediction\" class. Use functions like rbioFS_plsda_predict() to generate one.\n")

  ## plot
  if (multi_plot.ncol * multi_plot.nrow < nrow(pred.obj$predicted.value)){
    stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. Make sure they match the number of y groups.\n")
  }
  cat(paste("Plot being saved to file: ", deparse(substitute(pred.obj)),".plsda.classification.pdf...", sep = ""))  # initial message
  plt <- ggplot(pred.obj$probability.summary, aes(x = "", y = Probability, fill = Class)) +
    geom_col(width = 1, colour = "black", alpha = 0.8) +
    ggtitle(plot.Title) +
    xlab(NULL) +
    ylab(NULL) +
    facet_wrap(~Sample, nrow = multi_plot.nrow, ncol = multi_plot.ncol) +
    theme(strip.background = element_rect(fill = NA, colour = "black"),  # no strip background colour
          strip.text = element_text(face = "bold", size = multi_plot.stripLblSize),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
          axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
          axis.text.x = element_blank(),
          legend.position = multi_plot.legend.pos, legend.title = element_blank(),
          legend.text = element_text(size = plot.legendSize),
          legend.key = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())

  if (plot.probLabel.outside){
    plt <- plt +
      geom_label_repel(aes(label = precent.label, y = repel.label.pos), point.padding = unit(plot.probLabel.padding, "lines"),
                       show.legend = FALSE, nudge_x = plot.probLabel.outside.nudge, size = plot.probLabelSize)
  } else {
    plt <- plt +
      geom_label_repel(aes(label = precent.label, y = repel.label.pos), point.padding = unit(plot.probLabel.padding, "lines"),
                       show.legend = FALSE, size = plot.probLabelSize)
  }

  plt <- plt + coord_polar("y")

  # save
  grid.newpage()
  ggsave(filename = paste(deparse(substitute(pred.obj)),".plsda.classification.pdf", sep = ""), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  cat("Done!\n")

  invisible(plt)
}

