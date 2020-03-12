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


#' @title rbioClass_plsda
#'
#' @description PLS-DA modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class or a character vector.
#' @param ncomp Number of components to be used for modelling. Default is \code{length(unique(y)) - 1}.
#' @param method PLS-DA modelling method. Four PLSR algorithms are available: the kernel algorithm ("kernelpls"), the wide kernel algorithm ("widekernelpls"),
#'               SIMPLS ("simpls") and the classical orthogonal scores algorithm (also known as NIPALS) ("oscorespls"). Default is the popular \code{"simpls"}.
#' @param scale Logical, whether to scale the data or not. Default is \code{TRUE}.
#' @param validation Cross validation methods. Options are "none", "CV" (fold), "LOO" (leave-one-out). Default is \code{"CV"}.
#' @param segments Set only when \code{validation = "CV"}, the number of segement to be set. Default is \code{10}.
#' @param segments.type Method to set up the segments. Options are \code{"random", "consecutive", "interleaved"}. Default is \code{"random"}.
#' @param jackknife If to use jack-knife procedure. Default is \code{TRUE}.
#' @param ... Additional arguments for \code{mvr} function from \code{pls} pacakge.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details The center and scale processes are handled by \code{\link{center_scale}} function.
#'          Therefore, the \code{center} and \code{scale} arguments \code{pls::mvr} are set to \code{FALSE}.
#'
#'          The function always centers the data.
#'
#' @return Returns a PLS-DA model object, with classes "mvr" and "rbiomvr".
#'
#'        Additional items for \code{rbiomvr} object to \code{mvr} object from pls package:
#'
#'        \code{centerX}: centered X data if applicable, with scale option.
#'
#'        \code{dummy_y}: dummified y
#'
#'        \code{inputX}: raw input predictor data.
#'
#'        \code{inputY}: input group labels.
#'
#'        \code{validation_method}
#'
#'        \code{validation_segments_type}
#'
#'        \code{model.type}
#'
#' @details The data is always centered prior to modelling, i.e. \code{x - col.mean}. Thus no "center = TRUE/FALSE" option is provided.
#'
#'          For \code{ncomp} value, the default (full model) compares feature number with \code{class number - 1}, instead of observation number.
#'
#'          For sequencing data, the input x needs to be either tranformed with the function like \code{clr_ilr_transfo()} from \code{RBioArray} package,
#'          or normalized using methods like "TMM" or "RLE" implemented in \code{edgeR} pacakge.
#'
#' @importFrom pls plsr
#' @examples
#' \dontrun{
#' rbioClass_plsda(x, y, ncomp = 20)
#' }
#' @export
rbioClass_plsda <- function(x, y, ncomp = length(unique(y)) - 1, method = "simpls",
                            scale = TRUE,
                            validation = c("none", "CV", "LOO"),
                            segments = 10, segments.type = "random",
                            jackknife = TRUE, ...,
                            verbose = TRUE){
  ## check arguments
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x))) stop("x needs to be a matrix, data.frame or vector.")
  if (class(x) == "data.frame" | is.vector(x)){
    if (verbose) cat("x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }
  if (class(y) != "factor"){
    if (verbose) cat("y is converted to factor. \n")
    y <- factor(y, levels = unique(y))
  }
  if (is.null(ncomp))stop("please set the ncomp number.")

  ## data processing
  # X
  if (verbose) cat(paste0("data centered with the scale option ", ifelse(scale, "\"ON\" ", "\"OFF\" "), "prior to modelling..."))
  centered_X <- center_scale(x, scale = scale)  # center data with the option of scaling
  if (verbose) cat("DONE!\n")
  X <- centered_X$centerX
  # Y
  Y <- dummy(y)
  # constract dataframe for modelling
  model_dfm <- data.frame(n = paste("row", 1:nrow(Y), sep = ""))
  model_dfm$Y <- Y  # y is the y matrix as a whole, not the content of y
  model_dfm$X <- X  # x is the x matrix as a whole, not the content of x

  ## modelling
  if (validation == "LOO") {
    segments.type <- NA
    out_model <- plsr(Y ~ X, data = model_dfm, ncomp = ncomp,
                      method = method, scale = FALSE, center = FALSE,
                      validation = validation, jackknife = TRUE, ...)
  } else {
    out_model <- plsr(Y ~ X, data = model_dfm, ncomp = ncomp,
                      method = method, scale = FALSE, center = FALSE,
                      validation = validation, segments = segments, segments.type = segments.type, jackknife = TRUE, ...)
  }
  out_model$centerX <- centered_X
  out_model$dummy_y <- Y
  out_model$inputX <- x
  out_model$inputY <- y
  out_model$validation_method <- validation
  out_model$validation_segments_type <- segments.type
  out_model$model.type <- "classification"

  class(out_model) <- c("rbiomvr", "mvr")
  return(out_model)
}


#' @title rbioClass_plsda_tuplot
#'
#' @description T-U plot function for PLS-DA models.
#' @param object A \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
#' @param comps Integer vector. Components to plot. The index of the components are intergers. The vector length should be between 1 and the total number of components, inclusive. Can be Default is \code{c(1, 2)}.
#' @param multi_plot.ncol Set only when \code{length(comps) > 1}, number of columns on one figure page. Default is \code{length(comps)}.
#' @param multi_plot.nrow Set only when \code{length(comps) > 1}, number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.rightsideY If to show the right side y-axis. Only applicble when the length of \code{comps} is less than 2, inclusive. Default is \code{FALSE}. Note: the right side Y is ignored when \code{length(comps) > 1}
#' @param plot.sampleLabel.type If to show the sample labels on the graph. Options are \code{"none"}, \code{"direct"} and \code{"indirect"}. Default is \code{"none"}.
#' @param plot.sampleLabel.vector Set only when \code{plot.sampleLabel.type} is not set to \code{"none"}, a character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param plot.sampleLabelSize Only set when \code{plot.sampleLabel.type} is not \code{"none"}. The size of the sample label. Default is \code{2}.
#' @param plot.sampleLabel.padding Set only when \code{plot.sampleLabel.type = "indirect"}, the padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.Title tuplot title. Default is \code{NULL}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width tuplot width. Default is \code{170}.
#' @param plot.Height tuplot height. Default is \code{150}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#' rbioClass_plsda_tuplot(new_model, comps = c(1, 2, 3))
#' }
#' @export
rbioClass_plsda_tuplot <- function(object, comps = 1, multi_plot.ncol = length(comps), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                                plot.rightsideY = TRUE,
                                plot.sampleLabel.type = "none", plot.sampleLabel.vector = NULL,
                                plot.sampleLabelSize = 2, plot.sampleLabel.padding = 0.5,
                                plot.SymbolSize = 5, plot.Title = NULL,
                                plot.fontType = "sans",
                                plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                plot.legendSize = 9,
                                plot.Width = 170, plot.Height = 150,
                                verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr"))) stop("object needs to be either a \"rbiomvr\" class.")
  if (length(comps) > object$ncomp) stop("comps length exceeded the maximum comp length.")
  if (!all(comps %in% seq(object$ncomp))) stop("comps contain non-existant comp.")
  # if (!tolower(plot.sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\".")
  if (plot.rightsideY){
    if (length(comps) > 1){
      cat("right side y-axis ignored for multi-plot figure.\n")
    }
  }
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")
  plot.sampleLabel.type <- match.arg(tolower(plot.sampleLabel.type), c("none", "direct", "indirect"))


  ## extract and construt t-u plot dataframe
  tu_dfm <- data.frame(y = object$inputY)

  varpp_x <- 100 * object$Xvar / object$Xtotvar
  var_percentage_x <- varpp_x[paste0("Comp ", comps)] # extract the proportion of variance for the selected comps

  ## plot
  if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.tuplot.pdf...", sep = ""))  # initial message
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
  # grid.newpage()
  ggsave(filename = paste(deparse(substitute(object)),".plsda.tuplot.pdf", sep = ""), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  if (verbose) cat("Done!\n")

  ## return invinsible object
  invisible(plt)
}


#' @title rbioClass_plsda_q2r2()
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
#' @param plot.Width Plot width. Default is \code{170}.
#' @param plot.Height Plot height. Default is \code{150}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#' rbioClass_plsda_q2r2(object = new_model, multi_plot.ncol = 2, multi_plot.nrow = 2, intercept = TRUE)
#' }
#' @export
rbioClass_plsda_q2r2 <- function(object, intercept = TRUE, q2r2plot = TRUE,
                              plot.display.Title = TRUE,
                              multi_plot.ncol = length(dimnames(object$coefficients)[[2]]), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                              plot.rightsideY = TRUE, plot.fontType = "sans",
                              plot.SymbolSize = 2,
                              plot.xLabel = "Components", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                              plot.yLabel = "R2 & Q2", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                              plot.legendSize = 9,
                              plot.Width = 170, plot.Height = 150,
                              verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr", 'mvr'))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.")
  if (is.null(object$validation) || is.null(object$validation$coefficients)) stop("'object' was not fit with jackknifing enabled")  # from pls pacakge
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")
  multi_plot.legend.pos <- match.arg(tolower(multi_plot.legend.pos), c("bottom"))

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
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.q2r2plot.pdf...", sep = ""))  # initial message
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
    if (verbose) cat("Done!\n")
  }
}


#' @title randomiz.test
#'
#' @description Randomization function from \code{pls} pacakge. For \code{\link{rbioClass_plsda_ncomp_select}} function.
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


#' @title rbioClass_plsda_ncomp_select()
#'
#' @description Optimal number of components selection for PLS-DA model, with RMSEP plot funcitonality. Selection methods are modified based on \code{selectNcomp()} from \code{pls} pacakge.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ... Additional argument for \code{RMSEP} function from \code{pls} package.
#' @param ncomp.selection.method Optimal numbers of components selection method. Options are \code{"min"}, \code{"1err"}, and \code{"randomization"}. Default is \code{"1err"}.
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
#' @param plot.Width Plot width. Default is \code{170}.
#' @param plot.Height Plot height. Default is \code{150}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Prints the selected number of components for each y class. Returns RMSEP values for each y class to the environment as a \code{rbiomvr_ncomp_select} object,
#'         as well as a pdf file for the RMSEP plot if \code{rmsepplot = TRUE}.
#' @details The RMSEP figure shows both CV estimates and adjusted CV estimates, which is CV estimiates corrected for bias.
#'          Three methods are used for components number selection: \code{"min"} simply chooses the number of components to
#'          reach te minimum RMSEP; \code{"1err"} chooses the number of components when its RMSEP first reaches minimum as well as within one standard error;
#'          For "randomization", see the help file for \code{selectNcomp()} function from  \code{pls} pacakge.
#' @import ggplot2
#' @import foreach
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @importFrom pls RMSEP
#' @examples
#' \dontrun{
#' rbioClass_plsda_ncomp_select(new_model,  multi_plot.ncol = 2, multi_plot.nrow = 2, plot.optm.ncomp.line = T,
#'          ncomp.selection.method = "randomization", randomization.nperm = 999, randomization.alpha = 0.05)
#' }
#' @export
rbioClass_plsda_ncomp_select <- function(object, ...,
                                         ncomp.selection.method = c("1err", "min", "randomization"), randomization.nperm = 999, randomization.alpha = 0.05,
                                         rmsepplot = TRUE,
                                         plot.rightsideY = TRUE,
                                         plot.optm.ncomp.line = TRUE,
                                         multi_plot.ncol = length(dimnames(object$coefficients)[[2]]), multi_plot.nrow = 1,
                                         multi_plot.legend.pos = c("bottom", "top", "left", "right"),
                                         plot.display.Title = TRUE,
                                         plot.SymbolSize = 2,
                                         plot.fontType = "sans",
                                         plot.xLabel = "Components", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                         plot.yLabel = "RMSEP", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                         plot.legendSize = 9,
                                         plot.Width = 170, plot.Height = 150,
                                         verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("mvr", "rbiomvr"))) stop("object has to be a mvr or rbiomvr class.")
  if (!"validation" %in% names(object)) stop("PLS-DA model has to include Cross-Validation.")
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")
  # if (!tolower(ncomp.selection.method) %in% c("min", "1err", "randomization")) stop("ncomp.selection.method needs to be \"min\", \"1err\", or \"randomization\" exactly.")
  ncomp.selection.method <- match.arg(tolower(ncomp.selection.method), c("1err", "min", "randomization"))
  multi_plot.legend.pos <- match.arg(tolower(multi_plot.legend.pos), c("bottom", "top", "left", "right"))

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

      if (ncomp.selection.method == "1err"){
        dfm_cv$SD <- residsds
        min_rmsep_cv_idx <- which.min(dfm_cv$value)  # index for the minimum mean oob feature group
        sd_min_cv <- dfm_cv$SD[min_rmsep_cv_idx]  # oob SD for the feature group above
        cv_x <- min(with(dfm_cv, which(value <= (value[min_rmsep_cv_idx] + sd_min_cv))))  # 1err minimum selection

        dfm_adjcv$SD <- residsds
        min_rmsep_adjcv_idx <- which.min(dfm_adjcv$value)  # index for the minimum mean oob feature group
        sd_min_adjcv <- dfm_adjcv$SD[min_rmsep_adjcv_idx]  # oob SD for the feature group above
        adjcv_x <- min(with(dfm_adjcv, which(value <= (value[min_rmsep_adjcv_idx] + sd_min_adjcv))))  # 1erminimum selection
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
  if (verbose) cat(paste0("Optimial ncomp (selection method: ", ncomp.selection.method, ") : \n", sep = ""))
  if (verbose) print(ncompsel_mtx)

  # plotting
  if (rmsepplot){
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.rmsepplot.pdf...", sep = ""))  # initial message
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
    # grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.rmsepplot.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    if (verbose) cat("Done!\n")
  }

  ## return RMSEP values
  rmsep_dfm_list$ncomp_selected <- ncompsel_mtx
  class(rmsep_dfm_list) <- "rbiomvr_ncomp_select"
  assign(paste(deparse(substitute(object)), "_plsda_ncomp_select", sep = ""), rmsep_dfm_list, envir = .GlobalEnv)
}


#' @export
print.rbiomvr_ncomp_select <- function(x, ...){
  cat(paste0("PLS-DA ncomp selection results:\n"))
  cat("\n")
  print(x$ncomp_selected)
  cat("\n")
}


#' @title rbioClass_plsda_perm()
#'
#' @description Permutation test for PLS-DA models.
#' @param object A \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ncomp Model complexity, i.e. number of components to model. Default is \code{object$ncomp}, i.e. maximum complexity.
#' @param adjCV If to use adjusted CV, i.e. CV adjusted for unbalanced data. Default is \code{FALSE}.
#' @param perm.method Permutation method. Options are \code{"by_y"} and \code{"by_feature_per_y"}. Default is \code{"by_y"}. See details below.
#' @param nperm Number of permutations to run. Default is \code{999}.
#' @param perm.plot Wether to produce a plot or not. Default is \code{TRUE}.
#' @param ... Additional argument for \code{\link{rbioUtil_perm_plot}}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return The function returns \code{CSV} files for all intermediate permutation RMSEP values as well as the p-value resutls.
#'
#' The results are also exported to the environment as a \code{rbiomvr_perm} object with the following items:
#'
#' \code{perm.method} Permutation method.
#'
#' \code{nperm} The number of permutation runs.
#'
#' \code{perm.stats} The stats metric used for the permutation test.
#'
#' \code{adjCV} If the adjusted cross-validation stats is used.
#'
#' \code{p.value.summary} P values for the permutation test.
#'
#' \code{poerm.results} The intermediate permutation results, i.e. stats for each permutation test run in a data.frame. \code{nperm = 0} is the original stats.
#'
#' A scatter plot is also generaeted when \code{perm.plot = TRUE}.
#'
#' @details The function uses RMSEP as the stats for comparing original model with permutatsions.
#'
#'          Usually, we use the optimized PLS-DA model for \code{object}, which can be obtained from functions \code{\link{rbioClass_plsda}} and \code{\link{rbioClass_plsda_ncomp_select}}.
#'
#'          Data for permutation are object$centerX$centerX, meaning centered X are used if applicable.
#'
#'          Permutation methods are according to:
#'
#'             Ojala M, Garriga GC. 2010. Permutation test for studying classifier performance.
#'             J Mach Learn Res. 11: 1833 - 63.
#'
#'          For \code{perm.method = "by_y"}, labels (i.e. y) are permutated. A non-signifianct model (permutation p value > alpha, i.e. 0.05) in this case means the data is independent from the groups.
#'
#'          For \code{perm.method = "by_feature_per_by"}, X is first subset by label (i.e.y) before permutating data for each feature.
#'          Since the permutation is done for the features WITHIN the group,
#'          the test actually evaluates if the model will produce significantly different performannce from the permutation models with the original "betweeen-features" relation (if any) disturbed.
#'          Therefore, A non-significant result (permutation p value > alpha, i.e. 0.05) means either the features are independent, or the model doesn't consider correlation between the features.
#'
#' @import ggplot2
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pls RMSEP
#' @examples
#' \dontrun{
#' rbioClass_plsda_perm(object = new_model_optm, perm.method = "by_y", nperm = 999, adjCV = TRUE, parallelComputing = TRUE)
#' }
#' @export
rbioClass_plsda_perm <- function(object, ncomp = object$ncomp, adjCV = FALSE,
                                 perm.plot = TRUE, ...,
                                 perm.method = c("by_y", "by_feature_per_y"), nperm = 999,
                                 parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                                 verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("rbiomvr"))) stop("object has to be a \"rbiomvr\" class.")
  if (!"validation" %in% names(object) || is.null(object$validation)) stop("PLS-DA model has to include Cross-Validation.")
  if (!all(ncomp %in% seq(object$ncomp))) stop("ncomp contain non-existant comp.")
  # if (!perm.method %in% c("by_y", "by_feature_per_y")) stop("perm.method needs to be either \"by_y\" or \"by_feature_per_y\".")
  if (length(nperm) != 1) stop("nperm can only contain one integer.")
  if (nperm %% 1 != 0) stop("nperm can only be integer. \n")
  if (nperm < 1) stop("nperm can only take interger equal to or greater than 1.")
  perm.method <- match.arg(tolower(perm.method), c("by_y", "by_feature_per_y"))
  if (parallelComputing){
    clusterType <- match.arg(clusterType, c("PSOCK", "FORK"))
  }
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")

  ## calcuate RMSEP and construct original RMSEP data frame
  rmsep <- pls::RMSEP(object)
  rmsep_dfm <- foreach(i = 1:dim(rmsep$val)[2], .combine = "rbind") %do% {
    dfm <- as.data.frame(t(rmsep$val[, i, ncomp + 1]))  # +1 becuase intercept
    if (adjCV){
      stats <- dfm[, "adjCV"]
    } else {
      stats <- dfm[, "CV"]
    }
    dfm <- data.frame(comparison = dimnames(rmsep$val)[[2]][i], comps = object$ncomp, RMSEP = stats,
                      adjCV = adjCV, row.names = NULL)
  }

  ## permutation test
  # permutation test functions
  by_y_func <- function(i){
    set.seed(i)
    perm_y <- object$inputY[sample(1:length(object$inputY))]  # sample label permutation
    perm_model <- rbioClass_plsda(x = object$centerX$centerX, y = factor(perm_y, levels = unique(perm_y)),
                                  ncomp = ncomp, scale = FALSE, validation = object$validation_method,
                                  segments = length(object$validation$segments),
                                  segments.type = object$out_model$validation_segments_type,
                                  verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
    perm_rmsep <- pls::RMSEP(perm_model)  # permutation model RMSEP

    perm_rmsep_dfm <- foreach(j = 1:dim(perm_rmsep$val)[2], .combine = "rbind") %do% {
      dfm <- as.data.frame(t(perm_rmsep$val[, j, ncomp + 1]))
      if (adjCV){
        stats <- dfm[, "adjCV"]
      } else {
        stats <- dfm[, "CV"]
      }
      dfm <- data.frame(nperm = i, comparison = dimnames(perm_rmsep$val)[[2]][j],
                        comps = ncomp, RMSEP = stats, adjCV = adjCV, row.names = NULL)
    }
    return(perm_rmsep_dfm)
  }
  by_feature_per_y_func <- function(i){
    set.seed(i)
    perm_x <- foreach(m = dimnames(rmsep$val)[[2]], .combine = "rbind") %do% {
      sub_dat <- object$centerX$centerX[which(object$inputY == m), ]  # subsetting the centre-scaled X by label (Y)
      sub_dat_colnames <- colnames(sub_dat)
      perm_dat <- sub_dat[, sample(1:ncol(sub_dat))]
      # perm_dat <- foreach(n = 1:ncol(sub_dat), .combine = "cbind") %do% {
      #   col <- sub_dat[sample(1:nrow(sub_dat)), n, drop = FALSE]  # permutation for each feature
      #   col
      # }
      colnames(perm_dat) <- sub_dat_colnames
      perm_dat
    }

    perm_model <- rbioClass_plsda(x = perm_x, y = object$inputY,
                                  ncomp = ncomp, scale = FALSE, validation = object$validation_method,
                                  segments = length(object$validation$segments),
                                  segments.type = object$out_model$validation_segments_type,
                                  verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
    perm_rmsep <- pls::RMSEP(perm_model)  # permutation model RMSEP

    perm_rmsep_dfm <- foreach(j = 1:dim(perm_rmsep$val)[2], .combine = "rbind") %do% {
      dfm <- as.data.frame(t(perm_rmsep$val[, j, ncomp + 1]))
      if (adjCV){
        stats <- dfm[, "adjCV"]
      } else {
        stats <- dfm[, "CV"]
      }
      dfm <- data.frame(nperm = i, comparison = dimnames(perm_rmsep$val)[[2]][j],
                        comps = ncomp, RMSEP = stats, adjCV = adjCV, row.names = NULL)
    }
    perm_rmsep_dfm
  }

  if (verbose) cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
  if (verbose) cat(paste0("Running permutation test using ", perm.method, " method ","with ", nperm, " permutations (speed depending on hardware configurations)..."))
  # consolidated permutation model RMSEP data frame
  if (!parallelComputing){
    if (perm.method == "by_y"){  # permutate label
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %do% by_y_func(i)
    } else {  # subset by label first, then permutate data per feature
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %do% by_feature_per_y_func(i)
    }
  } else {  # parallel computing
    # set up cpu cluster
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # permutation test
    if (perm.method == "by_y"){  # permutate label
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind", .packages = c("foreach", "pls", "RBioFS")) %dopar% by_y_func(i)
    } else {  # subset by label first, then permutate data per feature
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind", .packages = c("foreach", "pls", "RBioFS")) %dopar% by_feature_per_y_func(i)
    }
  }

  # permutation model data frame
  perm_results_p_val <- foreach(i = dimnames(rmsep$val)[[2]], .combine = "rbind") %do% {
    tmp <- perm_dfm[perm_dfm$comparison == i, ]
    orig_rmsep <- data.frame(nperm = 0, rmsep_dfm[rmsep_dfm$comparison == i, ], row.names = NULL)
    tmp <- rbind(orig_rmsep, tmp)
    rownames(tmp) <- NULL
    write.csv(file = paste0(deparse(substitute(object)),"_", i ,"_vs_rest" , "_perm_intermediate.csv"), tmp, row.names = FALSE)  # export full permutation results to working directory
    orig.val <- tmp$RMSEP[which(tmp$nperm == 0)]
    perm.val <- tmp$RMSEP[which(tmp$nperm != 0)]
    p.val <- (length(which(perm.val <= orig.val)) + 1) / (nperm + 1)
    out <- data.frame(comparison = i, original.RMSEP = orig.val, p.value = p.val,
                      row.names = NULL, stringsAsFactors = FALSE)
  }
  if (verbose) cat("Done!\n")

  ## output
  if (verbose) cat("\n")
  if (verbose) cat("Permutation test restuls: \n")
  if (verbose) print(perm_results_p_val)
  if (verbose) cat("\n")

  # env
  perm_stats <- rbind(data.frame(nperm = rep(0, times = nrow(rmsep_dfm)), comparison = rmsep_dfm$comparison, RMSEP = rmsep_dfm$RMSEP, row.names = NULL),
                      data.frame(nperm = perm_dfm$nperm, comparison = perm_dfm$comparison, RMSEP = perm_dfm$RMSEP, row.names = NULL))

  out <- list(perm.method = perm.method, nperm = nperm, perm.stats = "RMSEP", adjCV = adjCV,
              p.value.summary = perm_results_p_val, perm.results =  perm_stats, model.type = object$model.type)

  class(out) <- "rbiomvr_perm"
  assign(paste(deparse(substitute(object)), "_perm", sep = ""), out, envir = .GlobalEnv)

  # directory
  if (verbose) cat(paste0("Permutation RMSEP values stored in csv files with suffix: _perm_intermediate.csv. \n"))
  if (verbose) cat(paste0("Permutation results stored in csv file: ", paste(deparse(substitute(object)), "_perm.csv", sep = ""), ". \n"))
  write.csv(file = paste(deparse(substitute(object)), "_perm.csv", sep = ""), perm_results_p_val, row.names = FALSE)

  ## plot
  if (perm.plot){
    plsda_permutation_test <- out
    rbioUtil_perm_plot(perm_res = plsda_permutation_test, ..., verbose = verbose)
  }
}


#' @export
print.rbiomvr_perm <- function(x, ...){
  cat(paste0("PLS permutation results with ", x$nperm, " permutations:\n"))
  cat("\n")
  print(x$p.value.summary)
  cat("\n")
}


#' @title rbioClass_plsda_scoreplot
#'
#' @description scoreplot function for PLS-DA models.
#' @param object A \code{rbiomvr} or \code{mvr} object. Make sure the object is generated with a \code{validation} section.
#' @param y Set when object class is \code{mvr}. Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class or a character vector.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Returns a pdf file for scoreplot.
#' @details When \code{length(comps) == 1}, the function generates a scatter plot plotting sample vs score for the comp of interest. When \code{length(comps) == 2}, the function generates a scatter plot plotting the two comps of interest against each other. When \code{length(comps) > 2}, the function generates a multi-panel correlation scoreplot matrix for the comps of interest - might be slow if the there are many comps.
#' @import ggplot2
#' @import ggrepel
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @examples
#' \dontrun{
#' rbioClass_plsda_scoreplot(new_model, comps = c(1, 2, 3), plot.ellipse = TRUE)
#' }
#' @export
rbioClass_plsda_scoreplot <- function(object, y = NULL, comps = c(1, 2),
                                   plot.rightsideY = FALSE,
                                   plot.Title = NULL,
                                   plot.sampleLabel.type = c("none", "direct", "indirect"), plot.sampleLabel.vector = NULL,
                                   plot.sampleLabelSize = 2, plot.sampleLabel.padding = 0.5,
                                   plot.SymbolSize = 2,
                                   plot.ellipse = FALSE, plot.ellipse_conf = 0.95,
                                   plot.mtx.densityplot = FALSE, plot.mtx.stripLblSize = 10,
                                   plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                                   plot.fontType = "sans",
                                   plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                   plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                   plot.legendSize = 9,
                                   plot.Width = 170, plot.Height = 150,
                                   verbose = TRUE){
  ## check variables
  if (any(class(object) == "rbiomvr")){
    y <- object$inputY
  } else if (!any(class(object) == "rbiomvr") & any(class(object) == "mvr")){
    if (is.null(y)) stop("for mvr class objects, please provide response factor vector y.")
    if (class(y) != "factor"){
      if (verbose) cat("y is converted to factor. \n")
      y <- factor(y, levels = unique(y))
    } else {
      y <- y
    }
  } else stop("object needs to be \"rbiomvr\" or \"mvr\" class, e.g. created from rbioClass_plsda() function.")
  if (length(comps) > object$ncomp)stop("comps length exceeded the maximum comp length.")
  if (!all(comps %in% seq(object$ncomp)))stop("comps contain non-existant comp.")
  # if (!tolower(plot.sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\".")
  plot.sampleLabel.type <- match.arg(tolower(plot.sampleLabel.type), c("none", "direct", "indirect"))
  if (tolower(plot.sampleLabel.type != "none")){
    if (!is.null(plot.sampleLabel.vector) & length(plot.sampleLabel.vector) != nrow(object$inputX)) {
      stop("plot.sampleLabel.vector has to be the same length as number of samples in the training set from object, i.e. nrow(object$inputX).")}
  }
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")

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

    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggplot(score_x, aes(x = sample, y = axis1)) +
      geom_line(aes(colour = group, linetype = group))

    if (plot.sampleLabel.type != "none"){  # labels
      if (tolower(plot.sampleLabel.type) == "direct"){
        scoreplt <- scoreplt +
          geom_text(aes(label = sample.label, colour = group), size = plot.SymbolSize)
      } else if (tolower(plot.sampleLabel.type) == "indirect") {
        scoreplt <- scoreplt +
          geom_point(size = plot.SymbolSize, aes(colour = group, shape = group)) +
          scale_shape_manual(values=1:nlevels(score_x$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          size = plot.sampleLabelSize, show.legend = FALSE)
      }
    } else {
      scoreplt <- scoreplt +
        geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) + # plot the sample score scatter plot
        scale_shape_manual(values=1:nlevels(score_x$group))
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

    # grid.newpage()
    if (plot.rightsideY){ # add the right-side y axis
      # extract gtable
      pltgtb <- rightside_y(scoreplt)
    } else { # no right side y-axis
      pltgtb <- scoreplt
    }

  } else if (length(comps) == 2){  # two components plot
    names(score_x)[1:2] <- c("axis1", "axis2")

    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
    scoreplt <- ggplot(score_x, aes(x = axis1, y = axis2))

    if (plot.sampleLabel.type != "none"){  # labels
      if (tolower(plot.sampleLabel.type) == "direct"){
        scoreplt <- scoreplt +
          geom_text(aes(label = sample.label, colour = group), size = plot.SymbolSize)
      } else if (tolower(plot.sampleLabel.type) == "indirect") {
        scoreplt <- scoreplt +
          geom_point(size = plot.SymbolSize, aes(colour = group, shape = group)) +
          scale_shape_manual(values=1:nlevels(score_x$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          size = plot.sampleLabelSize, show.legend = FALSE)
      }
    } else {
      scoreplt <- scoreplt +
        geom_point(aes(shape = group, colour = group), size = plot.SymbolSize) + # plot the sample score scatter plot
        scale_shape_manual(values=1:nlevels(score_x$group))
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

    # grid.newpage()
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
          scale_shape_manual(values=1:nlevels(data$group)) +
          geom_text_repel(aes(label = sample.label), point.padding = unit(plot.sampleLabel.padding, "lines"),
                          size = plot.sampleLabelSize, show.legend = FALSE)
      } else {
        g <- g + geom_point(...) +
          scale_shape_manual(values=1:nlevels(data$group))
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
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.scoreplot.pdf...", sep = ""))  # initial message
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

    # grid.newpage()
    pltgtb <- scoreplt
  }
  ggsave(filename = paste(deparse(substitute(object)),".plsda.scoreplot.pdf", sep = ""), plot = pltgtb,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  if (verbose) cat("Done!\n") # final message
  grid.draw(pltgtb)

  ## return invinsible object
  invisible(pltgtb)
}


#' @title rbioClass_plsda_jackknife
#'
#' @description Jack-Knife procedure for the \code{PLS} models, e.g. \code{PLS-DA} or \code{PLS-R}.
#' @param object A \code{mvr} or \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
#' @param ncomp Defaults is all the components the \code{mvr} object has.
#' @param use.mean Defaults is \code{FALSE}.
#' @param sig.p Alpha value for the jack-knife coffecient p values. Defaults is \code{0.05}.
#' @param plot.title Whether to display plot title on top of the plot. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.outlineCol The outline colour for the bar gars. Default is \code{"black"}.
#' @param plot.errorbar Set the type of errorbar. Options are standard error of the mean (\code{"sem"}), or standard deviation (\code{"sd"}), case insensitive. Default is \code{"sem"}.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#' rbioClass_plsda_jackknife(object = new_model_optm, use.mean = FALSE,
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
rbioClass_plsda_jackknife <- function(object, ncomp = object$ncomp, use.mean = FALSE, sig.p = 0.05,
                                   plot = TRUE,
                                   plot.title = FALSE, plot.titleSize = 10,
                                   plot.outlineCol = "black",
                                   plot.errorbar = c("sem", "sd"), plot.errorbarWidth = 0.2, plot.errorbarLblSize = 6,
                                   plot.fontType = "sans",
                                   plot.xLabel = NULL, plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
                                   plot.xTickBold = FALSE, plot.xAngle = 0,
                                   plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                                   plot.rightsideY = TRUE,
                                   plot.yLabel = NULL, plot.yLabelSize = 10, plot.yTickLblSize = 10, plot.yTickItalic = FALSE, plot.yTickBold = FALSE,
                                   plot.Width = 170, plot.Height = 150,
                                   verbose = TRUE){
  # check arguments
  if (!any(class(object) %in% c("mvr", "rbiomvr"))) stop("object has to be a mvr or rbiomvr class.")
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")
  plot.errorbar <- match.arg(tolower(plot.errorbar), c("sem", "sd"))

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
      if (verbose) cat(paste("Plot saved to file: ", deparse(substitute(object)), ".", names(plot_list)[i], ".jackknife.pdf...", sep = "")) # initial message
      DfPlt <- plot_list[[i]]

      if (tolower(plot.errorbar) %in% c("sem")){  # error bar
        err <- DfPlt$sem
      } else if (tolower(plot.errorbar) %in% c("sd")){
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
        scale_x_discrete(expand = c(0.05, 0.05)) +
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
      # grid.newpage()
      if (plot.rightsideY){ # add the right-side y axis
        pltgtb <- rightside_y(plt)
      } else { # no right side y-axis
        pltgtb <- plt
      }

      ## export the file and draw a preview
      ggsave(filename = paste(deparse(substitute(object)), ".", names(plot_list)[i], ".jackknife.pdf", sep = ""), plot = pltgtb,
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      if (verbose) cat("Done!\n") # final message
      grid.draw(pltgtb) # preview
    }
  }
  cat(paste0("jackknife is built on ", ncomp, " components.\n"))
  assign(paste(deparse(substitute(object)), "_plsda_jackknife_summary_list", sep = ""), plot_list, envir = .GlobalEnv)
}


#' @title rbioFS_plsda_vip
#'
#' @description VIP, or variable importance in projection, calcualtion and plotting for plsda models. This is another FS method, and can be used independently.
#' @param object A \code{mvr} or \code{rbiomvr} object. Make sure the model is built uisng \code{"oscorespls"} method.
#' @param vip.alpha Alpha value (threshold) for VIP values. Any VIP above this is considered important. Defaults is \code{1}.
#' @param comps Integer vector. Components to plot. The index of the components are intergers. The vector length should be between 1 and the total number of components, inclusive. Default is \code{c(1, 2)}.
#' @param bootstrap If to use boostrap for VIP calculation, so that standard deviation on VIP can be estimated. Default is \code{TRUE}.
#' @param boot.n Set only when \code{boostrap = TRUE}, nummbers of iterations for boostrap. Default is \code{50}.
#' @param boot.parallelComputing Set only when \code{boostrap = TRUE}, if to use parallel computering for bootstrap process. Default is \code{TRUE}.
#' @param boot.n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param boot.clusterType Set only when \code{boostrap = TRUE} and \code{boot.parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param plot If to generate a plot. Default is \code{TRUE}.
#' @param ... Additional arguments to \code{\link{rbioFS_plsda_vip_plot}}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Outputs a \code{rbiomvr_vip} class object to the environment.
#'
#' \code{rbiomvr_vip} items are:
#'
#' \code{vip_summary}
#'
#' \code{comps}: number of PLS-DA components calculated fro VIP
#'
#' \code{vip.alpha}: cutoff VIP value for a feature to be considered as important
#'
#' \code{features_above_alpha}: selected features according to vip.alpha, i.e. features with VIP >= vip.alpha
#'
#' \code{boostrap}
#'
#' \code{boot.n}: boostrap iterations
#'
#' \code{bootstrap.iteration.results}: raw VIP values for each bootstrap iteration
#'
#' @details Only works when the plsda model is fitted with the orthorgonal score algorithm, or NIPALS. Such model can be built using \code{\link{rbioClass_plsda}} with \code{method = "oscorespls"}.
#'
#' The \code{vip.alpha} of 1 is the most commonly accepted value. However it is also acceptable to set according to the data and objectives of the study.
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' rbioFS_plsda_vip(object = new_model_optm, boostrap = TRUE, boot.m = 50,
#'                  vip.alpha = 0.8,
#'                  plot = TRUE, plot.title = TRUE, plot.titleSize = 10,
#'                  plot.sig.line = TRUE,
#'                  plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
#'                  plot.errorbarLblSize = 6, plot.fontType = "sans", plot.xLabel = "Features",
#'                  plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
#'                  plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 1, plot.xvAligh = 0.2,
#'                  plot.rightsideY = TRUE,
#'                  plot.yLabel = "Coefficients", plot.yLabelSize = 10, plot.yTickLblSize = 10,
#'                  plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
#'                  plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
#'                  plot.Height = 150)
#' }
#' @export
rbioFS_plsda_vip <- function(object, vip.alpha = 1, comps = c(1, 2),
                             bootstrap = TRUE,
                             boot.n = 50, boot.parallelComputing = TRUE, boot.n_cores = parallel::detectCores() - 1, boot.clusterType = "PSOCK",
                             plot = TRUE,...,
                             verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c("rbiomvr", "mvr"))) stop("object needs to be either a \"rbiomvr\" or \"mvr\" class.")
  if (object$method != "oscorespls") stop("Object needs to fit using oscorespls (i.e. NIPALS) algorithm. Please re-fit using rbioClass_plsda with method = \"oscorespls\".") # only oscorespls algorithm is applicable to VIP
  if (length(comps) > object$ncomp)stop("comps length exceeded the maximum comp length.")
  if (!all(comps %in% seq(object$ncomp)))stop("comps contain non-existant comp.")
  if (bootstrap & (boot.n < 2 | boot.n %% 1 != 0)) stop("when boostrap = TRUE, boot.n needs to be an integer greater than 1.")
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")

  ## VIP process
  if (verbose) cat(paste0("Boostrap: ", ifelse(bootstrap, "ON\n", "OFF\n")))
  if (bootstrap){  # bootstrap VIP process
    # message
    if (verbose) cat(paste0("Parallel computing: ", ifelse(boot.parallelComputing, "ON\n", "OFF\n")))
    if (verbose) cat(paste0("Boostrap interation: ", boot.n, "\n"))
    if (verbose) cat("Boostrap VIP calculation (speed depending on hardware configuration)...")

    # set up data
    orig.dat <- data.frame(Y = object$inputY, object$inputX, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    sum.boot.vip_raw_list <- vector(mode = "list", length = boot.n)

    # bootstrap modelling and VIP calculation
    if (!boot.parallelComputing){  # non-parallel
      sum.boot.vip_raw_list[] <- foreach(i = 1:boot.n) %do% {
        # bootstrap resampling
        boot.idx <- foreach(j = levels(object$inputY), .combine = "c") %do% {
          group.idx <- which(object$inputY == j)
          group.boot.idx <- sample(group.idx, replace = TRUE)
          group.boot.idx
        }
        boot.dat <- orig.dat[boot.idx, ]

        # bootstrap modelling
        boot.X <- boot.dat[, -1]
        boot.Y <- boot.dat[, 1]
        boot.m <- rbioClass_plsda(x = boot.X, y = factor(boot.Y, levels = unique(boot.Y)),
                                  ncomp = object$ncomp, validation = "CV", method = "oscorespls",
                                  verbose = FALSE)

        # boostrap VIP
        boot.score <- boot.m$scores
        boot.lodw <- boot.m$loading.weights
        boot.Wnorm2 <- colSums(boot.lodw^2)

        # vip_raw_list contains VIP for all the components
        boot.vip_raw_list <- vector(mode = "list", length = dim(boot.m$Yloadings)[1])
        boot.vip_raw_list[] <- foreach(m = 1:dim(boot.m$Yloadings)[1]) %do% {
          ylod <- boot.m$Yloadings[m, ] # one y variable at a time
          boot.ss <- c(ylod)^2 * colSums(boot.score^2)
          boot.ssw <- sweep(boot.lodw^2, 2, boot.ss / boot.Wnorm2, "*")
          boot.vip <- sqrt(nrow(boot.ssw) * apply(boot.ssw, 1, cumsum) / cumsum(boot.ss))  # cumsum: cucmulative sum
          if (object$ncomp == 1) {
            boot.vip <- matrix(boot.vip, nrow = 1, dimnames = list('Comp 1', names(boot.vip)))
          }
          boot.vip <- boot.vip[comps, , drop = FALSE]
          return(boot.vip)
        }
        names(boot.vip_raw_list) <- dimnames(boot.m$Yloadings)[[1]]
        boot.vip_raw_list
      }
    } else {  # boot parallel computing
      # set up cluster
      n_cores <- boot.n_cores
      cl <- makeCluster(n_cores, type = boot.clusterType)
      registerDoParallel(cl)
      on.exit(stopCluster(cl)) # close connect when exiting the function

      # computing
      sum.boot.vip_raw_list[] <- foreach(i = 1:boot.n, .packages = c("foreach", "RBioFS")) %dopar% {
        # bootstrap resampling
        boot.idx <- foreach(j = levels(object$inputY), .combine = "c") %do% {
          group.idx <- which(object$inputY == j)
          group.boot.idx <- sample(group.idx, replace = TRUE)
          group.boot.idx
        }
        boot.dat <- orig.dat[boot.idx, ]

        # bootstrap modelling
        boot.X <- boot.dat[, -1]
        boot.Y <- boot.dat[, 1]
        boot.m <- rbioClass_plsda(x = boot.X, y = factor(boot.Y, levels = unique(boot.Y)),
                                  ncomp = object$ncomp, validation = "CV", method = "oscorespls",
                                  verbose = FALSE)

        # boostrap VIP
        boot.score <- boot.m$scores
        boot.lodw <- boot.m$loading.weights
        boot.Wnorm2 <- colSums(boot.lodw^2)

        # vip_raw_list contains VIP for all the components
        boot.vip_raw_list <- vector(mode = "list", length = dim(boot.m$Yloadings)[1])
        boot.vip_raw_list[] <- foreach(m = 1:dim(boot.m$Yloadings)[1]) %do% {
          ylod <- boot.m$Yloadings[m, ] # one y variable at a time
          boot.ss <- c(ylod)^2 * colSums(boot.score^2)
          boot.ssw <- sweep(boot.lodw^2, 2, boot.ss / boot.Wnorm2, "*")
          boot.vip <- sqrt(nrow(boot.ssw) * apply(boot.ssw, 1, cumsum) / cumsum(boot.ss))  # cumsum: cucmulative sum
          if (object$ncomp == 1) {
            boot.vip <- matrix(boot.vip, nrow = 1, dimnames = list('Comp 1', names(boot.vip)))
          }
          boot.vip <- boot.vip[comps, , drop = FALSE]
          return(boot.vip)
        }
        names(boot.vip_raw_list) <- dimnames(boot.m$Yloadings)[[1]]
        boot.vip_raw_list
      }
    }
    names(sum.boot.vip_raw_list) <- paste0("boot ", seq(boot.n))

    # sum.boot.vip_raw_list structure: boot$group$: (row)`comp 1`|`comp 2`|`comp 3`...
    # targeted format: group$`comp 1`: (colums)`boot 1`|`boot 2`|`boot 3`...
    group.comp.boot.vip_list <- vector(mode = "list", length = length(levels(object$inputY)))
    group.comp.boot.vip_list[] <- foreach(i = 1:length(levels(object$inputY))) %do% {
      group.comp.list <- vector(mode = "list", length = length(comps))
      group.comp.list[] <- foreach(j = comps) %do% {
        group.comp.boot_mtx <- foreach(m = 1:boot.n, .combine = "cbind") %do% {
          sum.boot.vip_raw_list[[m]][[i]][j, ]  # i: groups, j: boot, 1: comp
        }
        colnames(group.comp.boot_mtx) <- paste0("boot ", seq(boot.n))
        group.comp.boot_mtx
      }
      names(group.comp.list) <- paste0("comp ", comps)
      group.comp.list
    }
    names(group.comp.boot.vip_list) <- levels(object$inputY)

    # plot list
    boot.vip.plt_dat_list <- vector(mode = "list", length = length(levels(object$inputY)))
    boot.vip.plt_dat_list[] <- foreach(i = 1:length(levels(object$inputY))) %do% {
      boot.plt_list.comp <- vector(mode = "list", length = length(comps))
      boot.plt_list.comp[] <- foreach(j = comps) %do% {
        boot.mean <- rowMeans(group.comp.boot.vip_list[[i]][[j]])
        boot.sd <- matrixStats::rowSds(group.comp.boot.vip_list[[i]][[j]])
        boot.sem <- matrixStats::rowSds(group.comp.boot.vip_list[[i]][[j]])/sqrt(boot.n)
        outdfm <- data.frame(Features = rownames(group.comp.boot.vip_list[[i]][[j]]),
                             VIP = boot.mean,
                             sd = boot.sd,
                             sem = boot.sem,
                             row.names = NULL,
                             check.names = FALSE)
        outdfm <- outdfm[order(outdfm$VIP, decreasing = TRUE),]
        outdfm$Features <- factor(outdfm$Features, levels = unique(outdfm$Features))
        rownames(outdfm) <- NULL
        outdfm
      }
      names(boot.plt_list.comp) <- paste0("comp ", comps)
      boot.plt_list.comp
    }
    names(boot.vip.plt_dat_list) <- levels(object$inputY)
    vip_list <- boot.vip.plt_dat_list
  } else {  # non bootstrap VIP process
    # message
    if (verbose) cat("VIP calculation...")

    # VIP calculation
    score <- object$scores
    lodw <- object$loading.weights
    Wnorm2 <- colSums(lodw^2)

    # vip_raw_list contains VIP for all the components
    vip_raw_list <- vector(mode = "list", length = dim(object$Yloadings)[1])
    vip_raw_list[] <- foreach(i = 1:dim(object$Yloadings)[1]) %do% {
      ylod <- object$Yloadings[i, ] # one y variable at a time
      ss <- c(ylod)^2 * colSums(score^2)
      ssw <- sweep(lodw^2, 2, ss / Wnorm2, "*")
      vip <- sqrt(nrow(ssw) * apply(ssw, 1, cumsum) / cumsum(ss))  # cumsum: cucmulative sum
      vip <- vip[comps, , drop = FALSE]
      vip
    }
    names(vip_raw_list) <- dimnames(object$Yloadings)[[1]]

    vip_list <- vector(mode = "list", length = length(levels(object$inputY)))
    vip_list[] <- foreach(i = 1:length(levels(object$inputY))) %do% {
      plt_list.comp <- vector(mode = "list", length = length(comps))
      plt_list.comp[] <- foreach(j = comps) %do% {
        vip <- vip_raw_list[[i]][j, ]
        outdfm <- data.frame(Features = colnames(vip_raw_list[[i]]),
                             VIP = vip,
                             row.names = NULL,
                             check.names = FALSE)
        outdfm <- outdfm[order(outdfm$VIP, decreasing = TRUE),]
        outdfm$Features <- factor(outdfm$Features, levels = unique(outdfm$Features))
        rownames(outdfm) <- NULL
        outdfm
      }
      names(plt_list.comp) <- paste0("comp ", comps)
      plt_list.comp
    }
    names(vip_list) <- levels(object$inputY)
  }
  if (verbose) cat("Done!\n")

  ## output
  # important features lsit
  final.ipf_list <- vector(mode = "list", length = length(vip_list))
  final.ipf_list <- foreach(m = 1:length(vip_list)) %do% {
    ipf_list <- vector(mode = "list", length = length(vip_list[[m]]))
    ipf_list[] <- foreach(n = 1:length(vip_list[[m]])) %do% {
      as.character(vip_list[[m]][[n]][which(vip_list[[m]][[n]]$VIP > vip.alpha), 1])
    }
    names(ipf_list) <- names(vip_list[[m]])
    ipf_list
  }
  names(final.ipf_list) <- names(vip_list)

  # output
  out <- list(vip_summary = vip_list,
              vip.alpha = vip.alpha,
              features_above_alpha = final.ipf_list,
              comps = comps,
              bootstrap = bootstrap,
              boot.n = if (bootstrap) boot.n else NULL,
              bootstrap.iteration.results = if (bootstrap) group.comp.boot.vip_list else NULL,
              model.type = object$model.type)
  class(out) <- "rbiomvr_vip"
  assign(paste(deparse(substitute(object)), "_plsda_vip", sep = ""), out, envir = .GlobalEnv)

  ## plot
  if (plot){
    # message
    if (verbose) cat("\n")
    RBioFS::rbioFS_plsda_vip_plot(vip_obj = out, ...)
  }
}



#' @title rbioFS_plsda_vip_plot
#'
#' @description Plotting function for \code{\link{rbioFS_plsda_VIP}}.
#' @param vip_obj A \code{rbiomvr_vip}, generated by \code{\link{rbioFS_plsda_VIP}}.
#' @param plot.preview If to preview plots. Default is \code{TRUE}.
#' @param plot.title Whether to display plot title on top of the plot. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.sig.line Wether to display a horizontal line indicating the VIP threshold. Default is \code{TRUE}.
#' @param plot.outlineCol The outline colour for the bar gars. Default is \code{"black"}.
#' @param plot.errorbar Set the type of errorbar. Only applicable if the object model is built with bootstrapping. Options are standard error of the mean (\code{"SEM"}, \code{"standard error"}, \code{"standard error of the mean"}), or standard deviation (\code{"SD"}, \code{"standard deviation"}), case insensitive. Default is \code{"SEM"}.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Outputs pdf figure files to the working directory.
#' @details Only works with \code{rbiomvr_vip} objects. Set \code{plot.preview = FALSE} to speed up the process as preview rendering may be slow.
#'          The function also applies to the regression model (model.type = "regression").
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y
#' @importFrom scales rescale_none
#' @examples
#' \dontrun{
#' rbioFS_plsda_vip_plot(vip_obj = plsda_m_vip,
#'                       plot.title = TRUE, plot.titleSize = 10,
#'                       plot.sig.line = TRUE,
#'                       plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
#'                       plot.errorbarLblSize = 6, plot.fontType = "sans", plot.xLabel = "Features",
#'                       plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
#'                       plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 1, plot.xvAligh = 0.2, plot.rightsideY = TRUE,
#'                       plot.yLabel = "Coefficients", plot.yLabelSize = 10, plot.yTickLblSize = 10,
#'                       plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
#'                       plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
#'                       plot.Height = 150)
#' }
#' @export
rbioFS_plsda_vip_plot <- function(vip_obj, plot.preview = TRUE,
                                  plot.title = TRUE, plot.titleSize = 10,
                                  plot.sig.line = TRUE,
                                  plot.outlineCol = "black", plot.errorbar = "SEM", plot.errorbarWidth = 0.2,
                                  plot.errorbarLblSize = 6, plot.fontType = "sans",
                                  plot.xLabel = "Features", plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
                                  plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 0.95, plot.xvAlign = 0.5,
                                  plot.rightsideY = TRUE,
                                  plot.yLabel = "VIP", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                  plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
                                  plot.legendTtl = FALSE, plot.legendTtlSize = 9, plot.Width = 170,
                                  plot.Height = 150, verbose = TRUE){
  ## argument check
  if (!any(class(vip_obj) %in% c("rbiomvr_vip"))) stop("object needs to be a \"rbiomvr_vip\" class.")
  # if (vip_obj$model.type != "classification") stop("object needs to have model.type = \"classification\"")

  ## plot
  # set up bootstrap setting
  bootstrap <- vip_obj$bootstrap

  # set up data
  vip_list <- vip_obj$vip_summary

  # plot
  # grid.newpage()
  final.plt_list <- vector(mode = "list", length = length(vip_list))
  final.plt_list <- foreach(i = 1:length(vip_list)) %do% {
    plt_list <- foreach(j = 1:length(vip_list[[i]])) %do% {
      if (bootstrap){  #  detect if to use errorbar
        if (tolower(plot.errorbar) %in% c("sem", "standard error", "standard error of the mean")){  # error bar
          err <- vip_list[[i]][[j]]$sem
        } else if (tolower(plot.errorbar) %in% c("sd", "standard deviation")){
          err <- vip_list[[i]][[j]]$sd
        }
        y_axis_Mx <- max((vip_list[[i]][[j]]$VIP + err) * 1.05) * 1.2

      } else {
        y_axis_Mx <- max((vip_list[[i]][[j]]$VIP) * 1.05) * 1.2
      }
      y_axis_Mn <- 0

      baseplt <- ggplot(data = vip_list[[i]][[j]], aes(x = Features, y = VIP)) +
        geom_bar(position = "dodge", stat = "identity", color = plot.outlineCol) +
        # scale_x_discrete(expand = c(0.05, 0.05)) +
        scale_y_continuous(expand = c(0, 0), limits = c(y_axis_Mn, y_axis_Mx),
                           oob = rescale_none) +
        xlab(plot.xLabel) +
        ylab(paste0(plot.yLabel, " (", names(vip_list[[i]])[j], ")")) +
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

      if (bootstrap){
        plt <- baseplt +
          geom_errorbar(aes(ymin = VIP - err, ymax = VIP + err),
                        position = position_dodge(0.9), color = "black", width = plot.errorbarWidth)
      } else {
        plt <- baseplt
      }

      if (plot.title){
        plt <- plt + ggtitle(names(vip_list)[i])
      } else (
        plt <- plt + ggtitle(NULL)
      )

      if (plot.sig.line){
        plt <- plt +
          geom_hline(yintercept = vip_obj$vip.alpha, linetype = "dashed", colour = "red")
      }

      if (plot.xTickItalic & plot.xTickBold){
        plt <- plt +
          theme(axis.text.x = element_text(face = "bold.italic"))
      } else if (plot.xTickItalic & !plot.xTickBold){
        plt <- plt +
          theme(axis.text.x = element_text(face = "italic"))
      } else if (plot.xTickBold & !plot.xTickItalic){
        plt <- plt +
          theme(axis.text.x = element_text(face = "bold"))
      }

      if (plot.yTickItalic & plot.yTickBold){
        plt <- plt +
          theme(axis.text.y  = element_text(face = "bold.italic"))
      } else if (plot.yTickItalic & !plot.yTickBold){
        plt <- plt +
          theme(axis.text.y = element_text(face = "italic"))
      } else if (plot.yTickBold & !plot.yTickItalic){
        plt <- plt +
          theme(axis.text.y = element_text(face = "bold"))
      }

      plt <- RBioplot::rightside_y(plt)

      ## export the file and draw a preview
      if (verbose) cat(paste0("Plot saved to file: ", deparse(substitute(vip_obj)), ".", names(vip_list)[i], ".", names(vip_list[[i]])[j], ".vip.pdf..."))
      ggsave(filename = paste(deparse(substitute(vip_obj)), ".", names(vip_list)[i], ".", names(vip_list[[i]])[j], ".vip.pdf", sep = ""), plot = plt,
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      if (verbose) cat("Done!\n") # final message
      if (plot.preview) grid.draw(plt) # preview
    }
    plt_list
  }
  names(final.plt_list) <- names(vip_list)

  ## not-export
  invisible(final.plt_list)
}


#' @title rbioClass_plsda_roc_auc()
#'
#' @description ROC-AUC analysis and ploting for plsda model
#' @param object A \code{rbiomvr} object. Make sure the object is generated with a \code{validation} section.
#' @param newdata Newdata (test data) for ROC-AUC analysis, excluding labels, i.e. y. If missing, the function will use the transformed data from the model object.
#' @param newdata.label Newdata label vector (i.e. test data y). If missing, the function will use the training data and its corresponding labels.
#' @param center.newdata Only set when both \code{newdata} and \code{newdata.label} are set, if to center.scale newdata. Default is \code{TRUE}.
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
#' @param plot.lineSize Line size. Default is \code{1}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"1 - specificity"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL. Default is \code{"sensitivity"}.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width ROC width. Default is \code{170}.
#' @param plot.Height ROC height. Default is \code{150}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#' rbioClass_plsda_roc_auc(object = model_binary, rocplot = TRUE, plot.comps = 1:2)
#' }
#' @export
rbioClass_plsda_roc_auc <- function(object, newdata, newdata.label, center.newdata = TRUE,
                                    rocplot = TRUE,
                                    plot.comps = 1:object$ncomp,
                                    plot.smooth = FALSE,
                                    multi_plot.ncol = length(plot.comps), multi_plot.nrow = 1, multi_plot.legend.pos = c("bottom", "top", "left", "right"),
                                    plot.rightsideY = TRUE,
                                    plot.SymbolSize = 2, plot.lineSize = 1,
                                    plot.display.Title = TRUE, plot.titleSize = 10,
                                    plot.fontType = "sans",
                                    plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                    plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                    plot.legendSize = 9,
                                    plot.Width = 170, plot.Height = 150,
                                    verbose = TRUE){
  ## check arguments
  if (missing(newdata) || is.null(newdata) || missing(newdata.label) || is.null(newdata.label)) {
    cat("Note: newdata or newdata.label info isn't complete. Proceed with training data.\n")
    newdata <- object$centerX$centerX
    outcome <- object$inputY
  } else {
    if (class(newdata.label) != "factor"){
      if (verbose) cat("y is converted to factor. \n")
      newdata.label <- factor(newdata.label, levels = unique(newdata.label))
    }
    ## center data with the option of scaling
    if (center.newdata){
      if (verbose) cat("Data center.scaled using training data column mean and sd, prior to modelling.\n")
      centerdata <- t((t(newdata) - object$centerX$meanX) / object$centerX$columnSD)
      newdata <- centerdata
    }
    outcome <- newdata.label
  }
  if (!any(class(object) %in% c("rbiomvr"))) stop("object needs to be either a \"rbiomvr\" class.")
  if(plot.smooth) cat("ROC smooth: ON.\n") else cat("ROC smooth: OFF.\n")
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")
  multi_plot.legend.pos <- match.arg(tolower(multi_plot.legend.pos), c("bottom", "top", "left", "right"))

  ## calcuate ROC-AUC
  pred_raw <- predict(object, newdata = newdata)

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
      perf <- tryCatch(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth, ci = TRUE),
                       error = function(err){
                         cat("Curve not smoothable. Proceed without smooth.\n")
                         pROC::roc(controls = controls, cases = cases, smooth = FALSE, ci = TRUE)
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
      thresholds <- perf$thresholds
      mtx <- cbind(fpr, tpr, thresholds)
      if (length(levels(outcome)) == 2){
        df <- data.frame(mtx, group = rep(levels(outcome)[j], times = nrow(mtx)), row.names = NULL, check.names = FALSE)
      } else {
        df <- data.frame(mtx, group = rep(paste0(levels(outcome)[j], " (vs Others)"), times = nrow(mtx)), row.names = NULL, check.names = FALSE)
      }
      df <- df[order(df$tpr), ]  # order by tpr so that all the points will be connected on graph
      df
    }
    out
  }
  names(roc_dfm_list) <- paste0("comp ", 1:length(plot.comps))

  ## plot
  if (rocplot){
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.roc.pdf...", sep = ""))  # initial message

    # extract comp variance distribution
    varpp_x <- 100 * object$Xvar / object$Xtotvar
    boxdfm_x <- data.frame(comp_x = as.numeric(gsub("Comp", "", names(varpp_x))), varpp_x = varpp_x)
    var_percentage_x <- varpp_x[paste0("Comp ", plot.comps)] # extract the proportion of variance for the selected PCs
    comp_axis_lbl <- paste("comp ", plot.comps, " (", round(var_percentage_x, digits = 2), "%)", sep = "")

    plt_list <- vector(mode = "list", length = length(plot.comps))
    plt_list[] <- foreach(k = 1:length(plot.comps)) %do% {
      plt <- ggplot(data = roc_dfm_list[[k]], aes(x = fpr, y = tpr, group = group, colour = group)) +
        geom_line(aes(linetype = group), size = plot.lineSize) +
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
    # grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.roc.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    if (verbose) cat("Done!\n")
  }
  assign(paste(deparse(substitute(object)), "_plsda_roc_list", sep = ""), roc_dfm_list, envir = .GlobalEnv)
}


#' @title rbioClass_plsda_predict
#'
#' @description Prediction function for PLS-DA analysis. The function calculates the predicted value for unknown sample data using the input PLS-DA model.
#' @param object A \code{rbiomvr} object.
#' @param comps  Number of PLS-DA components used in the model. Default is \code{object$ncomp}.
#' @param newdata Input data to be classified. Make sure it is a \code{matrix} class and has the same variables as the model, i.e. same number of columns as the training data.
#' @param center.newdata If to center the newdata. When \code{TRUE}, it will also apply the same scaling option as the \code{object}. Default is \code{TRUE}.
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
#' @param plot.Width Predplot width. Default is \code{170}.
#' @param plot.Height Predplot height. Default is \code{150}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return  A \code{prediction} obejct, as well as pdf figure file for predicted values if \code{predplot = TRUE}.
#'
#' The items of the object are:
#'
#' \code{model.type}
#'
#' \code{classifier.class}
#'
#' \code{predited.value}
#'
#' \code{prob.method}  Method to caluculate posterier probability
#'
#' \code{probability.summary}
#'
#' \code{raw.newdata}
#'
#' \code{center.scale}
#'
#' \code{center.scaled.newdata}
#'
#' \code{newdata.y}
#'
#' @details Although optional, the \code{newdata} matrix should be centered prior to testing, with the same scaling setting as the input \code{rbiomvr} object. The option \code{center.newdata = FALSE} is
#' for the already centered the data matrix. This center.scale process should use training data's column mean and column standard deviation.
#'
#' Regarding \code{threshold}, the value should between \code{0} and \code{1}. It's the flank region around the dummified classification values \code{0} (i.e. "control") and \code{1} (i.e. "case").
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
#' For classification plot, the output \code{prediction} object should be used with function \code{\link{rbioClass_plsda_classplot()}}.
#'
#' If \code{sampleLabel.vector = NULL} or missing, the function uses row numbers as label.
#'
#' @import ggplot2
#' @import pls
#' @import ggrepel
#' @importFrom klaR NaiveBayes
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioClass_plsda_predict(object = new_model_optm, newdata = newdata, prob.method = "Bayes",
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
rbioClass_plsda_predict <- function(object, comps = object$ncomp, newdata, center.newdata = TRUE,
                                    prob.method = "Bayes",
                                    threshold = 0.2,
                                    predplot = TRUE,
                                    plot.sampleLabel.type = "none", plot.sampleLabel.vector = NULL,
                                    plot.sampleLabelSize = 2, plot.sampleLabel.padding = 0.5,
                                    multi_plot.ncol = length(levels(object$inputY)), multi_plot.nrow = 1, multi_plot.legend.pos = c("bottom", "top", "left", "right"),
                                    plot.rightsideY = TRUE,
                                    plot.SymbolSize = 2, plot.display.Title = TRUE, plot.titleSize = 10,
                                    plot.fontType = "sans", plot.unclassifiedColour = "gray", plot.classifiedColour = "red",
                                    plot.xLabel = "Samples", plot.xLabelSize = 10, plot.xTickLblSize = 6, plot.xTickItalic = FALSE,
                                    plot.xAngle = 0, plot.xhAlign = 0.5, plot.xvAlign = 0.5,
                                    plot.yLabel = "Predicted values", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                    plot.legendSize = 9,
                                    plot.Width = 170, plot.Height = 150,
                                    verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c("rbiomvr"))) stop("object needs to be a \"rbiomvr\" class.")
  if (object$model.type != "classification") stop("object needs to have model.type = \"classification\"")
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x)))stop("x needs to be a matrix, data.frame or vector.")
  if (class(newdata) == "data.frame" | is.null(dim(newdata))){
    if (verbose) cat("newdata converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (ncol(newdata) != ncol(object$inputX)) stop("newdata needs to have the same number of variables, i.e. columns, as the object.")
  if (!prob.method %in% c("softmax", "Bayes")) stop("Probability method should be either \"softmax\" or \"Bayes\".")
  # if (!tolower(plot.sampleLabel.type) %in% c("none", "direct", "indirect")) stop("sampleLabel.type argument has to be one of \"none\", \"direct\" or \"indirect\".")
  multi_plot.legend.pos <- match.arg(multi_plot.legend.pos, c("bottom", "top", "left", "right"))
  if (tolower(plot.sampleLabel.type != "none")){
    if (!is.null(plot.sampleLabel.vector) & length(plot.sampleLabel.vector) != nrow(newdata)) {
      stop("plot.sampleLabel.vector has to be the same length as nrow(newdata).")}
  }

  ## center data with the option of scaling
  if (center.newdata){
    if (verbose) cat("Data center.scaled using training data column mean and sd, prior to modelling.\n")
    centerdata <- t((t(newdata) - object$centerX$meanX) / object$centerX$columnSD)
    test <- centerdata
  } else {
    centerdata <- NULL
    test <- newdata
  }

  ## prediction and class assign
  rownames(test) <- 1:nrow(test)
  pred <- predict(object = object, ncomp = comps, newdata = test, type = "response")
  if (is.null(plot.sampleLabel.vector) | missing(plot.sampleLabel.vector)){
    sample.label <- seq(nrow(test))
  } else {
    sample.label <- plot.sampleLabel.vector
  }

  ## classification and probability calculation
  pred_mtx <- pred[, ,1]
  if (is.null(dim(pred_mtx))) {  # if only one sample
    pred_mtx <- t(as.matrix(pred_mtx))
  }
  if (!missing(plot.sampleLabel.vector) & !is.null(plot.sampleLabel.vector)){
    rownames(pred_mtx) <- plot.sampleLabel.vector
  }

  if (verbose) cat(paste("Classification probability calculation using ", prob.method, " method...", sep = ""))  # initial message
  if (prob.method == "Bayes"){ # Naive Bayes probability calculation and prediction for sample data
    training_mtx <- as.matrix(object$centerX$centerX)
    trainingpred <- predict(object = object, ncomp = comps, newdata = training_mtx, type = "response")
    bayes.prob <- klaR::NaiveBayes(x = trainingpred,
                                   grouping = object$inputY, usekernel = TRUE)
    bayes.prob$train.posterior <- predict(bayes.prob)$posterior  # calcuate posterior probability for
    bayes.prob$x <- NULL

    bayespred <- predict(object = bayes.prob, newdata = pred_mtx)
    prob <- bayespred$posterior
  } else {
    group <- colnames(pred_mtx)
    # calcuate probability
    prob_mtx <- apply(pred_mtx, 1, FUN = function(x) exp(x) / sum(exp(x)))
    rownames(prob_mtx) <- group
    prob <- t(prob_mtx)
  }

  prob_dfm <- foreach(i = 1:nrow(pred_mtx), .combine = "rbind") %do% {
    prob_dfm <- data.frame(Sample = rep(rownames(prob)[i], times = ncol(prob)),
                           Class = colnames(prob), Probability = prob[i, ], stringsAsFactors = FALSE, row.names = NULL)
    prob_dfm$Sample <- factor(prob_dfm$Sample, unique(prob_dfm$Sample))
    prob_dfm$repel.label.pos <- rev(cumsum(rev(prob_dfm$Probability)) - rev(prob_dfm$Probability) / 2)  # calculate the repel lable position, seemingly from bottom up
    prob_dfm$precent.label <- paste0(signif(prob_dfm$Probability, 4) * 100, "%")
    prob_dfm$Class <- factor(prob_dfm$Class, unique(prob_dfm$Class))
    return(prob_dfm)
  }
  if (verbose) cat(paste("Done!\n", sep = ""))  # initial message

  ## prediction plot
  if (predplot){
    # prepare prediction value list
    predlist <- vector(mode = "list", length = length(levels(object$inputY)))
    predlist[] <- foreach(i = 1:length(levels(object$inputY))) %do% {
      preddfm <- data.frame(sample = as.integer(rownames(test)), sample.label = sample.label, predicted.value = pred[, i,])
      preddfm$classification <- sapply(preddfm$predicted.value, FUN = function(x)ifelse(x > 1 - threshold & x < (1 + threshold), levels(object$inputY)[i], ifelse(x > - threshold & x < threshold, "rest", "undetermined")))
      preddfm$`Within threshold` <- ifelse(preddfm$classification == "undetermined", "N", "Y")
      preddfm$`Within threshold` <- factor(preddfm$`Within threshold`, levels = c("Y", "N"))
      return(preddfm)
    }
    names(predlist) <- levels(object$inputY)

    # plot
    if (plot.sampleLabel.type != "none" & is.null(plot.sampleLabel.vector)){  # message when no label vector is provided.
      cat("plot.sampleLabel.vector not provided. Proceed with row numbers as sampole labels.\n")
    }
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".plsda.predict.pdf...", sep = ""))  # initial message
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
    # grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".plsda.predict.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    if (verbose) cat("Done!\n")
  }

  ## export
  out <- list(model.type = "classification",
              classifier.class = class(object),
              predicted.value = pred_mtx,
              tot.predict.RMSE = NULL,
              probability.method = prob.method,
              probability.summary = prob_dfm,
              raw.newdata = newdata,
              center.scale = center.newdata,
              center.scaled.newdata = centerdata,
              newdata.y <- NULL)
  class(out) <- "prediction"
  assign(paste(deparse(substitute(object)), "_plsda_predict", sep = ""), out, envir = .GlobalEnv)
}
