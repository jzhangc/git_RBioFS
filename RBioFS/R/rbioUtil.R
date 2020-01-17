#' Title rbioUtil_classif_accuracy
#'
#' @description Plotting function for permutation test results.
#' @param object Classification model object. The object should be either \code{rbiosvm}, or \code{rbiomvr} classes.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.label Permutation test results object. The object should be either \code{rbiosvm_perm}, or \code{rbiomvr_perm} classes.
#' @param center.scale.newdata Logical, wether center and scale the newdata with training data mean and standard deviation. Default is \code{TRUE}.
#' @details accuracy = true predictions/total predictions
#' @return Numerical model accuracy value.
#' @examples
#'
#' \dontrun{
#' rbioUtil_classif_accuracy(perm_res = svm_model_perm)
#' }
#'
#' @export
rbioUtil_classif_accuracy <- function(object, newdata, newdata.label, center.scale.newdata = TRUE, verbose = TRUE) {
  # check arguments
  if (!any(class(object) %in% c('rbiosvm', "rbiomvr"))) stop("object needs to be one of \"rbiosvm\" or \"rbiomvr\" classes.")
  if (object$model.type != "classification") stop("The object should be a classification model type. ")
  if (!class(newdata) %in% c("data.frame", "matrix") & !is.null(dim(newdata))) stop("newdata needs to be a matrix, data.frame or vector.")
  if (class(newdata) == "data.frame" | is.null(dim(newdata))){
    if (verbose) cat("newdata converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (ncol(newdata) != ncol(object$inputX)) stop("test data should have the same number of variables as the training data.")


  ## process data
  if (center.scale.newdata){ # using training data mean and sd
    if (verbose) cat(paste0("Data center.scaled using training data column mean and sd, prior to modelling.\n"))
    centered_newdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
    test <- centered_newdata
  } else {
    centered_newdata <- NULL
    test <- newdata
  }

  # accuracy calculation
  newdata.label <- factor(newdata.label, levels = unique(newdata.label))
  confusion_mtx <- table(predict(object = object, newdata = newdata), newdata.label, dnn = c("Predicted", "Actual"))
  accu <- sum(diag(confusion_mtx))/sum(confusion_mtx)  # true prediction/total prediction

  return(accu)
}



#' Title rbioUtil_perm_plot
#'
#' @description Plotting function for permutation test results.
#' @param perm_res Permutation test results object. The object should be either \code{rbiosvm_perm}, or \code{rbiomvr_perm} classes.
#' @param ... Additional argument for the plot settings, see details.
#' @return A pdf file containing a scatter plot for permutation results.
#' @examples
#'
#' \dontrun{
#' rbioUtil_perm_plot(perm_res = svm_model_perm)
#' }
#'
#' @export
rbioUtil_perm_plot <- function(perm_res, ...){
  UseMethod("rbioUtil_perm_plot", perm_res)
}


#' Title rbioUtil_perm_plot.rbiosvm_perm
#'
#' @rdname rbioUtil_perm_plot
#' @method rbioUtil_perm_plot rbiosvm_perm
#' @param perm_res Permutation test results object. The object should be \code{rbiosvm_perm} class.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param ... Additional argument for the plot settings.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @export
rbioUtil_perm_plot.rbiosvm_perm <- function(perm_res, plot.SymbolSize = 2,
                                            plot.Width = 170, plot.Height = 150,
                                            ...,
                                            verbose = TRUE){
  # plot dfm
  perm_res_dfm <- perm_res$perm.results

  # plot
  if (perm_res$model.type == "regression"){
    baseplt <- ggplot(data = perm_res_dfm, aes(x = nperm, y = tot.RMSE)) +
      geom_line(linetype = "dashed") +
      geom_point(size = plot.SymbolSize) +
      geom_hline(yintercept = perm_res_dfm[1, 2], linetype = "dashed", colour = "red")

    plt <- rbioUtil_perm_plot.default(baseplt = baseplt, plot.yLabel = "tot.RMSE", ...)
  } else {
    baseplt <- ggplot(data = perm_res_dfm, aes(x = nperm, y = tot.accuracy)) +
      geom_line(linetype = "dashed") +
      geom_point(size = plot.SymbolSize) +
      geom_hline(yintercept = perm_res_dfm[1, 2], linetype = "dashed", colour = "red")

    plt <- rbioUtil_perm_plot.default(baseplt = baseplt, plot.yLabel = "tot.Accuracy", ...)
  }

  # save
  if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(perm_res)),".svm.perm.plot.pdf...", sep = ""))  # initial message
  # grid.newpage()
  ggsave(filename = paste(deparse(substitute(perm_res)),".svm.perm.plot.pdf", sep = ""), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  if(verbose) cat("Done! \n")
}


#' Title rbioUtil_perm_plot.rbiomvr_perm
#'
#' @rdname rbioUtil_perm_plot
#' @method rbioUtil_perm_plot rbiomvr_perm
#' @param perm_res Permutation test results object. The object should be \code{rbiomvr_perm} class.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param ... Additional argument for the plot settings.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @details The function automatically adjust for the regression or classification types for rbiomvr object.
#' @export
rbioUtil_perm_plot.rbiomvr_perm <- function(perm_res, plot.SymbolSize = 2,
                                            plot.Width = 170, plot.Height = 150,
                                            ...,
                                            verbose = TRUE){
  # plot dfm
  perm_res_dfm <- perm_res$perm.results

  # plot
  if (perm_res$model.type == "classification"){
    baseplt <- ggplot(data = perm_res_dfm, aes(x = nperm, y = RMSEP, group = comparison)) +
      geom_line(aes(linetype = comparison, color = comparison)) +
      geom_point(aes(shape = comparison, color = comparison), size = plot.SymbolSize) +
      geom_hline(yintercept = perm_res_dfm[perm_res_dfm$nperm == 0, 3], linetype = "dashed", colour = "red")

    plt <- rbioUtil_perm_plot.default(baseplt = baseplt, plot.yLabel = "RMSEP", ...)
  } else {
    baseplt <- ggplot(data = perm_res_dfm, aes(x = nperm, y = RMSEP)) +
      geom_line(aes(linetype = comparison)) +
      geom_point(aes(shape = comparison), size = plot.SymbolSize) +
      geom_hline(yintercept = perm_res_dfm[perm_res_dfm$nperm == 0, 3], linetype = "dashed", colour = "red")

    plt <- rbioUtil_perm_plot.default(baseplt = baseplt, plot.yLabel = "Total RMSEP", ...)
  }


  # save
  if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(perm_res)),".plsda.perm.plot.pdf...", sep = ""))  # initial message
  # grid.newpage()
  ggsave(filename = paste(deparse(substitute(perm_res)),".plsr.perm.plot.pdf", sep = ""), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  if(verbose) cat("Done! \n")
}


#' Title rbioUtil_perm_plot.default
#'
#' @rdname rbioUtil_perm_plot
#' @method rbioUtil_perm_plot default
#' @param baseplot A base plot generated by the \code{\link{rbioUtil_perm_plot}} method.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.xLabel X-axis label. Type with quotation marks. Could be NULL. Default is \code{"Permutation"}.
#' @param plot.xLabelSize X-axis label size. Default is \code{10}.
#' @param plot.xTickLblSize X-axis tick label size. Default is \code{10}.
#' @param plot.yLabel Y-axis label. Type with quotation marks. Could be NULL.
#' @param plot.yLabelSize Y-axis label size. Default is \code{10}.
#' @param plot.yTickLblSize Y-axis tick label size. Default is \code{10}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}.
#' @return A baseplot object for permutation results.
#' @import ggplot2
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @export
rbioUtil_perm_plot.default <- function(baseplt,
                                       plot.display.Title = TRUE, plot.titleSize = 10,
                                       plot.fontType = "sans",
                                       plot.xLabel = "Permutations", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                       plot.yLabel = "Accuracy", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                       plot.legendSize = 9,
                                       plot.rightsideY = TRUE){
  ## plot
  plt <- baseplt +
    ggtitle(ifelse(plot.display.Title, "Permutation test", NULL)) +
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

  if (plot.rightsideY){
    plt <- RBioplot::rightside_y(plt)
  }
  plt
}


#' @title rbioUtil_classplot
#'
#' @description Classification plot function for classification predction. The function uses \code{prediction} object generated from \code{\link{rbioClass_plsda_predict}} or \code{\link{rbioClass_svm_predict}} to generate classification probablity pie charts.
#' @param pred.object A \code{prediction} object, which can be obtained from funciton \code{\link{rbioClass_plsda_predict}}.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return  A \code{classification} obejct with classification probability summary for each sample, as well as pdf figure file fif \code{classplot = TRUE}.
#' @details The function operates in conjunction with the prediction function \code{\link{rbioClass_plsda_predict}}, to which the sample(s) of intestested is provided.
#' @import ggplot2
#' @import ggrepel
#' @importFrom grid grid.newpage grid.draw
#' @examples
#' \dontrun{
#' rbioUtil_classplot(pred.obj = new_model_optm_plsda_predict, multi_plot.ncol = 4, multi_plot.nrow = 4, plot.probLabelSize = 2)
#' }
#' @export
rbioUtil_classplot <- function(pred.obj,
                              multi_plot.ncol = nrow(pred.obj), multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
                              multi_plot.stripLblSize = 10,
                              plot.Title = NULL, plot.titleSize = 10,
                              plot.probLabelSize = 5, plot.probLabel.padding = 0,
                              plot.probLabel.outside = FALSE, plot.probLabel.outside.nudge = 1.5,
                              plot.fontType = "sans",
                              plot.legendSize = 9,
                              plot.Width = 170, plot.Height = 150,
                              verbose = TRUE){
  ## check arguments
  if (!any(class(pred.obj) %in% "prediction")) stop("pred.obj needs to be a  \"prediction\" class. Use functions like rbioClass_plsda_predict() to generate one.")
  if (pred.obj$model.type != "classification") stop("the function only supports \"classification\".")

  ## plot
  if (multi_plot.ncol * multi_plot.nrow < nrow(pred.obj$raw.newdata)){
    stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. Make sure they match the number of y groups.\n")
  }
  if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(pred.obj)),".plsda.classification.pdf...", sep = ""))  # initial message
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
  # grid.newpage()
  ggsave(filename = paste(deparse(substitute(pred.obj)),".plsda.classification.pdf", sep = ""), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  if (verbose) cat("Done!\n")

  invisible(plt)
}


#' @title uni_PreProc
#'
#' @description a legacy function borrowed from RBioArray. Data pre-processing function for the microarary data, which is now used for univariate reduction during svm nested CV.
#' @param rawlist Input data, either a list, \code{EList} or \code{MAList} object.
#' @param logTrans If to perfom a log transformation on the data or not. Default is \code{FALSE}.
#' @param logTransMethod If \code{logTrans = TRUE}, set which method to use for the transformation, \code{"log2"} or \code{"log10"}. Default is \code{"log2"}.
#' @param logTransObjT If \code{logTrans = TRUE}, set the file name for the output \code{csv} file containing the log transformed data.
#' @param logTransParallelComputing If \code{logTrans = TRUE}, set if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param bgMethod Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param normMethod Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @details The function does not use design matrix for array weight calculation.
#'          Therefore, the DE analysis based on the output from this function will yeild slightly different resutls from the \code{\link{rbioarray_filter_combine}}.
#'
#' @return Depending on the input type, the function outputs a \code{list}, \code{Elist} or \code{MAList} object with corrected and normalized expression values.
#'         If \code{logTrans = TRUE}, the function also outputs a \code{csv} file containing the log transformed data.
#' @import doParallel
#' @import foreach
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma backgroundCorrect normalizeBetweenArrays backgroundCorrect.matrix arrayWeights
#' @export
uni_PreProc <- function(rawlist, logTrans = FALSE, logTransMethod = "log2",
                              logTransObjT = "data", logTransParallelComputing = FALSE,
                              bgMethod = "auto", normMethod = "quantile", ...){
  if (class(rawlist) == "list"){
    ## log transform  or not
    if (logTrans){
      if (!logTransParallelComputing){
        # log transform
        mtx <- apply(rawlist$E, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))
      } else {
        # parallel computing
        # set up cpu cluster
        n_cores <- detectCores() - 1
        cl <- makeCluster(n_cores, type = "PSOCK")
        on.exit(stopCluster(cl)) # close connect when exiting the function
        # log transform
        mtx <- foreach(i = rawlist$E) %dopar% {
          out <- apply(i, c(1,2), FUN = ifelse(logTransMethod == "log10", log10, log2))
        }
      }
      tmpdata <- list(E = mtx, genes = rawlist$genes, targets = rawlist$targets)
      # store and export log transformed data into a csv file
      logTransOut <- data.frame(rawlist$genes, mtx)
      write.csv(logTransOut, file = paste(logTransObjT, "_log_transformed.csv", sep = ""), row.names = FALSE)
    } else {
      tmpdata <- rawlist
    }

    ## normalization
    BgC <- backgroundCorrect.matrix(tmpdata$E, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm) # array weight
    output <- list(E = Norm, genes = rawlist$genes, targets = rawlist$targets, ArrayWeight = Wgt)
  } else {
    ## normalization
    BgC <- backgroundCorrect(rawlist, method = bgMethod, ...) #background correction
    Norm <- normalizeBetweenArrays(BgC, method = normMethod) # quantile normalization
    Wgt <- arrayWeights(Norm)
    Norm$ArrayWeight <- Wgt
    output <- Norm
  }
  return(output)
}

