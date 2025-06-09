#' Title rbioUtil_fscount_plot
#'
#' @description Bar graph for the FS count numbers
#' @param object Object or data frame containing FS count information. This currently supports S3 classes \code{table}, \code{data.frame} and \code{rbiosvm_nestedcv}.
#' @param ... Arguments for S3 class methods.
#' @return ggplot 2 bar graph exported as PDF file
#' @export
rbioUtil_fscount_plot <- function(object, ...){
  # -- class check --
  if (!class(object) %in% c("rbiosvm_nestedcv", "data.frame")) stop("object needs to be either a \"rbiosvm_nestedcv\" or \"data.frame\" object")

  # -- use method --
  UseMethod("rbioUtil_fscount_plot", object)
}


#' Title rbioUtil_fscount_plot.rbiosvm_nestedcv
#' @rdname rbioUtil_fscount_plot
#' @method rbioUtil_fscount_plot rbiosvm_nestedcv
#' @param object A \code{rbiosvm_nestedcv} object.
#' @param export.name Name used to output the plot.
#' @param ... Additional arguments for the default method.
#'
#' @export
rbioUtil_fscount_plot.rbiosvm_nestedcv <- function(object, export.name, ...) {
  x <- data.frame(object$best.nested.fs.count)
  threshold <- object$fs.count.threshold
  fvar <- "best.nested.fs"
  cvar <- "Freq"
  rbioUtil_fscount_plot.default(dfm = x, threshold = threshold,
                                export.name = export.name,
                                feature_var = fvar, count_var = cvar, ...)
}


#' Title rbioUtil_fscount_plot.default
#' @rdname rbioUtil_fscount_plot
#' @method rbioUtil_fscount_plot default
#' @param dfm Input data.frame.
#' @param feature_var Variable name for feature names.
#' @param count_var Variable name for FS count.
#' @param export.name Name used to output the plot.
#' @param threshold An integer indicating the FS count threshold for consensus feature selection.
#' @param plot.preview If to preview plots. Default is \code{TRUE}.
#' @param plot.title Whether to display plot title on top of the plot. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.outlineCol The outline colour for the bar graph. Default is \code{"gray4"}.
#' @param plot.fillCol The fill colour for the bar. Default is \code{"lightcyan"}
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is available on the website: \url{http://kenstoreylab.com/?page_id=2448}.
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
#' @details
#'    The input data.frame should at least have variables for feature name and FS count.
#'
#' @import ggplot2
#' @importFrom scales rescale_none
#' @export
rbioUtil_fscount_plot.default <- function(dfm, feature_var, count_var,
                                          export.name = "plot",
                                          threshold = 2,
                                          plot.preview = TRUE,
                                          plot.title = TRUE, plot.titleSize = 10,
                                          plot.sig.line = TRUE,
                                          plot.outlineCol = "gray4", plot.fillCol = "lightcyan",
                                          plot.errorbarLblSize = 6, plot.fontType = "sans",
                                          plot.xLabel = "Features", plot.xLabelSize = 10, plot.xTickLblSize = 10, plot.xTickItalic = FALSE,
                                          plot.xTickBold = FALSE, plot.xAngle = 90, plot.xhAlign = 1, plot.xvAlign = 0.5,
                                          plot.rightsideY = TRUE,
                                          plot.yLabel = "FS counts", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                          plot.yTickItalic = FALSE, plot.yTickBold = FALSE, plot.legendSize = 9,
                                          plot.legendTtl = FALSE, plot.legendTtlSize = 9,
                                          plot.Width = 170, plot.Height = 150) {
  # -- data prep --
  d <- dfm[dfm$Freq >= threshold, ]
  if (!feature_var %in% names(dfm)) stop("feature variable not found.")
  if (!count_var %in% names(dfm) ) stop("count variable not found.")

  # plot
  baseplt <- ggplot(data = d, aes(x = .data[[feature_var]], y = .data[[count_var]])) +
    geom_bar(position = "dodge", stat = "identity", color = plot.outlineCol, fill = plot.fillCol) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    scale_y_continuous(limits = c(0, d$Freq), expand = expansion(add = c(0, 1)), oob = scales::rescale_none, sec.axis = dup_axis(),
                       breaks = scales::breaks_pretty()) +
    xlab(plot.xLabel) +
    ylab(plot.yLabel) +
    geom_hline(yintercept = threshold, col = "red", linetype = "dashed") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType),
          axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
          axis.title.y = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
          axis.title.y.right = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = plot.legendSize),
          axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType, angle = plot.xAngle,
                                     hjust = plot.xhAlign, vjust = plot.xvAlign),
          axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5),
          axis.ticks.x = if(plot.xTickLblSize == 0) element_blank())

  p <- baseplt +
    annotate(geom = "text", label = paste0("FS count threshold = ", threshold), colour = "red", y = max(d$Freq), x = nrow(d), hjust = 0.9) +
    annotate(geom = "segment", x = nrow(d), y = threshold+0.5, xend = nrow(d), yend = max(d$Freq)-0.5,
             arrow = arrow(type = "closed", length = unit(0.02, "npc")))

  # -- export --
  if (plot.preview) show(p)
  ggsave(filename = paste0(export.name, ".fscount.pdf"), plot = p,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)

  # return(p)
}



#' Title rbioUtil_classif_accuracy
#'
#' @description Classification accuracy calculation function
#' @param object Classification model object. The object should be either \code{rbiosvm}, or \code{rbiomvr} classes.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.label Permutation test results object. The object should be either \code{rbiosvm_perm}, or \code{rbiomvr_perm} classes.
#' @param center.scale.newdata Logical, whether center and scale the newdata with training data mean and standard deviation. Default is \code{TRUE}.
#' @details accuracy = true predictions/total predictions
#' @return A list includes:
#'         (1) Numerical model accuracy value based on the input test data.
#'         (2) Confusion matrix
#'         (3) Input data label
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
  confusion_mtx <- table(predict(object = object, newdata = test), newdata.label, dnn = c("Predicted", "Actual"))
  accu <- sum(diag(confusion_mtx))/sum(confusion_mtx)  # true prediction/total prediction

  # output
  out <- list(accuracy = accu, confusion = confusion_mtx, newdata.label = newdata.label)

  return(out)
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
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
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
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
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
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is available on the website: \url{http://kenstoreylab.com/?page_id=2448}.
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
  baseplt <- baseplt +
    ggtitle(ifelse(plot.display.Title, "Permutation test", NULL)) +
    xlab(plot.xLabel) +
    ylab(plot.yLabel)

  if (plot.rightsideY) {
    plt <- baseplt +
      scale_y_continuous(sec.axis = dup_axis()) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
            axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
            axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
            axis.title.y.right = element_blank(),
            legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
            legend.key = element_blank(),
            axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
            axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
  } else {
    plt <- baseplt +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
            plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
            axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
            axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
            legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
            legend.key = element_blank(),
            axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
            axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
  }

  # below: not needed for ggplot 3.5.0
  # if (plot.rightsideY){
  #   plt <- RBioplot::rightside_y(plt)
  # }

  plt
}

#' @title rbioUtil_classplot
#'
#' @description Classification plot function for classification prediction. The function uses \code{prediction} object generated from \code{\link{rbioClass_plsda_predict}} or \code{\link{rbioClass_svm_predict}} to generate classification probablity pie charts.
#' @param pred.object A \code{prediction} object, which can be obtained from function \code{\link{rbioClass_plsda_predict}}.
#' @param export.name String. Optional user defined export name prefix. Default is \code{NULL}.
#' @param multi_plot.ncol Number of columns on one figure page. Default is \code{nrow(pred.obj)}.
#' @param multi_plot.nrow Number of rows on one figure page. Default is \code{1}.
#' @param multi_plot.legend.pos The legend position. Only effective when multi-plot is generated. Options are \code{"bottom"}, \code{"top"}, \code{"left"} and \code{"right"}. Default is \code{"bottom"}.
#' @param plot.Title Plot title. Default is \code{FALSE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
#' @param plot.probLabelSize The size of the sample label. Default is \code{2}.
#' @param plot.probLabel.padding The padding between sample symbol and the label. Default is \code{0.5}.
#' @param plot.probLabel.outside If to put the probability label outside of the pies. Default is \code{"TRUE"}.
#' @param plot.probLabel.outside.nudge Set only when \code{plot.probLabel.outside = TRUE}, adjustment to nudge the starting position of each label. Default is \code{"1.5"}.
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is available on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
#' @param plot.legendSize Legend size. Default is \code{9}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return  A \code{classification} object with classification probability summary for each sample, as well as pdf figure file when \code{classplot = TRUE}.
#' @details The function operates in conjunction with the prediction function \code{\link{rbioClass_plsda_predict}}, to which the sample(s) of interested is provided.
#' @import ggplot2
#' @import ggrepel
#' @importFrom grid grid.newpage grid.draw
#' @examples
#' \dontrun{
#' rbioUtil_classplot(pred.obj = new_model_optm_plsda_predict, multi_plot.ncol = 4, multi_plot.nrow = 4, plot.probLabelSize = 2)
#' }
#' @export
rbioUtil_classplot <- function(pred.obj, export.name = NULL,
                              multi_plot.ncol = nrow(pred.obj$raw.newdata),
                              multi_plot.nrow = 1, multi_plot.legend.pos = "bottom",
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
  if (is.null(export.name)){
    export.name <- deparse(substitute(pred.obj))
  } else {
    export.name <- export.name
  }

  ## plot
  if (multi_plot.ncol * multi_plot.nrow < nrow(pred.obj$raw.newdata)){
    stop("multi_plot.ncol and multi_plot.nrow settings are incorrect. Make sure they match the number of y groups.\n")
  }
  if (verbose) cat(paste0("Plot being saved to file: ", export.name, ".plsda.classification.pdf..."))  # initial message
  plt <- ggplot(pred.obj$probability.summary, aes(x = "", y = Probability, fill = Class)) +
    geom_col(width = 1, colour = "black", alpha = 0.8) +
    ggtitle(plot.Title) +
    xlab(NULL) +
    ylab(NULL) +
    facet_wrap(~Sample, nrow = multi_plot.nrow, ncol = multi_plot.ncol) +
    theme(strip.background = element_rect(fill = NA, colour = "black"),  # no strip background colour
          strip.text = element_text(face = "bold", size = multi_plot.stripLblSize),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
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
  ggsave(filename = paste0(export.name, ".", pred.obj$classifier.class[2], ".classification.pdf"), plot = plt,
         width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
  grid.draw(plt)
  if (verbose) cat("Done!\n")

  invisible(plt)
}


#' @title convertToFactor
#'
#' @description Convert data.frame features into factors.
#' @param dfm data.frame. Input data.frame.
#' @param cols str or vector of strings. Column names to convert.
#' @return data.frame with columns converted into factors.
#' @details
#'  The output data.frame has the same dim() as the input data.frame.
#' @examples
#' \dontrun{
#'   convertToFactor(dfm = dfm, cols = c("gender", "regions"))
#' }
#' @export
convertToFactor <- function(dfm, cols) {
  if (!"data.frame" %in% class(dfm)) stop("input dt needs to be a data.frame.")
  if (all(!cols %in% names(dfm))) stop("one or more cols not found in the input dataframe.")

  out <- dfm
  for (col in cols) {
    out[, col] <- factor(dfm[, col], levels = unique(dfm[, col]))
  }
  return(out)
}


#' @title dummy
#'
#' @description dummification of categorical variables
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
  if (!is.factor(x)){  # use this one for the arguments
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


#' @title rbioUtil_onehotEncoding
#'
#' @description Convert data.frame features into factors.
#' @param data data.frame. Input data.frame.
#' @param cols str or vector of strings. Column names to encode.
#' @param export_mode str. Export as data.table or data.frame. Default is \code{"dfm"}.
#' @return data.frame with columns one hot encoded.
#' @details
#'  The output data.frame has all the data of the inoput data, as opposed to just encoded columns.
#' @importFrom mltools one_hot
#' @importFrom data.table data.table
#' @examples
#' \dontrun{
#'   rbioUtil_onehotEncoding(dfm = data, cols = c("gender"), export_mode = "dfm")
#' }
#' @export
rbioUtil_onehotEncoding <- function(data, cols, export_mode = c("dfm", "dt")) {
  if (length(export_mode) > 1) stop("only set one value for \"export_mode\"")
  export_mode <- match.arg(export_mode)

  dt <- convertToFactor(data, cols)
  dt <- data.table(dt)
  out <- one_hot(dt, cols = cols)

  switch(export_mode,
         dfm = {
           return(as.data.frame(out))
         },
         dt = {
           return(out)
         })
}


#' @title substr_right
#'
#' @description Extract n'th character from string from right.
#' @param x str. Input string.
#' @param n int. Optional user defined export name prefix. Default is \code{NULL}.
#' @return character extracted from string.
#' @examples
#' \dontrun{
#'   substr_right("abc", 2)
#' }
#' @export
substr_right <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


#' @title na_summary
#'
#' @description NA ration summary.
#' @param data data.frame, matrix list of data.frames or matrices. Input data.
#' @param by str. Summarize by row or column. Default is \code{"row"}.
#' @return Vector or list of vector of NA ratios.
#' @examples
#' \dontrun{
#'   na_summary(data = data, by = "row")
#' }
#' @export
na_summary <- function(data, by = c("row", "col")) {
  if (length(by) >1) stop("only set one value for \"by\". ")
  by <- match.arg(by)

  if ("list" %in% class(data)) {
    out <- vector(mode = "list", length = length(data))
    names(out) <- names(data)

    switch(by,
           row = {
             for (i in 1:length(out)) {
               out[[i]] <- rowSums(is.na(data[[i]]))/ncol(data[[i]])
             }
           },
           col = {
             for (i in 1:length(out)) {
               out[[i]] <- colSums(is.na(data[[i]]))/nrow(data[[i]])
             }
           }
    )

  } else if (any(c("data.frame", "matrix") %in% class(data))) {
    out <- vector()
    switch(by,
           row = {
             for (i in 1:length(out)) {
               out <- rowSums(is.na(data))/ncol(data)
             }
           },
           col = {
             out <- colSums(is.na(data))/nrow(data)
           }
    )
  }
  return(out)
}


#' @title rbio_shap_svm_label_prob
#' @description
#' Helper function to extract probability vectors from svm model predictions.
#' @param object A SVM model object, in \code{rbiosvm} class.
#' @param newdata Input data to be classified. Make sure it is a \code{matrix} class and has the same variables as the model, i.e. same number of columns as the training data.
#' @param col_idx Probability column (i.e. label or outcome) index or name for which the probability vector is extracted.
#' @param ... Additional parameters for the \code{e1071:::predict.svm} method.
#' @details
#'    1. \code{newdata} should not include outcome column.
#'    2. \code{newdata} should be 0-1 scaled.
#'    3. \code{newdata} will be standardized using training data's col mean and sd, if the input \code{rbiosvm} model \code{object} contains standardization information:
#'       \code{object$center.scaledX}.This means the model was trainined using standardization.
#' @return Vector or data frame of prediction probabilities.
#' @import e1071
rbio_shap_svm_label_prob <- function(object, newdata, col_idx = NULL, ...) {
  # --- arg check ---
  if (!any(class(object) %in% c("rbiosvm"))) stop("object has to be a \"rbiosvm\" class.\n")

  # --- prediction ---
  if (!is.null(object$center.scaledX)) {
    x <-  t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
  }
  pred <- e1071:::predict.svm(object, x, probability = TRUE, ...)
  prob_dfm <- attr(pred, "probabilities")

  # --- prob extraction ---
  if (is.null(col_idx)) {
    o <- prob_dfm
  } else {
    o <- prob_dfm[, col_idx]
  }
  return(o)
}


#' @title shap_g_theme_helper
#' @description
#'    A helper function to set the plot themes for the \code{\link{rbioClass_svm_shap_aggreated}} function.
#' @import ggplot2
#' @return \code{ggplot2} graph object.
shap_g_theme_helper <- function(g,
                                plot.SymbolSize = 2, plot.lineSize = 1,
                                plot.display.Title = TRUE, plot.title = "SHAP", plot.titleSize = 10,
                                plot.fontType = "sans",
                                plot.xLabel = "SHAP value", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                plot.yLabel = "Features", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                plot.legendPosition = "right", plot.legendSize = 9,
                                plot.Width = 170, plot.Height = 150) {
  g <- g + ggtitle(ifelse(plot.display.Title, plot.title, NULL)) +
    xlab(plot.xLabel) +
    ylab(plot.yLabel) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
          axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
          legend.position = plot.legendPosition, legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
          legend.key = element_blank(),
          axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
          axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 1))

  g
}


#' @title sv_importance_shapviz
#' @description
#' A alternative version of \code{sv_importance} function from \code{shapviz} package
#'
#' @details
#'    Mainly to add custom graphic settings to the original function
#' @import shapviz
#' @return \code{ggplot2} graph object.
sv_importance_shapviz <- function(object, kind = c("bar", "beeswarm", "both", "no"),
                                  max_display = 15L, fill = "#fca50a", bar_width = 2/3,
                                  bee_width = 0.4, bee_adjust = 0.5,
                                  viridis_args = getOption("shapviz.viridis_args"),
                                  color_bar_title = "Feature value",
                                  show_numbers = FALSE, format_fun = format_max,
                                  number_size = 3.2, sort_features = TRUE, ...) {
  stopifnot("format_fun must be a function" = is.function(format_fun))
  kind <- match.arg(kind)
  imp <- shapviz:::.get_imp(get_shap_values(object), sort_features = sort_features)

  if (kind == "no") {
    return(imp)
  }

  # Deal with too many features
  if (ncol(object) > max_display) {
    imp <- imp[seq_len(max_display)]
  }
  ord <- names(imp)
  object <- object[, ord]  # not required for kind = "bar"

  # ggplot will need to work with data.frame
  imp_df <- data.frame(feature = factor(ord, rev(ord)), value = imp)
  is_bar <- kind == "bar"
  if (is_bar) {
    p <- ggplot2::ggplot(imp_df, ggplot2::aes(x = value, y = feature)) +
      ggplot2::geom_bar(fill = fill, width = bar_width, stat = "identity", colour = "black", alpha = 0.1, ...) +
      ggplot2::labs(x = "mean(|SHAP value|)", y = ggplot2::element_blank())
  } else {
    # Prepare data.frame for beeswarm plot
    S <- get_shap_values(object)
    X <- shapviz:::.scale_X(get_feature_values(object))
    df <- transform(
      as.data.frame.table(S, responseName = "value"),
      feature = factor(Var2, levels = rev(ord)),
      color = as.data.frame.table(X)$Freq
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature))
    if (kind == "both") {
      p <- p +
        ggplot2::geom_bar(
          data = imp_df, fill = fill, width = bar_width, stat = "identity",  colour = "black", alpha = 0.1,
        )
    }
    p <- p +
      ggplot2::geom_vline(xintercept = 0, color = "darkgray") +
      ggplot2::geom_point(
        ggplot2::aes(color = color),
        position = shapviz:::position_bee(width = bee_width, adjust = bee_adjust),
        ...
      ) +
      shapviz:::.get_color_scale(
        viridis_args = viridis_args,
        bar = !is.null(color_bar_title),
        ncol = length(unique(df$color))   # Special case of constant feature values
      ) +
      ggplot2::labs(
        x = "SHAP value", y = ggplot2::element_blank(), color = color_bar_title
      ) +
      ggplot2::theme(legend.box.spacing = grid::unit(0, "pt"))
  }
  if (show_numbers) {
    p <- p +
      ggplot2::geom_text(
        data = imp_df,
        ggplot2::aes(
          x = if (is_bar) value + max(value) / 60 else
            min(df$value) - diff(range(df$value)) / 20,
          label = format_fun(value)
        ),
        hjust = !is_bar,
        size = number_size
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = 0.05 + c(0.12 *!is_bar, 0.09 * is_bar))
      )
  }
  p
}
