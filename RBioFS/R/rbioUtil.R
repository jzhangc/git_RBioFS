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

