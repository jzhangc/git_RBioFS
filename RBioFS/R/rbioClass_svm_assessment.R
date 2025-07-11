#' @title rbioReg_svm_r2
#'
#' @description Support Vector Regression (SVR) R2 calculation
#' @param object A \code{rbiosvm} object.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.y For regression model only, the vector for the new data's continuous outcome variable. Default is \code{NULL}
#' @return RMSE value with either new data or training data.
#'
#' @details
#'
#' The R2 is calculated as following: 1 - rss/tss.
#' rss: residual sum of squares: sum((yhat-y)^2)
#' tss: total sum of squares: sum((y-mean(y))^2)
#'
#' @examples
#' \dontrun{
#'  test_r2 <- rbioReg_svm_r2(object = svm_m, newdata = svm_test[,-1], newdata.y = svm_test$y)
#' }
#' @export
rbioReg_svm_r2 <- function(object, newdata=NULL, newdata.y=NULL){
  # argument check
  if (!any(class(object) %in% "rbiosvm")) stop("The input object needs to be a \"rbiosvm\" class.")
  if (object$model.type != "regression") stop("The input model type needs to be \"regression\".")
  if (is.null(newdata)) {
    stop("newdata cannot be NULL.")
  } else {
    if (is.null(newdata.y)) stop("newdata.y needs to be specified if newdata is available.")
    if (nrow(newdata) != length(newdata.y)) stop("sample size needs to match the length of newdata.y.")
  }

  # computation
  center_scale_newdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
  pred <- predict(object, center_scale_newdata)
  err <- newdata.y - pred
  ss_res <- sum(err^2)  # sum of residual squares
  ss_tot <- sum((newdata.y - mean(newdata.y))^2)  # sum of squares
  out_r2 <- 1 - ss_res/ss_tot

  # output
  return(out_r2)
}


#' @title rbioClass_svm_roc_auc
#'
#' @description ROC-AUC analysis and ploting for SVM model, for classification model only.
#' @param object A \code{rbiosvm} object.
#' @param fileprefix String. A file prefix to use for export file name, instead of the objecte name. Default is \code{NULL}.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.label The correspoding label vector to the data. Make sure it is a \code{factor} object. Defaults is \code{NULL}.
#' @param center.scale.newdata Logical, whether center and scale the newdata with training data mean and standard deviation. Default is \code{TRUE}.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.smooth If to smooth the curves. Uses binormal method to smooth the curves. Default is \code{FALSE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Prints AUC values in the console. And a pdf file for ROC plot. The function also exports a ROC results list as a \code{svm_roc_auc} class to the environment.
#'
#'         Items of the \code{svm_roc_auc} class:
#'         \code{model.type}
#'         \code{y.threshold}
#'         \code{regression.categories}
#'         \code{svm.roc_object}
#'         \code{svm.roc_dataframe}
#'         \code{input.newdata}
#'         \code{input.newdata.y}
#'         \code{input.newdata.label}
#'         \code{newdata.center.scaled}
#'
#' @details Uses pROC module to calculate ROC. The function supports more than two groups or more than one threshold for classification and regression model.
#'
#'          When \code{newdata} is not provided, the function uses the training data from the input SVM object, and the training data is automatically cneter.scaled.
#'
#'          Although optional, the \code{newdata} matrix should use training data's column mean and column standard deviation to center.scale prior to ROC-AUC analysis.
#'          The option \code{center.scaled.newdata = FALSE} is used when the whole (training and test sets) data were center.scaled before SVM training and testing.
#'
#'          For regression models, it is ok to provide multiple thresholds for \code{y.threshold}.
#'
#'          Along with the data frame, the \code{svm_roc_auc} output includes the \code{roc} object from the \code{pROC} package, with which the stats can be done to compare ROCs.
#'
#'          Also: a ROC curve with AUC == 1 is always 1-1 for 95% CI and can be misleading.
#'
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioClass_svm_roc_auc(object = svm_m, newdata = svm_test[, -1],
#'                       newdata.label = factor(svm_test$y, levels = unique(svm_test$y)))
#' }
#' @export
rbioClass_svm_roc_auc <- function(object, fileprefix = NULL,
                                  newdata = NULL, newdata.label = NULL,
                                  center.scale.newdata = TRUE,
                                  rocplot = TRUE,
                                  plot.smooth = FALSE,
                                  plot.SymbolSize = 2, plot.lineSize = 1,
                                  plot.display.Title = TRUE, plot.titleSize = 10,
                                  plot.fontType = "sans",
                                  plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                  plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                  plot.legendSize = 9, plot.rightsideY = TRUE,
                                  plot.Width = 170, plot.Height = 150,
                                  verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c('rbiosvm'))) stop("object needs to be \"rbiosvm\" class.")
  if (object$model.type != "classification") stop("object needs to be a \"classification\" model")
  if (is.null(newdata)) {
    cat("No newdata input, proceed with training data.\n\n")
    newdata <- object$inputX
    newdata.label <- object$inputY
    if (any(class(newdata.label) != "factor")){
      if (verbose) cat("newdata.label is converted to factor. \n")
      newdata.label <- factor(newdata.label, levels = unique(newdata.label))
    }
    center.scale.newdata <- TRUE # automatically center.scale training data X when no new data is provided
  }
  if (!any(class(newdata) %in% c("data.frame", "matrix")) & !is.null(dim(newdata))) stop("newdata needs to be a matrix, data.frame or vector.")
  if (any(class(newdata) == "data.frame") | is.null(dim(newdata))){
    if (verbose) cat("newdata converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (ncol(newdata) != ncol(object$inputX)) stop("test data should have the same number of variables as the training data.")
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data while center.scale.newdata = TRUE.")
  }

  ## process data
  if (center.scale.newdata){ # using training data mean and sd
    if (verbose) cat(paste0("Data center.scaled using training data column mean and sd, prior to modelling.\n"))
    centered_newdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
    test <- centered_newdata
  } else {
    centered_newdata <- NULL
    test <- newdata
  }

  ## ROC-AUC calculation
  pred <- predict(object, newdata = test, decision.values = TRUE, probability = TRUE)  # prediction
  pred_prob <- attr(pred, "probabilities")
  # if (object$model.type == "regression"){
  #   pred <- cut(pred, rg, labels = reg.group.names)
  # }
  outcome <- newdata.label  # origial label

  # auc
  roc_auc_list <- vector(mode = "list", length = length(levels(outcome)))
  roc_auc_list[] <- foreach(i = 1:length(levels(outcome))) %do% {
    response <- outcome
    predictor <- pred_prob[, levels(response)[i]]  # probability of the current outcome
    levels(response)[-i] <- "others"
    # predictor <- dummy(pred)
    # predictor <- as.matrix(predictor[, i], ncol = 1)
    splt <- split(predictor, response)  # split function splist array according to a factor
    controls <- splt$others
    cases <- splt[[levels(outcome)[i]]]
    if (length(cases) && length(controls)) {  # check if there are any zero controls or zero cases
      perf <- tryCatch(suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth, ci= TRUE))),
                       error = function(err){
                         cat("Curve not smoothable. Proceed without smooth.\n")
                         suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = FALSE, ci= TRUE)))
                       })
      if (length(levels(outcome)) == 2){
        cat(paste0("AUC - ", levels(outcome)[i], ": ", perf$auc, "\n"))
      } else {
        cat(paste0(" AUC - ", levels(outcome)[i], " (vs Others): ", perf$auc, "\n"))
      }
      perf
    } else {
      return(NULL)
    }
  }
  names(roc_auc_list) <- unique(outcome)

  roc_dfm <- foreach(i = 1:length(levels(outcome)), .combine = "rbind") %do% {
    perf <- roc_auc_list[[i]]
    fpr <- 1 - perf$specificities
    tpr <- perf$sensitivities
    thresholds <- perf$thresholds
    mtx <- cbind(fpr, tpr, thresholds)
    if (length(levels(outcome)) == 2){
      df <- data.frame(mtx, group = rep(levels(outcome)[i], times = nrow(mtx)), row.names = NULL, check.names = FALSE)
    } else {
      df <- data.frame(mtx, group = rep(paste0(levels(outcome)[i], " (vs Others)"), times = nrow(mtx)), row.names = NULL, check.names = FALSE)
    }
    df <- df[order(df$tpr), ]
    return(df)
  }

  ## return
  if (any(sapply(roc_auc_list, is.null))) {
    cat("Either no control or no case obesevated. ROC-AUC fails. \n")
    out <- NULL
  } else {
    out <- list(model.type = object$model.type,
                svm.roc_object = roc_auc_list,
                svm.roc_dataframe = roc_dfm,
                input.newdata = newdata,
                input.newdata.label = newdata.label,
                newdata.center.scaled = centered_newdata)
    class(out) <- "svm_roc_auc"
    if (is.null(fileprefix)) {  # export
      assign(paste0(deparse(substitute(object)), "_svm_roc_auc"), out, envir = .GlobalEnv)
    } else {
      assign(paste0(as.character(fileprefix), "_svm_roc_auc"), out, envir = .GlobalEnv)
    }

    ## plotting
    if (rocplot){
      if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".svm.roc.pdf...", sep = ""))  # initial message

      plt <- ggplot(data = roc_dfm, aes(x = fpr, y = tpr, group = group, colour = group)) +
        geom_line(aes(linetype = group), linewidth = plot.lineSize) +
        geom_point(aes(shape = group), size = plot.SymbolSize) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggtitle(ifelse(plot.display.Title, "ROC", NULL)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel)

      if (plot.rightsideY) {
        plt <- plt +
          scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA), sec.axis = dup_axis()) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                axis.title.y.right = element_blank(),
                legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                legend.key = element_blank(),
                axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
      } else {
        plt <- plt +
          scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA)) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                legend.key = element_blank(),
                axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
      }

      # # below: not needed for ggplot2 3.5.0
      # if (plot.rightsideY){
      #   plt <- RBioplot::rightside_y(plt)
      # }

      # save
      # grid.newpage()
      if (is.null(fileprefix)) {
        ggsave(filename = paste(deparse(substitute(object)),".svm.roc.pdf", sep = ""), plot = plt,
               width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      } else {
        ggsave(filename = paste(as.character(fileprefix),".svm.roc.pdf", sep = ""), plot = plt,
               width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      }
      grid.draw(plt)
      if (verbose) cat("Done!\n")
    }
  }
}


#' @title rbioClass_svm_roc_auc_inter
#'
#' @description Interpolated ROC-AUC analysis and ploting for SVM model, for classification model only.
#' @param object A \code{rbiosvm} object.
#' @param fileprefix String. A file prefix to use for export file name, instead of the objecte name. Default is \code{NULL}.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.label The corresponding label vector to the data. Make sure it is a \code{factor} object. Defaults is \code{NULL}.
#' @param center.scale.newdata Logical, whether center and scale the newdata with training data mean and standard deviation. Default is \code{TRUE}.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.smooth If to smooth the curves. Uses binormal method to smooth the curves. Default is \code{FALSE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Prints AUC values in the console. And a pdf file for ROC plot. The function also exports a ROC results list as a \code{svm_roc_auc_inter} class to the environment.
#'         The \code{svm_roc_auc_inter} class is a subclass to \code{svm_roc_auc}, meaning it contains all the items from the latter:
#'
#'         \code{model.type}
#'         \code{y.threshold}
#'         \code{regression.categories}
#'         \code{svm.roc_object}
#'         \code{svm.roc_dataframe}
#'         \code{input.newdata}
#'         \code{input.newdata.y}
#'         \code{input.newdata.label}
#'         \code{newdata.center.scaled}
#'
#'         with the additional interpolation item:
#'         \code{svm.interporlated_roc_object}
#'
#' @details Uses pROC module to calculate ROC. The function supports more than two groups or more than one threshold for classification and regression model.
#'
#'          When \code{newdata} is not provided, the function uses the training data from the input SVM object, and the training data is automatically cneter.scaled.
#'
#'          Although optional, the \code{newdata} matrix should use training data's column mean and column standard deviation to center.scale prior to ROC-AUC analysis.
#'          The option \code{center.scaled.newdata = FALSE} is used when the whole (training and test sets) data were center.scaled before SVM training and testing.
#'
#'          For regression models, it is ok to provide multiple thresholds for \code{y.threshold}.
#'
#'          Along with the data frame, the \code{svm_roc_auc} output includes the \code{roc} object from the \code{pROC} package, with which the stats can be done to compare ROCs.
#'
#'          Also: a ROC curve with AUC == 1 is always 1-1 for 95% CI and can be misleading.
#'
#'          Interporlation is computed with the \code{approx} function, with FPR set at 0-1, with 100 intervals.
#'
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioClass_svm_roc_auc_itner(object = svm_m, newdata = svm_test[, -1],
#'                       newdata.label = factor(svm_test$y, levels = unique(svm_test$y)))
#' }
#' @export
rbioClass_svm_roc_auc_inter <- function(object, fileprefix = NULL,
                                        newdata = NULL, newdata.label = NULL,
                                        center.scale.newdata = TRUE,
                                        rocplot = TRUE,
                                        plot.smooth = FALSE,
                                        plot.SymbolSize = 2, plot.lineSize = 1,
                                        plot.display.Title = TRUE, plot.titleSize = 10,
                                        plot.fontType = "sans",
                                        plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                        plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                        plot.legendSize = 9, plot.rightsideY = TRUE,
                                        plot.Width = 170, plot.Height = 150,
                                        verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c('rbiosvm'))) stop("object needs to be \"rbiosvm\" class.")
  if (object$model.type != "classification") stop("object needs to be a \"classification\" model")
  if (is.null(newdata)) {
    cat("No newdata input, proceed with training data.\n\n")
    newdata <- object$inputX
    newdata.label <- object$inputY
    if (any(class(newdata.label) != "factor")){
      if (verbose) cat("newdata.label is converted to factor. \n")
      newdata.label <- factor(newdata.label, levels = unique(newdata.label))
    }
    center.scale.newdata <- TRUE # automatically center.scale training data X when no new data is provided
  }
  if (!any(class(newdata) %in% c("data.frame", "matrix")) & !is.null(dim(newdata))) stop("newdata needs to be a matrix, data.frame or vector.")
  if (any(class(newdata) == "data.frame") | is.null(dim(newdata))){
    if (verbose) cat("newdata converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (ncol(newdata) != ncol(object$inputX)) stop("test data should have the same number of variables as the training data.")
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data while center.scale.newdata = TRUE.")
  }

  ## process data
  if (center.scale.newdata){ # using training data mean and sd
    if (verbose) cat(paste0("Data center.scaled using training data column mean and sd, prior to modelling.\n"))
    centered_newdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
    test <- centered_newdata
  } else {
    centered_newdata <- NULL
    test <- newdata
  }

  ## ROC-AUC calculation
  pred <- predict(object, newdata = test, decision.values = TRUE, probability = TRUE)  # prediction
  pred_prob <- attr(pred, "probabilities")
  outcome <- newdata.label  # origial label

  # auc
  roc_auc_list <- vector(mode = "list", length = length(levels(outcome)))
  roc_auc_list[] <- foreach(i = 1:length(levels(outcome))) %do% {
    response <- outcome
    predictor <- pred_prob[, levels(response)[i]]  # probability of the current outcome
    levels(response)[-i] <- "others"
    # predictor <- dummy(pred)
    # predictor <- as.matrix(predictor[, i], ncol = 1)
    splt <- split(predictor, response)  # split function splist array according to a factor
    controls <- splt$others
    cases <- splt[[levels(outcome)[i]]]
    if (length(cases) && length(controls)) {  # check if there are any zero controls or zero cases
      perf <- tryCatch(suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth, ci= TRUE))),
                       error = function(err){
                         cat("Curve not smoothable. Proceed without smooth.\n")
                         suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = FALSE, ci= TRUE)))
                       })
      if (length(levels(outcome)) == 2){
        cat(paste0("AUC - ", levels(outcome)[i], ": ", perf$auc, "\n"))
      } else {
        cat(paste0(" AUC - ", levels(outcome)[i], " (vs Others): ", perf$auc, "\n"))
      }
      perf
    } else {
      return(NULL)
    }
  }
  names(roc_auc_list) <- unique(outcome)

  roc_dfm <- foreach(i = 1:length(levels(outcome)), .combine = "rbind") %do% {
    perf <- roc_auc_list[[i]]
    fpr <- 1 - perf$specificities
    tpr <- perf$sensitivities
    thresholds <- perf$thresholds
    mtx <- cbind(fpr, tpr, thresholds)
    if (length(levels(outcome)) == 2){
      df <- data.frame(mtx, group = rep(levels(outcome)[i], times = nrow(mtx)), row.names = NULL, check.names = FALSE)
    } else {
      df <- data.frame(mtx, group = rep(paste0(levels(outcome)[i], " (vs Others)"), times = nrow(mtx)), row.names = NULL, check.names = FALSE)
    }
    df <- df[order(df$tpr), ]
    return(df)
  }

  ## return
  if (any(sapply(roc_auc_list, is.null))) {
    cat("Either no control or no case obesevated. ROC-AUC fails. \n")
    out <- NULL
  } else {
    # construct a  list containing data.frames of AUC results over folds per group
    # for every fold: list = roc dfm w cols: tpr, fpr, folds, auc dfm w cols: auc, folds
    # use interpolation for CV ROC: set fpr seq(0, 1, length.out = 100), and interpolated tpr
    # to make sure all ROCs from different CV classifiers to
    # have the same length
    # interpolation for ROC thus is also used to plot ROCs for different types classifiers
    # for comparison
    auc_res_group <- names(roc_auc_list)
    fpr_init <- seq(0, 1, length.out = 100)
    auc_inter_list <- vector(mode = "list", length = length(auc_res_group))
    for (i in 1:length(auc_res_group)){
      group_label <- auc_res_group[[i]]
      if (verbose) cat(paste0("processing group label: ",  group_label, "\n"))

      # original roc-auc
      tpr <- roc_auc_list[[i]]$sensitivities
      fpr <- 1 - roc_auc_list[[i]]$specificities
      auc <- roc_auc_list[[i]]$auc
      thresholds <- roc_auc_list[[i]]$thresholds

      tpr_inter <- approx(fpr, tpr, xout = fpr_init, ties = "mean")$y
      fpr_inter <- fpr_init

      # interporlated roc-auc
      auc_inter_list[[i]] <- list(og_roc_auc = list(tpr = tpr, fpr = fpr, thresholds = thresholds),
                                  tpr_inter = tpr_inter,
                                  fpr_inter = fpr_inter)
      # break
    }
    names(auc_inter_list) <- auc_res_group

    # export
    out <- list(model.type = object$model.type,
                svm.roc_object = roc_auc_list,
                svm.interporlated_roc_object = auc_inter_list,
                svm.roc_dataframe = roc_dfm,
                input.newdata = newdata,
                input.newdata.label = newdata.label,
                newdata.center.scaled = centered_newdata)
    class(out) <- c("svm_roc_auc", "svm_roc_auc_inter")


    if (is.null(fileprefix)) {  # export
      assign(paste0(deparse(substitute(object)), "_svm_roc_auc_inter"), out, envir = .GlobalEnv)
    } else {
      assign(paste0(as.character(fileprefix), "_svm_roc_auc_inter"), out, envir = .GlobalEnv)
    }

    ## plotting
    if (rocplot){
      if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".svm.roc.pdf...", sep = ""))  # initial message

      auc_inter_list$high$tpr_inter
      roc_inter_dfm <- foreach(i = 1:length(levels(outcome)), .combine = "rbind") %do% {
        roc_inter <- auc_inter_list[[i]]
        fpr <- roc_inter$fpr_inter
        tpr <- roc_inter$tpr_inter
        mtx <- cbind(fpr, tpr)
        if (length(levels(outcome)) == 2){
          init_df <- data.frame(fpr = 0, tpr = 0,  group = levels(outcome)[i])
          end_df <- data.frame(fpr = 1, tpr = 1,  group = levels(outcome)[i])
          df <- data.frame(mtx, group = rep(levels(outcome)[i], times = nrow(mtx)), row.names = NULL, check.names = FALSE)
          df <- rbind(init_df, df)
          df <- rbind(df, end_df)
        } else {
          init_df <- data.frame(fpr = 0, tpr = 0,  group = paste0(levels(outcome)[i], " (vs Others)"))
          end_df <- data.frame(fpr = 1, tpr = 1,  group = paste0(levels(outcome)[i], " (vs Others)"))
          df <- data.frame(mtx, group = rep(paste0(levels(outcome)[i], " (vs Others)"), times = nrow(mtx)), row.names = NULL, check.names = FALSE)
          df <- rbind(init_df, df)
          df <- rbind(df, end_df)
        }
        df <- df[order(df$tpr), ]
        return(df)
      }

      plt <- ggplot(data = roc_inter_dfm, aes(x = fpr, y = tpr, group = group, colour = group)) +
        scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, NA)) +
        geom_line(aes(linetype = group), linewidth = plot.lineSize) +
        geom_point(aes(shape = group), size = plot.SymbolSize) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggtitle(ifelse(plot.display.Title, "ROC", NULL)) +
        xlab(plot.xLabel) +
        ylab(plot.yLabel)

      if (plot.rightsideY) {
        plt <- plt +
          scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA), sec.axis = dup_axis()) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                axis.title.y.right = element_blank(),
                legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                legend.key = element_blank(),
                axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
      } else {
        plt <- plt +
          scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA)) +
          theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                legend.key = element_blank(),
                axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
      }

      if (is.null(fileprefix)) {
        ggsave(filename = paste(deparse(substitute(object)),".svm.roc_inter.pdf", sep = ""), plot = plt,
               width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      } else {
        ggsave(filename = paste(as.character(fileprefix),".svm.roc_inter.pdf", sep = ""), plot = plt,
               width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
      }
      grid.draw(plt)
      if (verbose) cat("Done!\n")
    }
  }
}


#' @title svm_cv_rocauc_helper
#'
#' @description ROC-AUC helper function for cross-validated SVM models
#' @param object A \code{rbiosvm_cv} or \code{rbiosvm_nestedcv} object.
#' @param fileprefix String. A file prefix to use for export file name, instead of the object name. Default is \code{NULL}.
#' @param roc.smooth If to smooth the curves. Uses binomial method to smooth the curves. Default is \code{FALSE}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Returns a \code{cv_auc_res} (also a \code{list}) class object for follow up use.
#'
#' @details Uses pROC module to calculate ROC. The function supports more than two groups or more than one threshold
#'          for classification and regression model.
#'
#'          The test data included in the input \code{rbiosvm_nestedcv} is already normalized with training data information.
#'          So there is no need for additional data transformation.
#'
#'          Also: a ROC curve with AUC == 1 is always 1-1 for 95% CI and can be misleading.
#'
#'          The \code{cv_auc_res} class includes the following items:
#'          \code{cv_res_list}: a list containing CV ROC-AUC results organized per CV folds
#'          \code{smooth_roc}: if the ROC calculation is smoothed
#'          \code{input_model_class}: input CV model type, \code{rbiosvm_cv} (CV) or \code{rbiosvm_nestedcv} (nested CV)
#'
#' @import foreach
#' @importFrom pROC roc
#' @examples
#' \dontrun{
#' rbioClass_svm_cv_roc_auc(object = svm_nested_cv)
#' }
#' @export
svm_cv_rocauc_helper <- function(object, roc.smooth = FALSE, verbose = TRUE) {
  # ---- arguments check ----
  if (!any(class(object) %in% c('rbiosvm_nestedcv', 'rbiosvm_cv'))) stop("object needs to be \"rbiosvm_nestedcv\" or \"rbiosvm_cv\" classes.")
  input_class_idx <- c('rbiosvm_nestedcv', 'rbiosvm_cv') %in% class(object)
  input_class <- c('rbiosvm_nestedcv', 'rbiosvm_cv')[input_class_idx]

  # ---- read in data ----
  if (class(object) == 'rbiosvm_nestedcv') {
    cv_model_list <- object$nested.cv.models
  } else {
    cv_model_list <- object$cv.models
  }

  if (object$model.type == "regression"){
    stop('ROC-AUC only applies to classfication models.')
  }

  # --- check the validity of and, if needed, process the model list ---
  for (m in names(cv_model_list)) {
    if (verbose) cat(paste0("processing model: ", m))
    if ("simpleError" %in% class(cv_model_list[[m]])) {
      cv_model_list[[m]] <- NULL
    }
    if (verbose) cat("\n")
  }

  if (length(cv_model_list) < 1) stop("No valid model found in the object.")

  # ---- intermediate function ----
  # x = cv_model_list[[1]]
  cv_auc_func <- function(x){
    # data
    cv_test <- x$cv_test_data
    cv_test_x <- cv_test[, !names(cv_test) %in% "y"]
    # cv_test_y <- cv_test$y
    cv_test_y <- droplevels(cv_test$y) # to deal with fragments without all the classes

    # model
    cv_m <- x$cv_svm_model

    # pred
    pred <- predict(cv_m, newdata = cv_test_x, decision.values = TRUE, probability = TRUE)  # prediction
    pred_prob <- attr(pred, "probabilities")

    # ROC-AUC
    roc_auc_list <- vector(mode = "list", length = length(levels(cv_test_y)))
    roc_auc_list[] <- foreach(i = 1:length(levels(cv_test_y))) %do% {
      response <- cv_test_y
      predictor <- pred_prob[, levels(response)[i]]  # probability of the current outcome
      levels(response)[-i] <- "others"
      splt <- split(predictor, response)  # split function splist array according to a factor
      controls <- splt$others
      cases <- splt[[levels(cv_test_y)[i]]]
      if (length(cases) && length(controls)) {
        perf <- tryCatch(suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = roc.smooth, ci= TRUE))),
                         error = function(err){
                           cat("Curve not smoothable. Proceed without smooth.\n")
                           suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = FALSE, ci= TRUE)))
                           roc.smooth <- FALSE
                         })
        if (verbose){
          if (length(levels(cv_test_y)) == 2){
            cat(paste0("AUC - ", levels(cv_test_y)[i], ": ", perf$auc, "\n"))
          } else {
            cat(paste0(" AUC - ", levels(cv_test_y)[i], " (vs Others): ", perf$auc, "\n"))
          }
        }

        perf
      } else {
        return(NULL)
      }
    }
    names(roc_auc_list) <- unique(cv_test_y)

    roc_dfm <- foreach(i = 1:length(levels(cv_test_y)), .combine = "rbind") %do% {
      perf <- roc_auc_list[[i]]
      fpr <- 1 - perf$specificities
      tpr <- perf$sensitivities
      thresholds <- perf$thresholds
      mtx <- cbind(fpr, tpr, thresholds)
      if (length(levels(cv_test_y)) == 2){
        df <- data.frame(mtx, group = rep(levels(cv_test_y)[i], times = nrow(mtx)), row.names = NULL, check.names = FALSE)
      } else {
        df <- data.frame(mtx, group = rep(paste0(levels(cv_test_y)[i], " (vs Others)"), times = nrow(mtx)), row.names = NULL, check.names = FALSE)
      }
      df <- df[order(df$tpr), ]
      return(df)
    }

    # return
    if (any(sapply(roc_auc_list, is.null))) {
      out <- NULL
    } else {
      out <- list(svm.roc_object = roc_auc_list,
                  svm.roc_dataframe = roc_dfm)
    }
    return(out)
  }

  # ---- compute and return ----
  auc_res_list <- vector(mode = "list", length = length(cv_model_list))
  auc_res_list[] <- foreach(i = 1:length(cv_model_list)) %do% {
    out <- cv_auc_func(cv_model_list[[i]])
    out
  }
  names(auc_res_list) <- names(cv_model_list)


  if (any(sapply(auc_res_list, is.null))){
    cat("NULL element removed from the CV ROC-AUC list.\n")
    auc_res_list <- auc_res_list[-which(sapply(auc_res_list, is.null))]
  }

  # output results
  out_list <- vector(mode = "list", length = 3)
  names(out_list) <- c("cv_res_list", "smooth_roc", "input_model_class")

  out_list$cv_res_list <- auc_res_list
  out_list$smooth_roc <- roc.smooth
  out_list$input_model_class <- input_class

  # the cv_auc_res class is to make sure the followup functions
  # takes the correct input
  class(out_list) <- c("list", "cv_auc_res")
  return(out_list)
}


#' @title rbioClass_svm_cv_roc_auc
#'
#' @description ROC-AUC analysis and ploting for SVM cross-validation (CV) models, for classification only.
#' @param object A \code{rbiosvm_nestedcv} or \code{rbiosvm_cv} object.
#' @param fileprefix String. A file prefix to use for export file name, instead of the object name. Default is \code{NULL}.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.smooth If to smooth the curves. Uses binomial method to smooth the curves. Default is \code{FALSE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Prints AUC values in the console. And a pdf file for ROC plot.
#'         The function also exports a ROC results list with an ROC object for each CV fold.
#'
#' @details Uses pROC module to calculate ROC. The function supports more than two groups or more than one threshold
#'          for classification and regression model.
#'
#'          The test data included in the input \code{rbiosvm_nestedcv} is already normalized with training data information.
#'          So there is no need for additional data transformation.
#'
#'          Also: a ROC curve with AUC == 1 is always 1-1 for 95% CI and can be misleading.
#'
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioClass_svm_cv_roc_auc(object = svm_nested_cv)
#' }
#' @export
rbioClass_svm_cv_roc_auc <- function(object, fileprefix = NULL,
                                     rocplot = TRUE,
                                     plot.smooth = FALSE,
                                     plot.lineSize = 1,
                                     plot.display.Title = TRUE, plot.titleSize = 10,
                                     plot.fontType = "sans",
                                     plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                     plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                     plot.legendSize = 9, plot.rightsideY = TRUE,
                                     plot.Width = 170, plot.Height = 150,
                                     verbose = TRUE){
  # ---- argements check ----
  if (!any(class(object) %in% c('rbiosvm_nestedcv', 'rbiosvm_cv'))) stop("object needs to be \"rbiosvm_nestedcv\" or \"rbiosvm_cv\" classes.")

  # ---- read in data ----
  if (class(object) == 'rbiosvm_nestedcv') {
    cv_model_list <- object$nested.cv.models
  } else {
    cv_model_list <- object$cv.models
  }

  if (object$model.type == "regression"){
    stop('ROC-AUC only applies to classfication models.')
  }

  # --- check the validity of and, if needed, process the model list ---
  for (m in names(cv_model_list)) {
    if (verbose) cat(paste0("processing model: ", m))
    if ("simpleError" %in% class(cv_model_list[[m]])) {
      cv_model_list[[m]] <- NULL
    }
    if (verbose) cat("\n")
  }

  if (length(cv_model_list) < 1) stop("No valid model found in the object.")

  # ---- intermediate function ----
  cv_auc_func <- function(x){
    # data
    cv_test <- x$cv_test_data
    cv_test_x <- cv_test[, !names(cv_test) %in% "y"]
    # cv_test_y <- cv_test$y
    cv_test_y <- droplevels(cv_test$y) # to deal with fragments without all the classes

    # model
    cv_m <- x$cv_svm_model

    # pred
    pred <- predict(cv_m, newdata = cv_test_x, decision.values = TRUE, probability = TRUE)  # prediction
    pred_prob <- attr(pred, "probabilities")

    # ROC-AUC
    roc_auc_list <- vector(mode = "list", length = length(levels(cv_test_y)))
    roc_auc_list[] <- foreach(i = 1:length(levels(cv_test_y))) %do% {
      response <- cv_test_y
      predictor <- pred_prob[, levels(response)[i]]  # probability of the current outcome
      levels(response)[-i] <- "others"
      splt <- split(predictor, response)  # split function splist array according to a factor
      controls <- splt$others
      cases <- splt[[levels(cv_test_y)[i]]]
      if (length(cases) && length(controls)) {
        perf <- tryCatch(suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth, ci= TRUE))),
                         error = function(err){
                           cat("Curve not smoothable. Proceed without smooth.\n")
                           suppressWarnings(suppressMessages(pROC::roc(controls = controls, cases = cases, smooth = FALSE, ci= TRUE)))
                         })
        if (verbose){
          if (length(levels(cv_test_y)) == 2){
            cat(paste0("AUC - ", levels(cv_test_y)[i], ": ", perf$auc, "\n"))
          } else {
            cat(paste0(" AUC - ", levels(cv_test_y)[i], " (vs Others): ", perf$auc, "\n"))
          }
        }

        perf
      } else {
        return(NULL)
      }
    }
    names(roc_auc_list) <- unique(cv_test_y)

    roc_dfm <- foreach(i = 1:length(levels(cv_test_y)), .combine = "rbind") %do% {
      perf <- roc_auc_list[[i]]
      fpr <- 1 - perf$specificities
      tpr <- perf$sensitivities
      thresholds <- perf$thresholds
      mtx <- cbind(fpr, tpr, thresholds)
      if (length(levels(cv_test_y)) == 2){
        df <- data.frame(mtx, group = rep(levels(cv_test_y)[i], times = nrow(mtx)), row.names = NULL, check.names = FALSE)
      } else {
        df <- data.frame(mtx, group = rep(paste0(levels(cv_test_y)[i], " (vs Others)"), times = nrow(mtx)), row.names = NULL, check.names = FALSE)
      }
      df <- df[order(df$tpr), ]
      return(df)
    }

    # return
    if (any(sapply(roc_auc_list, is.null))) {
      out <- NULL
    } else {
      out <- list(svm.roc_object = roc_auc_list,
                  svm.roc_dataframe = roc_dfm)
    }
    return(out)
  }

  # ---- compute and return ----
  auc_res_list <- vector(mode = "list", length = length(cv_model_list))
  auc_res_list[] <- foreach(i = 1:length(cv_model_list)) %do% {
    out <- cv_auc_func(cv_model_list[[i]])
    out
  }
  names(auc_res_list) <- names(cv_model_list)

  if (any(sapply(auc_res_list, is.null))){
    cat("NULL element removed from the CV ROC-AUC list.\n")
    auc_res_list <- auc_res_list[-which(sapply(auc_res_list, is.null))]
  }

  if (class(object) == 'rbiosvm_nestedcv') {  # export
    if (is.null(fileprefix)) {
      assign(paste(deparse(substitute(object)), "_svm_nestedcv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    } else {
      assign(paste(as.character(fileprefix), "_svm_nestedcv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    }
  } else {
    if (is.null(fileprefix)) {
      assign(paste(deparse(substitute(object)), "_svm_cv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    } else {
      assign(paste(as.character(fileprefix), "_svm_cv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    }
  }

  # ---- plotting ----
  if (rocplot){
    if (length(auc_res_list) < 1){
      cat("ROC data empty. No plots can be generated. \n")
    } else {
      plot_dfm <- foreach(i = 1:length(auc_res_list), .combine = "rbind") %do% {
        dfm <- auc_res_list[[i]]$svm.roc_dataframe
        dfm$cv_fold <- names(auc_res_list)[i]
        dfm
      }
      plot_dfm$cv_fold <- factor(plot_dfm$cv_fold, levels = unique(plot_dfm$cv_fold))

      for (i in 1:length(unique(plot_dfm$group))){
        if (verbose) cat(paste0("Plot being saved to file: ", deparse(substitute(object)),".cv_roc.", unique(plot_dfm$group)[i], ".pdf..."))  # initial message

        plot_dfm_g <- plot_dfm[plot_dfm$group %in% unique(plot_dfm$group)[i], ]

        plt <- ggplot(data = plot_dfm_g, aes(x = fpr, y = tpr, group = cv_fold, colour = cv_fold)) +
          geom_line(aes(linetype = cv_fold), linewidth = plot.lineSize) +
          geom_abline(intercept = 0) +
          ggtitle(ifelse(plot.display.Title, "ROC", NULL)) +
          xlab(plot.xLabel) +
          ylab(plot.yLabel)
        if (plot.rightsideY) {
          plt <- plt +
            scale_y_continuous(sec.axis = dup_axis()) +
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                  plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                  axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                  axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                  axis.title.y.right = element_blank(),
                  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                  legend.key = element_blank(),
                  axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                  axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
        } else {
          plt <- plt +
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                  plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                  axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                  axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                  legend.key = element_blank(),
                  axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                  axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
        }

        # # below not needed for ggplot 2 3.5.0
        # if (plot.rightsideY){
        #   plt <- RBioplot::rightside_y(plt)
        # }

        # save
        # grid.newpage()
        if (is.null(fileprefix)){
          ggsave(filename = paste0(deparse(substitute(object)),".cv_roc.", unique(plot_dfm$group)[i], ".pdf"), plot = plt,
                 width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
        } else {
          ggsave(filename = paste0(as.character(fileprefix),".cv_roc.", unique(plot_dfm$group)[i], ".pdf"), plot = plt,
                 width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
        }
        grid.draw(plt)
        if (verbose) cat("Done!\n")
      }
    }
  }
}


#' @title rbioClass_svm_cv_roc_auc_v2
#'
#' @description ROC-AUC analysis and ploting for SVM cross-validation (CV) models, for classification only.
#' @param object A \code{rbiosvm_nestedcv} or \code{rbiosvm_cv} object.
#' @param fileprefix String. A file prefix to use for export file name, instead of the object name. Default is \code{NULL}.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.smooth If to smooth the curves. Uses binomial method to smooth the curves. Default is \code{FALSE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Prints AUC values in the console. And a pdf file for ROC plot.
#'         The function also exports a ROC results list with an ROC object for each CV fold.
#'
#' @details
#'          This "v2" function updates the original function \code{rbioClass_svm_cv_roc_auc} with the \code{svm_cv_rocauc_helper} function integrated,
#'          thereby substantially reducing code redundancy.The "v2" function functions exactly the same as the original function.
#'
#'          Eventually, the original function will be replaced by the "v2" function after testing.
#'
#'          Uses pROC module to calculate ROC. The function supports more than two groups or more than one threshold
#'          for classification and regression model.
#'
#'          The test data included in the input \code{rbiosvm_nestedcv} is already normalized with training data information.
#'          So there is no need for additional data transformation.
#'
#'          Also: a ROC curve with AUC == 1 is always 1-1 for 95% CI and can be misleading.
#'
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioClass_svm_cv_roc_auc_v2(object = svm_nested_cv)
#' }
#' @export
rbioClass_svm_cv_roc_auc_v2 <- function(object, fileprefix = NULL,
                                        rocplot = TRUE,
                                        plot.smooth = FALSE,
                                        plot.lineSize = 1,
                                        plot.display.Title = TRUE, plot.titleSize = 10,
                                        plot.fontType = "sans",
                                        plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                        plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                        plot.legendSize = 9, plot.rightsideY = TRUE,
                                        plot.Width = 170, plot.Height = 150,
                                        verbose = TRUE){
  # ---- argements check ----
  if (!any(class(object) %in% c('rbiosvm_nestedcv', 'rbiosvm_cv'))) stop("object needs to be \"rbiosvm_nestedcv\" or \"rbiosvm_cv\" classes.")
  auc_res <- svm_cv_rocauc_helper(object = object, roc.smooth = plot.smooth, verbose = verbose)
  auc_res_list <- auc_res$cv_res_list

  if (auc_res$input_model_class == 'rbiosvm_nestedcv') {  # export
    if (is.null(fileprefix)) {
      assign(paste(deparse(substitute(object)), "_svm_nestedcv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    } else {
      assign(paste(as.character(fileprefix), "_svm_nestedcv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    }
  } else {
    if (is.null(fileprefix)) {
      assign(paste(deparse(substitute(object)), "_svm_cv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    } else {
      assign(paste(as.character(fileprefix), "_svm_cv_roc_auc", sep = ""), auc_res_list, envir = .GlobalEnv)
    }
  }

  # ---- plotting ----
  if (rocplot){
    if (length(auc_res_list) < 1){
      cat("ROC data empty. No plots can be generated. \n")
    } else {
      plot_dfm <- foreach(i = 1:length(auc_res_list), .combine = "rbind") %do% {
        dfm <- auc_res_list[[i]]$svm.roc_dataframe
        dfm$cv_fold <- names(auc_res_list)[i]
        dfm
      }
      plot_dfm$cv_fold <- factor(plot_dfm$cv_fold, levels = unique(plot_dfm$cv_fold))

      for (i in 1:length(unique(plot_dfm$group))){
        if (verbose) cat(paste0("Plot being saved to file: ", deparse(substitute(object)),".cv_roc.", unique(plot_dfm$group)[i], ".pdf..."))  # initial message

        plot_dfm_g <- plot_dfm[plot_dfm$group %in% unique(plot_dfm$group)[i], ]

        plt <- ggplot(data = plot_dfm_g, aes(x = fpr, y = tpr, group = cv_fold, colour = cv_fold)) +
          geom_line(aes(linetype = cv_fold), linewidth = plot.lineSize) +
          geom_abline(intercept = 0) +
          ggtitle(ifelse(plot.display.Title, "ROC", NULL)) +
          xlab(plot.xLabel) +
          ylab(plot.yLabel)
        if (plot.rightsideY) {
          plt <- plt +
            scale_y_continuous(sec.axis = dup_axis()) +
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                  plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                  axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                  axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                  axis.title.y.right = element_blank(),
                  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                  legend.key = element_blank(),
                  axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                  axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
        } else {
          plt <- plt +
            theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
                  plot.title = element_text(face = "bold", size = plot.titleSize, family = plot.fontType, hjust = 0.5),
                  axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
                  axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
                  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.legendSize),
                  legend.key = element_blank(),
                  axis.text.x = element_text(size = plot.xTickLblSize, family = plot.fontType),
                  axis.text.y = element_text(size = plot.yTickLblSize, family = plot.fontType, hjust = 0.5))
        }

        if (is.null(fileprefix)){
          ggsave(filename = paste0(deparse(substitute(object)),".cv_roc.", unique(plot_dfm$group)[i], ".pdf"), plot = plt,
                 width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
        } else {
          ggsave(filename = paste0(as.character(fileprefix),".cv_roc.", unique(plot_dfm$group)[i], ".pdf"), plot = plt,
                 width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
        }
        grid.draw(plt)
        if (verbose) cat("Done!\n")
      }
    }
  }
}


#' @title rbioClass_svm_cv_roc_auc_mean
#'
#' @description CV ROC-AUC mean analysis and ploting for SVM model, for classification model only.
#' @param object A \code{cv_auc_res} object.
#' @param fileprefix String. A file prefix to use for export file name, instead of the objecte name. Default is \code{NULL}.
#' @param roc.smooth If to smooth the curves. Uses binomial method to smooth the curves. Default is \code{FALSE}.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.rightsideY If to show the right side y-axis. Default is \code{FALSE}.
#' @param plot.Width Scoreplot width. Default is \code{170}.
#' @param plot.Height Scoreplot height. Default is \code{150}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @returns
#'        1. PDF of the mean ROC-AUC plot
#'        2. A \code{CV_rocauc_mean} object to the environment
#'
#' @details
#'        The function construct a  list containing data.frames of AUC results over folds per group
#'        for every fold: list = roc dfm w cols: tpr, fpr, folds, auc dfm w cols: auc, folds
#'
#'        The function also uses interpolation for CV ROC: set fpr seq(0, 1, length.out = 100), and interpolated tpr,
#'        to make sure all ROCs from different CV classifiers to have the same length
#'        interpolation for ROC thus is also used to plot ROCs for different types classifiers for comparison
#'
#'        The ROC-AUC calcualte details follows the same as \code{\link{rbioClass_svm_roc_auc}}
#'
#'        \code{v_rocauc_mean} items:
#'          1. \code{cv_auc_mean_res}: mean CV ROC-AUC calculation results organized by outcome (or \code{y}) groups
#'          2. \code{roc_plot}: mean ROC-AUC plot
#'
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom matrixStats rowMeans2 rowSds
#'
#' @export
rbioClass_svm_cv_roc_auc_mean <- function(object, fileprefix = NULL,
                                          roc.smooth = FALSE,
                                          rocplot = TRUE,
                                          plot.lineSize = 1,
                                          plot.display.Title = TRUE, plot.titleSize = 10,
                                          plot.fontType = "sans",
                                          plot.xLabel = "1 - specificity", plot.xLabelSize = 15, plot.xTickLblSize = 10,
                                          plot.yLabel = "sensitivity", plot.yLabelSize = 15, plot.yTickLblSize = 10,
                                          plot.legendSize = 9, plot.rightsideY = TRUE,
                                          plot.Width = 170, plot.Height = 150,
                                          verbose = TRUE) {
  # ---- check the class ----
  # if (!'cv_auc_res' %in% class(object)) stop("object needs to be \"cv_auc_res\" classes.")
  if (!any(class(object) %in% c('rbiosvm_nestedcv', 'rbiosvm_cv'))) stop("object needs to be \"rbiosvm_nestedcv\" or \"rbiosvm_cv\" classes.")
  auc_res <- svm_cv_rocauc_helper(object = object, roc.smooth = roc.smooth, verbose = verbose)
  auc_res_list <- auc_res$cv_res_list

  # ---- CV ROC-AUC curve ----
  auc_res_group <- unique(auc_res_list$cv_fold_1$svm.roc_dataframe$group)

  # construct a  list containing data.frames of AUC results over folds per group
  # for every fold: list = roc dfm w cols: tpr, fpr, folds, auc dfm w cols: auc, folds
  # use interpolation for CV ROC: set fpr seq(0, 1, length.out = 100), and interpolated tpr
  # to make sure all ROCs from different CV classifiers to
  # have the same length
  # interpolation for ROC thus is also used to plot ROCs for different types classifiers
  # for comparison
  fpr_mean <- seq(0, 1, length.out = 100)
  auc_group_list <- vector(mode = "list", length = length(auc_res_group))

  for (i in 1:length(auc_res_group)){
    group_label <- auc_res_group[[i]]
    group <- gsub(" (vs Others)", "", auc_res_group[[i]], fixed = TRUE)
    if (verbose) cat(paste0("processing group: ", group, ", group label: ",  group_label, "\n"))
    cv_tpr_list <- vector(mode = "list",  length = length(auc_res_list))
    names(cv_tpr_list) <- paste0("cv_fold_", 1:length(cv_tpr_list))
    for (j in 1:length(cv_tpr_list)) {
      roc <- auc_res_list[[j]]$svm.roc_dataframe[auc_res_list[[j]]$svm.roc_dataframe$group %in% group_label, ]
      # cv_tpr_list[[paste0("cv_fold_", j)]] <- approx(roc$fpr, roc$tpr, xout = fpr_mean, ties = "mean")$y
      if (nrow(roc) == 0) {
        cv_tpr_list[[paste0("cv_fold_", j)]] <- NULL
      } else {
        cv_tpr_list[[paste0("cv_fold_", j)]] <- approx(roc$fpr, roc$tpr, xout = fpr_mean, ties = "mean")$y
      }
    }
    cv_tpr <- do.call(cbind, cv_tpr_list)
    cv_fpr <- fpr_mean

    cv_auc <- foreach(m = 1:length(auc_res_list), .combine = "rbind") %do% {
      tryCatch(
        dfm <- data.frame(auc = as.numeric(auc_res_list[[m]]$svm.roc_object[[group]]$auc), cv_fold = names(auc_res_list)[m]),
        error = function(e) {
          dfm <- NULL
        }
      )
    }
    auc_group_list[[i]] <- list(cv_tpr = cv_tpr, cv_fpr = cv_fpr, cv_auc = cv_auc)
    # break
  }
  names(auc_group_list) <- auc_res_group

  # calculate mean tpr and mean auc
  cv_roc_mean_list <- vector(mode = "list", length = length(auc_group_list))
  for (i in 1:length(cv_roc_mean_list)) {
    group_label <- auc_res_group[[i]]
    group <- gsub(" (vs Others)", "", auc_res_group[[i]], fixed = TRUE)
    if (verbose) cat(paste0("processing group: ", group, ", group label: ",  group_label, "\n"))

    mean_tpr <- rowMeans2(auc_group_list[[group_label]]$cv_tpr)
    mean_tpr[length(mean_tpr)] <- 1.0
    sd_tpr <- rowSds(auc_group_list[[group_label]]$cv_tpr)
    tprs_upper <- min(mean_tpr + sd_tpr, 1)
    tprs_lower <- max(mean_tpr - sd_tpr, 0)

    mean_auc <- mean(auc_group_list[[group_label]]$cv_auc$auc)
    sd_auc <- sd(auc_group_list[[group_label]]$cv_auc$auc)

    cv_roc_mean_list[[i]] <- list(
      mean_fpr = fpr_mean,
      mean_tpr = rowMeans2(auc_group_list[[group_label]]$cv_tpr),
      sd_tpr = rowSds(auc_group_list[[group_label]]$cv_tpr),
      tprs_upper = pmin(rowMeans2(auc_group_list[[group_label]]$cv_tpr) + rowSds(auc_group_list[[group_label]]$cv_tpr), 1),
      tprs_lower = pmax(rowMeans2(auc_group_list[[group_label]]$cv_tpr) - rowSds(auc_group_list[[group_label]]$cv_tpr), 0),
      mean_auc = mean(auc_group_list[[group_label]]$cv_auc$auc),
      sd_auc = sd(auc_group_list[[group_label]]$cv_auc$auc)
    )
  }
  names(cv_roc_mean_list) <- names(auc_group_list)

  # ------ plotting ------
  if (rocplot) {
    plot_dfm <- foreach(i = 1:length(cv_roc_mean_list), .combine = "rbind") %do% {
      init_d <- data.frame(mean_fpr = 0, mean_tpr = 0,
                           tprs_upper = 0, tprs_lower = 0)
      end_d <- data.frame(mean_fpr = 1, mean_tpr = 1,
                          tprs_upper = 1, tprs_lower = 1)
      d <- data.frame(mean_fpr = cv_roc_mean_list[[i]]$mean_fpr, mean_tpr = cv_roc_mean_list[[i]]$mean_tpr,
                      tprs_upper = cv_roc_mean_list[[i]]$tprs_upper, tprs_lower = cv_roc_mean_list[[i]]$tprs_lower)
      d <- rbind(init_d, d)
      d <- rbind(d, end_d)
      d$group <- as.factor(rep(names(cv_roc_mean_list)[i], times = nrow(d)))
      d
    }

    if (length(unique(plot_dfm$group)) < 2) {
      p <- ggplot(data = plot_dfm, aes(x = .data[["mean_fpr"]], y = .data[["mean_tpr"]], colour = group, group = group, shape = group)) +
        geom_point(size = 2) +
        geom_line() +
        geom_ribbon(aes(ymin = .data[["tprs_lower"]], ymax = .data[["tprs_upper"]],
                        fill = group), alpha = 0.15, colour = NA) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, NA)) +
        scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA), sec.axis = dup_axis()) +
        ggtitle("CV ROC-AUC") +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        theme_bw() +
        theme(legend.key = element_blank(),
              plot.margin=unit(c(1,0.5,1,0.5),"cm"),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              axis.title.y.right = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none",
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5))
    } else {
      p <- ggplot(data = plot_dfm, aes(x = .data[["mean_fpr"]], y = .data[["mean_tpr"]], colour = group, group = group, shape = group)) +
        geom_point(size = 2) +
        geom_line() +
        geom_ribbon(aes(ymin = .data[["tprs_lower"]], ymax = .data[["tprs_upper"]],
                        fill = group), alpha = 0.15, colour = NA) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, NA)) +
        scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, NA), sec.axis = dup_axis()) +
        ggtitle("CV ROC-AUC") +
        xlab(plot.xLabel) +
        ylab(plot.yLabel) +
        theme_bw() +
        theme(legend.key = element_blank(),
              plot.margin=unit(c(1,0.5,1,0.5),"cm"),
              axis.title.x = element_text(face = "bold", size = plot.xLabelSize, family = plot.fontType),
              axis.title.y = element_text(face = "bold", size = plot.yLabelSize, family = plot.fontType),
              axis.title.y.right = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5))
    }

    # -- export plot --
    if (is.null(fileprefix)){
      ggsave(filename = paste0(deparse(substitute(object)),".cv_roc_mean.pdf"), plot = p,
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    } else {
      ggsave(filename = paste0(as.character(fileprefix),".cv_roc_mean.pdf"), plot = p,
             width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    }
    grid.draw(p)
  } else {
    p <- NULL
  }

  # ----- export results to the environment ----
  out_list <- vector(mode = "list", length = 2)
  names(out_list) <- c("cv_auc_mean_res", "roc_plot")

  out_list$cv_auc_mean_res <- auc_group_list
  for (i in 1:length(out_list$cv_auc_mean_res)) {
    # cat(paste0("processing...", i))
    out_list$cv_auc_mean_res[[names(out_list$cv_auc_mean_res)[i]]]$mean_auc <- cv_roc_mean_list[[names(cv_roc_mean_list)[i]]]$mean_auc
    out_list$cv_auc_mean_res[[names(out_list$cv_auc_mean_res)[i]]]$sd_auc <- cv_roc_mean_list[[names(cv_roc_mean_list)[i]]]$sd_auc
  }
  out_list$roc_plot <- p

  class(out_list) <- "cv_rocauc_mean"

  if (auc_res$input_model_class == 'rbiosvm_nestedcv') {  # export
    if (is.null(fileprefix)) {
      assign(paste(deparse(substitute(object)), "_svm_nestedcv_rocauc_mean", sep = ""), out_list, envir = .GlobalEnv)
    } else {
      assign(paste(as.character(fileprefix), "_svm_nestedcv_rocauc_mean", sep = ""), out_list, envir = .GlobalEnv)
    }
  } else {
    if (is.null(fileprefix)) {
      assign(paste(deparse(substitute(object)), "_svm_cv_rocauc_mean", sep = ""), out_list, envir = .GlobalEnv)
    } else {
      assign(paste(as.character(fileprefix), "_svm_cv_rocauc_mean", sep = ""), out_list, envir = .GlobalEnv)
    }
  }
}


#' @title rbioClass_svm_perm()
#'
#' @description Permutation test for SVM models.
#' @param object A \code{rbiosvm} object. Make sure the object is generated with a \code{tot.accuracy} section.
#' @param perm.method Permutation method. Options are \code{"by_y"} and \code{"by_feature_per_y"}. Default is \code{"by_y"}. See details below.
#' @param nperm Number of permutations to run. Default is \code{999}.
#' @param perm.plot whether to produce a plot or not. Default is \code{TRUE}.
#' @param ... Additional argument for \code{\link{rbioUtil_perm_plot}}.
#' @param parallelComputing whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return The function returns \code{CSV} files for all intermediate permutation accuracy values as well as the p-value resutls.
#'
#'
#' The final results are also exported to the environment as a \code{rbiosvm_perm} object with the following items:
#'
#' \code{perm.method} Permutation method.
#'
#' \code{nperm} The number of permutation runs.
#'
#' \code{performance.type} The stats metric used for the permutation test.
#'
#' \code{original.performance} The original accuracy for the SVM model.
#'
#' \code{perm.results} The intermediate permutation results, i.e. stats for each permutation test run in a data.frame. \code{nperm = 0} is the original stats.
#'
#' \code{p.value} P value for the permutation test.
#'
#' \code{model.type} SVM model type.
#'
#' \code{run.time} run time
#'
#' A scatter plot is also generaeted when \code{perm.plot = TRUE}.
#'
#' @details The function uses RMSEP as the stats for comparing original model with permutatsions.
#'
#' Data for permutation are object$centerX$centerX, meaning centered X are used in applicable.
#'
#' Permutation methods are according to:
#'
#' Ojala M, Garriga GC. 2010. Permutation test for studying classifier performance. J Mach Learn Res. 11: 1833 - 63.
#'
#' For \code{perm.method = "by_y"}, labels (i.e. y) are permutatedted. A non-signifianct model (permutation p value > alpha, i.e. 0.05) in this case means the data is independent from the groups.
#'
#' For \code{perm.method = "by_feature_per_by"}, X is first subset by label (i.e.y) before permutating data for each feature, i.e. by column. Since the permutation is done for the features WITHIN the group,
#' the test actually evaluates if the model will produce significantly different performannce from the permutation models with the original "betweeen-features" relation (if any) disturbed.
#' Therefore, A non-significant result (permutation p value > alpha, i.e. 0.05) means either the features are independent, or the model doesn't consider correlation between the features.
#'
#' @import ggplot2
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' rbioClass_plsda_perm()
#' }
#' @export
rbioClass_svm_perm <- function(object,
                               perm.method = c("by_y", "by_feature_per_y"), nperm = 999,
                               perm.plot = TRUE, ...,
                               parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                               verbose = TRUE){
  ## initiate the run time
  start_time <- Sys.time()

  ## check arguments
  if (!any(class(object) %in% c("rbiosvm"))) stop("object has to be a \"rbiosvm\" class.\n")
  if (object$model.type == "regression") {
    if (!"tot.MSE" %in% names(object) || is.null(object$tot.MSE)) stop("Regression SVM model object has to include tot.accuracy value from Cross-Validation for permutation test. \n")
  } else {
    if (!"tot.accuracy" %in% names(object) || is.null(object$tot.accuracy)) stop("Classification SVM model object has to include tot.accuracy value from Cross-Validation for permutation test. \n")
  }
  # if (!perm.method %in% c("by_y", "by_feature_per_y")) stop("perm.method needs to be either \"by_y\" or \"by_feature_per_y\". \n")
  if (length(nperm) != 1) stop("nperm can only contain one integer. \n")
  if (nperm %% 1 != 0) stop("nperm can only be integer. \n")
  if (nperm < 1) stop("nperm can only take interger equal to or greater than 1. \n")
  perm.method <- match.arg(tolower(perm.method), c("by_y", "by_feature_per_y"))
  if (object$model.type == "regression" && perm.method == "by_feature_per_y") stop("perm.method cannot be \"by_feature_per_y\" when SVM model type is regression.")

  if (parallelComputing){
    clusterType <- match.arg(clusterType, c("PSOCK", "FORK"))
  }

  ## Perm test
  # calcuate permutation error and construct original error data frame
  if (object$model.type == "regression"){
    orig_perfm <- object$tot.MSE
    orig_perfm <- sqrt(orig_perfm)  # tot.RMSE
  } else {
    orig_perfm <- object$tot.accuracy
    # orig_accu <- object$tot.accuracy
  }

  # perm functions
  by_y_func <- function(i){
    set.seed(i)
    perm_y <- object$inputY[sample(1:length(object$inputY))]  # sample label permutation
    if (is.null(object$center.scaledX$centerX)) {
      perm_model <- rbioClass_svm(x = object$inputX, y = perm_y,
                                  center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                  tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                  verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred

    } else {
      perm_model <- rbioClass_svm(x = object$center.scaledX$centerX, y = perm_y,
                                  center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                  tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                  verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
    }

    if (object$model.type == "regression"){
      perm_perfm <- perm_model$tot.MSE
      perm_perfm <- sqrt(perm_perfm)  # RMSE
    } else {
      perm_perfm <- perm_model$tot.accuracy  # permutation model accuracy
    }
    perm_perfm_dfm <- data.frame(nperm = i, stats = perm_perfm, row.names = NULL)
    return(perm_perfm_dfm)
  }
  by_feature_per_y_func <- function(i){
    set.seed(i)
    perm_x <- foreach(m = unique(levels(object$inputY)), .combine = "rbind") %do% {
      if (is.null(object$center.scaledX$centerX)){
        sub_dat <- object$inputX[which(object$inputY == m), ]  # subsetting the centre-scaled X by label (Y)
      } else {
        sub_dat <- object$center.scaledX$centerX[which(object$inputY == m), ]  # subsetting the centre-scaled X by label (Y)
      }

      sub_dat_colnames <- colnames(sub_dat)
      perm_dat <- sub_dat[, sample(1:ncol(sub_dat))]
      # perm_dat <- foreach(n = 1:ncol(sub_dat), .combine = "cbind") %do% {
      #   col <- sub_dat[sample(1:nrow(sub_dat)), n, drop = FALSE]  # permutation for each feature
      #   col
      # }
      colnames(perm_dat) <- sub_dat_colnames
      perm_dat
    }

    perm_model <- rbioClass_svm(x = perm_x, y = object$inputY,
                                center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred

    if (object$model.type == "regression"){
      perm_perfm <- perm_model$tot.MSE
      perm_perfm <- sqrt(perm_perfm)  # RMSE
    } else {
      perm_perfm <- perm_model$tot.accuracy  # permutation model accuracy
    }
    perm_perfm_dfm <- data.frame(nperm = i, stats = perm_perfm, row.names = NULL)
    return(perm_perfm_dfm)
  }

  # permutation test: consolidated permutation model error data frame
  if (verbose){
    cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
    cat(paste0("Running permutation test using ", perm.method, " method ","with ", nperm, " permutations (speed depending on hardware configurations)..."))
  }

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
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %dopar% by_y_func(i)
    } else {  # subset by label first, then permutate data per feature
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %dopar% by_feature_per_y_func(i)
    }
  }

  perm_perfm_val <- perm_dfm$stats
  if (object$model.type == "regression"){
    p.val <- (length(which(perm_perfm_val <= orig_perfm)) + 1) / (nperm + 1)  # MSE the smaller the better
  } else {
    p.val <- (length(which(perm_perfm_val >= orig_perfm)) + 1) / (nperm + 1)  # accuracy the bigger the better
  }
  if (verbose) cat("Done!\n")

  ## output
  # end time
  end_time <- Sys.time()
  runtime <- end_time - start_time

  # output
  if (object$model.type == "regression"){
    perfm_type <- "tot.RMSE"
    names(perm_dfm)[2] <- "tot.RMSE"
  } else {
    perfm_type <- "tot.accuracy"
    names(perm_dfm)[2] <- "tot.accuracy"
  }
  complete_perfm <- rbind(c(0, orig_perfm), perm_dfm)
  out <- list(perm.method = perm.method,
              nperm = nperm,
              performance.type = perfm_type,
              original.performance = orig_perfm,
              perm.results = complete_perfm,
              p.value = p.val,
              model.type = object$model.type,
              run.time = paste0(signif(runtime[[1]], 4), " ", attributes(runtime)[2]))
  class(out) <- "rbiosvm_perm"
  assign(paste(deparse(substitute(object)), "_perm", sep = ""), out, envir = .GlobalEnv)

  # export to the directory
  if (verbose) cat(paste0("Permutation stats results stored in csv file: ", paste(deparse(substitute(object)), "_perm.csv", sep = ""), ". \n"))
  write.csv(file = paste0(deparse(substitute(object)), ".perm.csv"), complete_perfm, row.names = FALSE)

  ## plot
  if (perm.plot){
    svm_permutation_test <- out
    rbioUtil_perm_plot(perm_res = svm_permutation_test, ...)
  }

  # print results
  if (verbose){
    cat("\n")
    cat("Permutation test restuls (p value):")
    cat(p.val)
    cat("\n")
  }
}


#' @export
print.rbiosvm_perm <- function(x, ...){
  cat(paste0("SVM permutation results with ", x$nperm, " permutations:\n"))
  cat("\n")
  cat("p-value: \n")
  print(x$p.value)
  cat("\n")
}


#' @title rbioReg_svm_rmse
#'
#' @description Support Vector Regression (SVR) RMSE calculation
#' @param object A \code{rbiosvm} object.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.y For regression model only, the vector for the new data's continuous outcome variable. Default is \code{NULL}
#' @return RMSE value with either new data or training data.
#'
#' @details
#'
#' With newdata, the newdata is center_scaled using training data parameter (SD etc).
#' Without newdata, the RMSE is based on the CV RMSE.
#'
#' @examples
#' \dontrun{
#'  test_rmse <- rbioReg_svm_rmse(object = svm_m, newdata = svm_test[,-1], newdata.y = svm_test$y)
#' }
#' @export
rbioReg_svm_rmse <- function(object, newdata=NULL, newdata.y=NULL){
  # argument check
  if (!any(class(object) %in% "rbiosvm")) stop("The input object needs to be a \"rbiosvm\" class.")
  if (object$model.type != "regression") stop("The input model type needs to be \"regression\".")
  if (!is.null(newdata)) {
    if (is.null(newdata.y)) stop("newdata.y needs to be specified if newdata is available.")
    if (nrow(newdata) != length(newdata.y)) stop("sample size needs to match the length of newdata.y.")
  }
  if (!any(object$tune.method %in% c("cross", "boot"))) stop("The tune.method should be either \"cross\" or \"boot\".")

  # computation
  if (is.null(newdata)){
    out_rmse <- sqrt(object$tot.MSE)
  } else {
    center_scale_newdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
    pred <- predict(object, center_scale_newdata)
    err <- newdata.y - pred
    out_rmse <- sqrt(mean(err^2))
  }

  return(out_rmse)
}


#' @title rbioClass_svm_shap_aggregated
#' @description Aggregated SHAP analysis for SVM model.
#' @param model Input SVM model object and should be a \code{rbiosvm} class.
#' @param X Input data, matrix or  data.frame, with column: features, and row: records. It should not contain outcome column.
#' @param bg_X Background data, \code{matrix} or \code{data.frame}, column: features, and row: records. It should not contain outcome column.
#' @param bg_n Number of sample to draw from the background data.
#' @param parallelComputing whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param plot If to generate aggregated SHAP plot.
#' @param plot.filename.prefix Prefix string to add to the output plot file.
#' @param plot.type The type of the SHAP plot, options are \code{"both", "bee", "bar"}.
#' @param plot.n Number of features to display in the SHAP plot.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.Width Plot width, in \code{mm} unit. Default is \code{15}
#' @param plot.Height Plot height, in \code{mm} unit.. Default is \code{10}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details
#'    1. This function relies on \code{kernelshap} and \code{shapviz} logics.
#'    2. (To be tested) The function should work with regression SVM models as well.
#'    3. When set \code{bg_x = NULL} (default), the function uses the input \code{X} as the background data.
#'    4. The output plot file is suffixed with "_shap_aggregated_plot.pdf".
#' @return
#'    The function outputs a \code{rbio_shap} class object.
#' @import shapviz
#' @importFrom kernelshap kernelshap
#' @importFrom ggpubr ggarrange
#' @export
rbioClass_svm_shap_aggregated <- function(model, X, bg_X = NULL, bg_n = 200L,
                                          parallelComputing = FALSE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                                          plot = TRUE, plot.filename.prefix = NULL,
                                          plot.type = c("both", "bee", "bar"), plot.n = 15L,
                                          plot.SymbolSize = 2, plot.lineSize = 1,
                                          plot.display.Title = TRUE, plot.titleSize = 10,
                                          plot.fontType = "sans",
                                          plot.bee.colorscale = "A", plot.bar.color = "blue",
                                          plot.xLabel = "SHAP value", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                          plot.yLabel = "Features", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                          plot.legendPosition = "right", plot.legendSize = 9,
                                          plot.Width = 170, plot.Height = 150,
                                          verbose = TRUE) {
  # --- arg check ---
  if (!any(class(model) %in% c("rbiosvm"))) stop("object has to be a \"rbiosvm\" class.\n")

  # --- shap calculation ---
  if (parallelComputing) {
    pc <- TRUE
    # set up cpu cluster
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function
  } else {
    pc <- FALSE
  }

  if (model$model.type == "classification") {
    y_labels <- as.character(unique(model$inputY))
    ks_list <- vector(mode = "list", length = length(y_labels))

    for (i in 1:length(y_labels)) {
      l <- y_labels[i]
      ks_list[[i]] <- kernelshap(model, X = X,
                                 bg_X = bg_X,
                                 bg_n = bg_n,
                                 pred_fun = function(model, X) rbio_shap_svm_label_prob(model, X, col_idx = l),
                                 parallel = pc,
                                 verbose = verbose)
      names(ks_list)[i] <- l
      # break
    }
  } else {
    y_labels <- model$inputY
    ks_list <- vector(mode = "list", length = 1)
    ks_list[[1]] <- kernelshap(model, X = X,
                               bg_X = bg_X,
                               pred_fun = as.numeric(e1071:::predict.svm(model, X)),
                               parallel = pc,
                               verbose = verbose) # needed to be updated with prediction scaling

    names(ks_list) <- "y"
  }

  # --- plotting ---
  if (plot) {
    g_list <- vector(mode = "list", length = length(ks_list))
    opt <- list(begin = 0.05, end = 0.85, option = plot.bee.colorscale)
    for (i in 1:length(ks_list)) {
      ks <- ks_list[[i]]
      s <- shapviz(ks)
      g <- sv_importance_shapviz(s, kind = plot.type, max_display = plot.n,
                                 show_numbers = TRUE, viridis_args = opt, fill = "green")
      g <- shap_g_theme_helper(g,
                               plot.SymbolSize = plot.SymbolSize, plot.lineSize = plot.lineSize,
                               plot.display.Title = plot.display.Title, plot.titleSize = plot.titleSize,
                               plot.title = paste0("Class: ", names(ks_list[i])),
                               plot.fontType = plot.fontType,
                               plot.xLabel = plot.xLabel, plot.xLabelSize = plot.xLabelSize, plot.xTickLblSize = plot.xTickLblSize,
                               plot.yLabel = plot.yLabel, plot.yLabelSize = plot.yLabelSize, plot.yTickLblSize = plot.yTickLblSize,
                               plot.legendPosition = plot.legendPosition, plot.legendSize = plot.legendSize,
                               plot.Width = plot.Width, plot.Height = plot.Height)
      g_list[[i]] <- g
      names(g_list)[i] <- names(ks_list)[i]
    }
    g_aggreg <- ggarrange(plotlist = g_list,
                          common.legend = TRUE, legend = "right")
    plot_filename <- paste0(plot.filename.prefix, "_shap_aggregated_plot.pdf")
    ggsave(plot = g_aggreg, filename = plot_filename,
           device = "pdf", units = "mm",
           width = plot.Width, height = plot.Height, dpi= 600)
    if (verbose) cat(paste0("\nplot saved to file: ", plot_filename))
  } else {
    g_list <- NULL
  }

  # --- output ---
  o <- list(
    shap_ks = ks_list,
    shap_plot = g_list,
    shap_type = "aggregated",
    plot_type = plot.type
  )
  class(o) <- "rbio_shap"
  return(o)
}


#' @title rbioClass_svm_shap_individual
#' @description Individualized SHAP analysis for SVM model.
#' @param model Input SVM model object and should be a \code{rbiosvm} class.
#' @param X Input data, matrix or  data.frame, with column: features, and row: records. It should not contain outcome column.
#' @param bg_X Background data, \code{matrix} or \code{data.frame}, column: features, and row: records. It should not contain outcome column.
#' @param bg_n Number of sample to draw from the background data.
#' @param parallelComputing whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param plot If to generate aggregated SHAP plot.
#' @param plot.filename.prefix Prefix string to add to the output plot file.
#' @param plot.type The type of the individual SHAP plot, options are \code{"waterfall", "force"}.
#' @param plot.n Number of features to display in the SHAP plot.
#' @param ... Additional arguments for the \code{shapviz::sv_waterfall} or \code{shapviz::sv_force} functions.
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param plot.Width Plot width, in \code{mm} unit. Default is \code{15}
#' @param plot.Height Plot height, in \code{mm} unit.. Default is \code{10}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @details
#'    1. This fucntion relies on \code{kernelshap} and \code{shapviz} logic.
#'    2. (To be tested) The function should work with regression SVM models as well.
#' @return
#'    The function outputs a \code{rbio_shap} class object.
#' @import shapviz
#' @importFrom kernelshap kernelshap
#' @importFrom ggpubr ggarrange
#' @export
rbioClass_svm_shap_individual <- function(model, X, bg_X, bg_n = 200L,
                                          y = NULL,
                                          parallelComputing = FALSE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                                          plot.filename.prefix = NULL,
                                          plot.type = c("waterfall", "force"), plot.n = 15L, ...,
                                          plot.SymbolSize = 2, plot.lineSize = 1,
                                          plot.display.Title = TRUE, plot.titleSize = 10,
                                          plot.fontType = "sans",
                                          plot.bee.colorscale = "A", plot.bar.color = "blue",
                                          plot.xLabel = "Prediction value", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                          plot.yLabel = "Features", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                          plot.legendPosition = "right", plot.legendSize = 9,
                                          plot.Width = 170, plot.Height = 150,
                                          verbose = TRUE) {
  # --- arg check ---
  if (!any(class(model) %in% c("rbiosvm"))) stop("object has to be a \"rbiosvm\" class.\n")
  plot.type <- match.arg(plot.type)

  # --- shap calculateion ---
  if (parallelComputing) {
    pc <- TRUE
    # set up cpu cluster
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function
  } else {
    pc <- FALSE
  }

  if (model$model.type == "classification") {
    p <- rbioClass_svm_predict(model, newdata = X, verbose = verbose)
    l <- as.character(p$probability.summary$Class[which.max(p$probability.summary$Probability)])
    ks <- kernelshap(model, X = X,
                     bg_X = bg_X,
                     bg_n = bg_n,
                     pred_fun = function(model, X) RBioFS:::rbio_shap_svm_label_prob(model, X, col_idx = l),
                     parallel = pc,
                     verbose = verbose)
  } else {
    l <- as.numeric(e1071:::predict.svm(model, X))  # needed to be updated with scaling
    ks <- kernelshap(model, X = X,
                     bg_X = bg_X,
                     bg_n = bg_n,
                     pred_fun = as.numeric(e1071:::predict.svm(model, X)),
                     parallel = pc,
                     verbose = verbose)
  }

  if (!is.null(y)) {
    y <- as.character(y)
    if (l != y) {
      warning(paste0("prediction: ", l, " failed to match the input y: ", y, "\n"))
    }
    plot_title <- paste0("Sample prediction: ", l, " (input y: ", y, ")")
  } else {
    plot_title <- paste0("Sample prediction: ", l)
  }

  # --- plotting ---
  s <- shapviz(ks)
  if (plot.type == "waterfall") {
    g <- sv_waterfall(s, max_display = plot.n, ...)
  } else {
    g <- sv_force(s, max_display = plot.n, ...)
  }

  switch(plot.type,
         waterfall = {g <- sv_waterfall(s, max_display = plot.n, ...)},
         force = {g <- sv_force(s, max_display = plot.n, ...)}
  )

  g <- RBioFS:::shap_g_theme_helper(g,
                                    plot.SymbolSize = plot.SymbolSize, plot.lineSize = plot.lineSize,
                                    plot.display.Title = plot.display.Title, plot.titleSize = plot.titleSize,
                                    plot.title = plot_title,
                                    plot.fontType = plot.fontType,
                                    plot.xLabel = plot.xLabel, plot.xLabelSize = plot.xLabelSize, plot.xTickLblSize = plot.xTickLblSize,
                                    plot.yLabel = plot.yLabel, plot.yLabelSize = plot.yLabelSize, plot.yTickLblSize = plot.yTickLblSize,
                                    plot.legendPosition = plot.legendPosition, plot.legendSize = plot.legendSize,
                                    plot.Width = plot.Width, plot.Height = plot.Height)
  plot_filename <- paste0(plot.filename.prefix, "_shap_individual_waterfall.pdf")
  ggsave(plot = g, filename = plot_filename,
         device = "pdf", units = "mm",
         width = plot.Width, height = plot.Height, dpi= 600)
  if (verbose) cat(paste0("\nplot saved to file: ", plot_filename))

  # --- output ---
  o <- list(
    shap_ks = ks,
    shap_plot = g,
    type = "individual",
    plot_type = plot.type
  )
  class(o) <- "rbio_shap"
  return(o)
}
