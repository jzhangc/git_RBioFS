#' @title rbioFS_svm
#'
#' @description Support Vector Machine (SVM) modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param center.scale Logical, wether center and scale the data, i.e. subtracting mean and deviding by standard deviation. Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param svm.cross.k Fold of cross validation. Default is \code{10}.
#' @param tune.method Parameter tuning method. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{svm} function from \code{e1071} pacakge.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will be affect error or warning messeages.
#' @return Returns a SVM model object, with classes "svm" and "rbiosvm".
#'
#' Additional items for \code{rbiosvm} object to \code{svm} object from e1071 package:
#'
#' \code{inputX}: raw input predictor data.
#'
#' \code{inputY}: input group labels.
#'
#' \code{center.scaledX}: centered X data with scaling if applicable.
#'
#' \code{class.weight}: class weight, \code{1} if data is balanced.
#'
#' @details Parameter tuning is for gamma (not applicable when \code{kernel = "linear"}) and cost.
#'
#' The function automatically detects if the data is unbalanced or not, and applies class weight accordingly. For unbalanced data, the
#' class weight the inverse of the draw probability: \code{sample size / class size}.
#'
#' When \code{tune.method = "cross"} and the sample size set for \code{tune.cross.k}, the method is effectively "leave-one-out (LOO)" cross-validation.
#' Same goes for \code{svm.cross.k} argument.
#'
#' @importFrom e1071 svm tune tune.control
#' @importFrom matrixStats colSds
#' @examples
#' \dontrun{
#' svm_model <- rbioFS_svm(x = training_set[, -1], y = training_set[, 1], kernel = "radial", center = TRUE, scale = FALSE)
#' }
#' @export
rbioFS_svm <- function(x, y, center.scale = TRUE,
                       kernel = "radial", svm.cross.k = 10,
                       tune.method = "cross", tune.cross.k = 10, tune.boot.n = 10, ...,
                       verbose = TRUE){
  ## check arguments
  if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.")
  if (class(x) == "data.frame"){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    x <- as.matrix(x)
  }
  if (class(y) != "factor"){
    if (verbose) cat("y is converted to factor. \n")
    y <- factor(y, levels = unique(y))
  }
  if (!kernel %in% c("radial", "linear", "polynomial", "sigmoid")) stop("kernel needs to be exactly one of \"radial\", \"linear\", \"polynomial\", or \"sigmoid\". \n")

  ## data processing
  if (center.scale){
    if (verbose) cat(paste0("Data centered with the scaling prior to modelling..."))
    centered_X <- center_scale(x, scale = TRUE)  # center data with the option of scaling
    if (verbose) cat("DONE!\n")
    X <- centered_X$centerX
  } else {
    X <- x
    centered_X <- NULL
  }
  y <- y

  ## weight evaluation
  if (length(unique(table(y))) != 1){  # test if the sample is balanced
    wgt <- length(y) / table(y)
  } else {  # balanced sample, each class weight is 1
    wgt <- unique(table(y)) / table(y)
  }

  ## svm
  # tune parameters
  gamma_start <- ifelse(is.vector(X), 1, 1 / ncol(X))  # starting gamma value per svm() default settings
  if (verbose) cat(paste0("Grid searching for parameter optimization with ", tune.method, " method (speed depending on hardware configuration)..."))
  svm_tuned <- tune(svm, kernel = kernel,
                    train.x = X, train.y = y, ranges = list(gamma = gamma_start * 2^(-5:5), cost = 2^(-10:8)), class.weights = wgt,
                    tunecontrol = tune.control(sampling = tune.method, cross = tune.cross.k, nboot = tune.boot.n))
  if (verbose) cat("DONE!\n")

  # svm
  if (verbose) cat("SVM modelling...")
  m <- svm(x = X, y = y, kernel = kernel, cost = svm_tuned$best.parameters$cost, gamma = svm_tuned$best.parameters$gamma,
           scale = FALSE, class.weight = wgt, coef0 = ifelse(is.null(svm_tuned$coef.0), 0, svm_tuned$coef.0),
           cross = svm.cross.k,...)
  if (verbose) cat("Done!\n")

  # return
  m$inputX <- x
  m$scaled <- NULL
  m$x.scale <- NULL
  m$y.scale <- NULL
  m$center.scaledX <- centered_X
  m$class.weights <- wgt
  class(m) <- c("svm", "rbiosvm")
  return(m)
}


#' @title rbioFS_svm_roc_auc()
#'
#' @description ROC-AUC analysis and ploting for SVM model
#' @param object A \code{rbiosvm} object.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.label The correspoding label vector to the data. Make sure it is a \code{factor} object.
#' @param center.scale.newdata Logical, wether center and scale the newdata with training data mean and standard deviation. Default is \code{TRUE}.
#' @param rocplot If to generate a ROC plot. Default is \code{TRUE}.
#' @param plot.smooth If to smooth the curves. Uses binormal method to smooth the curves. Default is \code{FALSE}.
#' @param plot.comps Number of comps to plot. Default is \code{1:object$ncomp}
#' @param plot.display.Title If to show the name of the y class. Default is \code{TRUE}.
#' @param plot.titleSize The font size of the plot title. Default is \code{10}.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will be affect error or warning messeages.
#' @return Prints AUC values in the console. And a pdf file for ROC plot. The function also exports a ROC results list to the environment.
#' @details Uses pROC module to calculate ROC.
#'
#' Although optional, the \code{newdata} matrix should use training data's column mean and column standard deviation to enter.scale prior to ROC-AUC analysis.
#' The option \code{center.scaled.newdata = FALSE} is used when the whole (training and test sets) data were center.scaled before SVM training and testing.
#'
#' @import ggplot2
#' @import foreach
#' @importFrom pROC roc
#' @importFrom matrixStats colSds
#' @importFrom GGally ggpairs
#' @importFrom grid grid.newpage grid.draw
#' @importFrom RBioplot rightside_y multi_plot_shared_legend
#' @examples
#' \dontrun{
#' rbioFS_plsda_roc_auc(object = model_binary, rocplot = TRUE, plot.comps = 1:2)
#' }
#' @export
rbioFS_svm_roc_auc <- function(object, newdata, newdata.label,
                               center.scale.newdata = TRUE,
                               rocplot = TRUE,
                               plot.smooth = FALSE,
                               plot.rightsideY = TRUE,
                               plot.SymbolSize = 2, plot.display.Title = TRUE, plot.titleSize = 10,
                               plot.fontType = "sans",
                               plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                               plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                               plot.legendSize = 9,
                               plot.Width = 170, plot.Height = 150,
                               verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c('rbiosvm'))) stop("object needs to be \"rbiosvm\" class.\n")
  if (class(newdata) == "data.frame"){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    newdata <- as.matrix(newdata)
  }
  if (is.vector(newdata)) {
    if (!is.vector(object$inputX)) stop("test data should have the same dimension as the training data.")
  } else {
    if (ncol(newdata) != ncol(object$inputX)) stop("test data should have the same dimension as the training data.")
  }
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data.")
  }

  ## process data
  if (center.scale.newdata){ # using training data mean and sd
    if (verbose) cat(paste0("Data center.scaled using training data column mean and sd, prior to modelling..."))
    centered_newdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
    test <- centered_newdata
    if (verbose) cat("DONE!\n")
  } else {
    centered_newdata <- NULL
    test <- newdata
  }

  ## ROC-AUC calculation
  pred <- predict(object, newdata = test)  # prediction

  outcome <- newdata.label  # origial label
  roc_dfm <- foreach(j = 1:length(levels(outcome)), .combine = "rbind") %do% {
    response <- outcome
    levels(response)[-j] <- "others"
    predictor <- dummy(pred)
    predictor <- as.matrix(predictor[, j], ncol = 1)
    splt <- split(predictor, response)  # split function splist array according to a factor
    controls <- splt$others
    cases <- splt[[levels(outcome)[j]]]
    perf <- tryCatch(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth),
                     error = function(err){
                       cat("Curve not smoothable. Proceed without smooth.\n")
                       pROC::roc(controls = controls, cases = cases, smooth = FALSE)
                     })
    if (length(levels(outcome)) == 2){
      cat(paste0("AUC - ", levels(outcome)[j], ": ", perf$auc, "\n"))
    } else {
      cat(paste0(" AUC - ", levels(outcome)[j], " (vs Others): ", perf$auc, "\n"))
    }

    fpr <- 1 - perf$specificities
    tpr <- perf$sensitivities
    mtx <- cbind(fpr, tpr)
    if (length(levels(outcome)) == 2){
      df <- data.frame(mtx, group = rep(levels(outcome)[j], times = nrow(mtx)), row.names = NULL)
    } else {
      df <- data.frame(mtx, group = rep(paste0(levels(outcome)[j], " (vs Others)"), times = nrow(mtx)), row.names = NULL)
    }
    df <- df[order(df$tpr), ]
    return(df)
  }

  ## plotting
  if (rocplot){
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".svm.roc.pdf...", sep = ""))  # initial message

    plt <- ggplot(data = roc_dfm, aes(x = fpr, y = tpr, group = group, colour = group)) +
      geom_line(aes(linetype = group)) +
      geom_point(aes(shape = group), size = plot.SymbolSize) +
      geom_abline(intercept = 0) +
      ggtitle(ifelse(plot.display.Title, "ROC", NULL)) +
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

    # save
    grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".svm.roc.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    if (verbose) cat("Done!\n")
  }
  ## return
  out <- list(svm.roc = roc_dfm, input.newdata = newdata, input.newdata.label = newdata.label,
              newdata.center.scaled = centered_newdata)
  assign(paste(deparse(substitute(object)), "_svm_roc_list", sep = ""), out, envir = .GlobalEnv)
}
