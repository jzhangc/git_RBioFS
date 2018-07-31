#' @title rbioClass_svm
#'
#' @description Support Vector Machine (SVM) modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param center.scale Logical, wether center and scale the data, i.e. subtracting mean (col.mean) and deviding by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param svm.cross.k Fold of cross validation. Default is \code{10}.
#' @param tune.method Parameter tuning method. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{svm} function from \code{e1071} pacakge.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#' \code{svm.cross.k}: fold number for cross validation
#'
#' \code{tune.method}: grid search method
#'
#' \code{tune.cross.k}: fold number for grid search when using cross validation method, i.e. \code{tune.method = "cross"}.
#'
#' \code{tune.boot.n}: n number for grid search when using bootstrap, i.e. \code{tune.method = "boot"}.
#'
#' @details Model is trained with probability calculation enabled, so that \code{\link{rbioClass_svm_predict}} will be able calculate prediction probabilities.
#'
#' Parameter tuning is for gamma (not applicable when \code{kernel = "linear"}) and cost.
#'
#' The function automatically detects if the data is unbalanced or not, and applies class weight accordingly. For unbalanced data, the
#' class weight the inverse of the draw probability: \code{sample size / class size}.
#'
#' When \code{tune.method = "cross"} and the sample size set for \code{tune.cross.k}, the method is effectively "leave-one-out (LOO)" cross-validation.
#' Same goes for \code{svm.cross.k} argument.
#'
#' The center.scale is \code{(x -  col.mean)/col.sd}. Raw data for SVM should be center.scaled. The training data column mean and column sd are used by the
#' center.scale process for the test data for ROC-AUC analysis. Usually, we center.scale the training data and use the col.mean and col.sd for center.scaling the test data.
#'
#' The option \code{center.scale = FALSE} is for prior center.scaled whole data (training + test sets), i.e center scale the data then split for training and test sets.
#' Such process is only used for pure mathematical analysis of the data set, instead of the generalized use for SVM model evaluation and utility (i.e. classify unknown data).
#'
#'
#' @importFrom e1071 svm tune tune.control
#' @importFrom matrixStats colSds
#' @examples
#' \dontrun{
#' svm_model <- rbioClass_svm(x = training_set[, -1], y = training_set[, 1], kernel = "radial", center = TRUE, scale = FALSE)
#' }
#' @export
rbioClass_svm <- function(x, y, center.scale = TRUE,
                       kernel = "radial", svm.cross.k = 10,
                       tune.method = "cross", tune.cross.k = 10, tune.boot.n = 10, ...,
                       verbose = TRUE){
  ## check arguments
  if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x))) stop("x needs to be a matrix, data.frame or vector.")
  if (class(x) == "data.frame" | is.vector(x)){
    if (verbose) cat("x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }
  if (class(y) != "factor"){
    if (verbose) cat("y is converted to factor. \n")
    y <- factor(y, levels = unique(y))
  }
  if (!kernel %in% c("radial", "linear", "polynomial", "sigmoid")) stop("kernel needs to be exactly one of \"radial\", \"linear\", \"polynomial\", or \"sigmoid\".")

  ## data processing
  if (center.scale){
    if (verbose) cat(paste0("Data centered with the scaling prior to modelling.\n"))
    centered_X <- center_scale(x, scale = TRUE)  # center data with the option of scaling
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
  m <- svm(x = X, y = y, kernel = kernel, probability = TRUE,
           cost = svm_tuned$best.parameters$cost, gamma = svm_tuned$best.parameters$gamma,
           scale = FALSE, class.weight = wgt, coef0 = ifelse(is.null(svm_tuned$coef.0), 0, svm_tuned$coef.0),
           cross = svm.cross.k,...)
  if (verbose) cat("Done!\n")

  # return
  m$inputX <- x
  m$inputY <- y
  m$scaled <- NULL
  m$x.scale <- NULL
  m$y.scale <- NULL
  m$center.scaledX <- centered_X
  m$class.weights <- wgt
  m$svm.cross.k <- svm.cross.k
  m$tune.method <- tune.method
  m$tune.cross.k <- if(tune.method == "cross") tune.cross.k else NULL
  m$tune.boot.n <- if(tune.method == "boot") tune.boot.n else NULL
  class(m) <- c("svm", "rbiosvm")
  return(m)
}


#' @title rbioClass_svm_ncv_fs
#'
#' @description Nested cross-validation assessment for SVM classification. It evaluates the overall performace of SVM modelling given the training data.
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param center.scale Logical, wether center and scale the data, i.e. subtracting mean (col.mean) and deviding by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param cross.k Fold of nested cross validation, i.e. outer loop. Default is \code{10}.
#' @param tune.method Parameter tuning method, i.e. innter loop. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{rbioClass_svm}.
#' @param fs.method Feature selection method. Only \code{"rf"} (i.e. random forest) is supported so far. Default is \code{"rf"}.
#' @param rf.ifs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the initial feature selection step of recursive feature selection. Default is \code{501}.
#' @param rf.sfs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the sequential forward selection step of recursive feature selection. Default is \code{501}.
#' @param fs.count.cutoff A integer for feature vote cutoff. Default is outer loop cross-valiation fold \code{cross.k}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Returns a SVM model object, with classes "svm" and "rbiosvm".
#'
#' Additional items for \code{rbiosvm_nestedcv}:
#'
#' \code{cv.fold}: number of (outer) cross-validation fold
#'
#' \code{randomized.sample.index}: randomized sample order for (outer) cross-validation
#'
#' \code{tot.nested.accuracy.summary}: total (i.e. mean) nested cross-validation accouracy
#'
#' \code{nested.accuracy}: accuracy for each cross-validation iteration
#'
#' \code{fs.method}: feature selection method
#'
#' \code{selected.features}
#'
#' \code{nested.fs.count}: total vote counts for each selected feature
#'
#' \code{tune.method}: (inner) loop method for SVM grid search
#'
#' \code{tune.cross.k}: fold number for (inner) loop if cross-validation is chosen
#'
#' \code{tune.boot.n}: iternation number for (inner) loop if bootstrap is chosen
#'
#' @details TBA
#'
#' For now, RBioFS implementation of two-step random forest feature selection is used to select features based on nested cross-validation.
#' Resulted features from each nested cross-validation round are voted. Features with votes equal or greater than the cutoff are reported as selected features.
#' It is also a good idea to set \code{fs.count.cutoff} as \code{cross.k - 1}. Notably, \code{fs.count.cutoff = 1} is "no threshold",
#' meaning maximum number of features selected from the nested cross-validation are reported.
#'
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' svm_nestedcv <- rbioClass_svm_ncv_fs(x = mydata[, -c(1:2)],
#'                                      y = factor(mydata$Conditions, levels = unique(mydata$Conditions)),
#'                                      center.scale = TRUE,
#'                                      cross.k = 5,
#'                                      tune.method = "cross",
#'                                      tune.cross.k = 5,
#'                                      fs.method = "rf",
#'                                      rf.ifs.ntree = 1001, rf.sfs.ntree = 1001,
#'                                      parallelComputing = TRUE, clusterType = "FORK")
#'
#' }
#' @export
rbioClass_svm_ncv_fs <- function(x, y, center.scale = TRUE,
                                 kernel = "radial",
                                 cross.k = 10,
                                 tune.method = "cross",
                                 tune.cross.k = 10, tune.boot.n = 10, ...,
                                 fs.method = "rf", rf.ifs.ntree = 1001, rf.sfs.ntree = 1001,
                                 fs.count.cutoff = cross.k,
                                 parallelComputing = TRUE, clusterType = "PSOCK",
                                 verbose = TRUE){
  ## check arguments
  if (!fs.method %in% c("rf")) stop("So far, fs.method has to be \"rf\". More methods will be implemented")
  if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
  if (cross.k > nrow(x)) stop("Cross-validation fold setting cross.k exceeded limit. Hint: max at total sample number.\n")
  if (class(x) == "data.frame"){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }
  if (class(y) != "factor"){
    if (verbose) cat("y is converted to factor. \n")
    y <- factor(y, levels = unique(y))
  }
  if (!kernel %in% c("radial", "linear", "polynomial", "sigmoid")) stop("kernel needs to be exactly one of \"radial\", \"linear\", \"polynomial\", or \"sigmoid\".")
  if (fs.count.cutoff %% 1 != 0 | fs.count.cutoff < 1 | fs.count.cutoff > cross.k) stop("fs.count.cutoff should be an integer between 1 and cross.k.")

  ## nested cv
  # processing training data
  dfm <- data.frame(y, x, check.names = FALSE)
  random_sample_idx <- sample(nrow(dfm))
  dfm_randomized <- dfm[random_sample_idx, ]  # randomize samples
  fold <- cut(seq(1:nrow(dfm_randomized)), breaks = cross.k, labels = FALSE)  # create a factor object with fold indices

  # nested cv
  if (verbose) cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
  if (verbose) cat(paste0("Data center.scale: ", ifelse(center.scale, " ON\n", " OFF\n")))
  if (verbose) cat("Nested cross-validation with feature selection (speed depending on hardware configuration)...")
  nested.cv.list <- vector(mode = "list", length = cross.k)
  if (!parallelComputing){  # single core
    nested.cv.list[] <- foreach(i = 1:cross.k) %do% {
      training <- dfm_randomized[which(fold != i, arr.ind = TRUE), ]
      # fs
      if (center.scale){
        fs_training <- center_scale(training[, -1], scale = TRUE)$centerX  # without y
      } else {
        fs_training <- training[, -1]
      }
      rbioFS_rf_initialFS(objTitle = "svm_nested", x = fs_training, targetVar = training$y, nTimes = 50,
                          nTree = rf.ifs.ntree, multicore = FALSE, plot = FALSE)
      # fs <- svm_nested_initial_FS$feature_initial_FS
      rbioFS_rf_SFS(objTitle = "svm_nested", x = svm_nested_initial_FS$training_initial_FS, targetVar = training$y, nTimes = 50,
                    nTree = rf.sfs.ntree, mTry = "recur_default", multicore = FALSE, plot = FALSE)
      if (length(svm_nested_SFS$selected_features) > 1){
        fs <- svm_nested_SFS$selected_features
      } else {
        fs <- svm_nested_initial_FS$feature_initial_FS
      }

      # cv svm
      m <- rbioClass_svm(x = training[, -1][, fs], y = training$y, center.scale = center.scale,
                         svm.cross.k = 0, tune.method = tune.method,
                         tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE)

      # processing test data
      test <- dfm_randomized[which(fold == i, arr.ind = TRUE), ][, c("y", fs)]  # preseve y and selected fetures
      if (center.scale){ # using training data mean and sd
        centered_newdata <- t((t(test[, -1]) - m$center.scaledX$meanX) / m$center.scaledX$columnSD)
        test[, -1] <- centered_newdata
      } else {
        centered_newdata <- NULL
      }
      pred <- predict(m, newdata = test[, -1])
      accu <- sum(diag(table(pred, test$y))) / length(test$y)  # accuracy = total TP / total (TP: true positive)

      # foreach output
      tmp_out <- list(selected.features = fs, nested.cv.accuracy = accu)
      tmp_out
    }
  } else {  # multiple cores
    # set clusters
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores, clusterType = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # computing
    nested.cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS")) %dopar% {
      training <- dfm_randomized[which(fold != i, arr.ind = TRUE), ]
      # fs
      if (center.scale){
        fs_training <- center_scale(training[, -1], scale = TRUE)$centerX  # without y
      } else {
        fs_training <- training[, -1]
      }
      rbioFS_rf_initialFS(objTitle = "svm_nested", x = fs_training, targetVar = training$y, nTimes = 50,
                          nTree = rf.ifs.ntree, multicore = FALSE, plot = FALSE)
      # fs <- svm_nested_initial_FS$feature_initial_FS
      rbioFS_rf_SFS(objTitle = "svm_nested", x = svm_nested_initial_FS$training_initial_FS, targetVar = training$y, nTimes = 50,
                    nTree = rf.sfs.ntree, mTry = "recur_default", multicore = FALSE, plot = FALSE)
      if (length(svm_nested_SFS$selected_features) > 1){
        fs <- svm_nested_SFS$selected_features
      } else {
        fs <- svm_nested_initial_FS$feature_initial_FS
      }

      # cv svm
      m <- rbioClass_svm(x = training[, -1][, fs], y = training$y, center.scale = center.scale,
                         svm.cross.k = 0, tune.method = tune.method,
                         tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE)

      # processing test data
      test <- dfm_randomized[which(fold == i, arr.ind = TRUE), ][, c("y", fs)]  # preseve y and selected fetures
      if (center.scale){ # using training data mean and sd
        centered_newdata <- t((t(test[, -1]) - m$center.scaledX$meanX) / m$center.scaledX$columnSD)
        test[, -1] <- centered_newdata
      } else {
        centered_newdata <- NULL
      }
      pred <- predict(m, newdata = test[, -1])
      accu <- sum(diag(table(pred, test$y))) / length(test$y)  # accuracy = total TP / total (TP: true positive)

      # foreach output
      tmp_out <- list(selected.features = fs, nested.cv.accuracy = accu)
      tmp_out
    }
  }
  names(nested.cv.list) <- paste0("cv_fold_", c(1:cross.k))

  nested.accu <- foreach(i = 1:cross.k, .combine = "c") %do% {
    nested.cv.list[[i]]$nested.cv.accuracy
  }
  nested.fs <- foreach(i = 1:cross.k, .combine = "c") %do% {
    nested.cv.list[[i]]$selected.features
  }
  fs.count <- sort(table(nested.fs), decreasing = TRUE)
  if (verbose) cat("Done!\n")

  ## output
  tot.nested.acc.summary <- c(mean(nested.accu), sd(nested.accu), sd(nested.accu)/sqrt(cross.k))
  names(tot.nested.acc.summary) <- c("tot.nested.accuracy", "sd", "sem")
  selected.features <- names(fs.count[which(fs.count >= fs.count.cutoff)])

  # display
  if (verbose) cat("\n")
  if (verbose) cat("Nested cross-validation accuracy summary: \n")
  if (verbose) print(tot.nested.acc.summary)
  if (verbose) cat("\n")
  if (verbose) cat("Nested cross-validation selected features: \n")
  if (verbose) print(selected.features)

  # export to environment
  out <- list(cv.fold = fold,
              randomized.sample.index = random_sample_idx,
              tot.nested.accuracy.summary = tot.nested.acc.summary,
              nested.accuracy = nested.accu,
              fs.method = fs.method,
              selected.features = selected.features,
              nested.fs.count = fs.count,
              tune.method = tune.method,
              tune.cross.k = if(tune.method == "cross") tune.cross.k else NULL,
              tune.boot.n = if(tune.method == "boot") tune.boot.n else NULL)

  class(out) <- "rbiosvm_nestedcv"
  return(out)
}


#' @title rbioClass_svm_roc_auc()
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
#' @param plot.fontType The type of font in the figure. Default is "sans". For all options please refer to R font table, which is avaiable on the website: \url{http://kenstoreylab.com/?page_id=2448}.
#' @param plot.SymbolSize Symbol size. Default is \code{2}.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#' rbioClass_plsda_roc_auc(object = model_binary, rocplot = TRUE, plot.comps = 1:2)
#' }
#' @export
rbioClass_svm_roc_auc <- function(object, newdata, newdata.label,
                                  center.scale.newdata = TRUE,
                                  rocplot = TRUE,
                                  plot.smooth = FALSE,
                                  plot.SymbolSize = 2, plot.display.Title = TRUE, plot.titleSize = 10,
                                  plot.fontType = "sans",
                                  plot.xLabel = "1 - specificity", plot.xLabelSize = 10, plot.xTickLblSize = 10,
                                  plot.yLabel = "sensitivity", plot.yLabelSize = 10, plot.yTickLblSize = 10,
                                  plot.legendSize = 9, plot.rightsideY = TRUE,
                                  plot.Width = 170, plot.Height = 150,
                                  verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c('rbiosvm'))) stop("object needs to be \"rbiosvm\" class.")
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x)))stop("x needs to be a matrix, data.frame or vector.")
  if (class(newdata) == "data.frame" | is.null(dim(newdata))){
    if (verbose) cat("newdata converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (ncol(newdata) != ncol(object$inputX)) stop("test data should have the same dimension as the training data.")
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


#' @title rbioClass_svm_perm()
#'
#' @description Permutation test for SVM models.
#' @param object A \code{rbiosvm} object. Make sure the object is generated with a \code{tot.accuracy} section.
#' @param perm.method Permutation method. Options are \code{"by_y"} and \code{"by_feature_per_y"}. Default is \code{"by_y"}. See details below.
#' @param nperm Number of permutations to run. Default is \code{999}.
#' @param perm.plot Wether to produce a plot or not. Default is \code{TRUE}.
#' @param ... Additional argument for \code{\link{rbioUtil_perm_plot}}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return The function returns \code{CSV} files for all intermediate permuatation accuracy values as well as the p-value resutls.
#'
#'
#' The final results are also exported to the environment as a \code{rbiosvm_perm} object with the following items:
#'
#' \code{original.accuracy} The original accuracy for the SVM model.
#'
#' \code{p.value} P value for the permutation test.
#'
#' \code{perm.method} Permutation method.
#'
#' \code{perm.stats} The stats metric used for the permutation test.
#'
#' \code{nperm} The number of permutation runs.
#'
#' \code{perm.results} The intermediate permutation results, i.e. stats for each permutation test run in a data.frame. \code{nperm = 0} is the original stats.
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
#' For \code{perm.method = "by_y"}, labels (i.e. y) are permuated. A non-signifianct model (permutation p value > alpha, i.e. 0.05) in this case means the data is independent from the groups.
#'
#' For \code{perm.method = "by_feature_per_by"}, X is first subset by label (i.e.y) before permutating data for each feature. Since the permutation is done for the features WITHIN the group,
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
                               perm.method = "by_y", nperm = 999,
                               perm.plot = TRUE, ...,
                               parallelComputing = TRUE, clusterType = "PSOCK",
                               verbose = TRUE){
  ## check arguments
  if (!any(class(object) %in% c("rbiosvm"))) stop("object has to be a \"rbiosvm\" class.\n")
  if (!"tot.accuracy" %in% names(object) || is.null(object$tot.accuracy)) stop("SVM model has to include tot.accuracy value from Cross-Validation.\n")
  if (!perm.method %in% c("by_y", "by_feature_per_y")) stop("perm.method needs to be either \"by_y\" or \"by_feature_per_y\". \n")
  if (length(nperm) != 1) stop("nperm can only contain one integer. \n")
  if (nperm %% 1 != 0) stop("nperm can only be integer. \n")
  if (nperm < 1) stop("nperm can only take interger equal to or greater than 1. \n")

  ## test
  ## calcuate RMSEP and construct original RMSEP data frame
  orig_accu <- object$tot.accuracy

  ## permutation test
  if (verbose) cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
  if (verbose) cat(paste0("Running permutation test using ", perm.method, " method ","with ", nperm, " permutations (speed depending on hardware configurations)..."))
  # consolidated permutation model RMSEP data frame
  if (!parallelComputing){
    if (perm.method == "by_y"){  # permutate label
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %do% {
        perm_y <- object$inputY[sample(1:length(object$inputY))]  # sample label permutation
        perm_model <- rbioClass_svm(x = object$center.scaledX$centerX, y = factor(perm_y, levels = unique(perm_y)),
                                    center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                    tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                    verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
        perm_accu <- perm_model$tot.accuracy  # permutation model accuracy

        perm_accu_dfm <- data.frame(nperm = i, stats = perm_accu, row.names = NULL)
        perm_accu_dfm
      }

    } else {  # subset by label first, then permutate data per feature
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind") %do% {
        perm_x <- foreach(m = unique(levels(object$inputY)), .combine = "rbind") %do% {
          sub_dat <- object$center.scaledX$centerX[which(object$inputY == m), ]  # subsetting the centre-scaled X by label (Y)
          perm_dat <- foreach(n = 1:ncol(sub_dat), .combine = "cbind") %do% {
            col <- sub_dat[sample(1:nrow(sub_dat)), n, drop = FALSE]  # permutation for each feature
            col
          }
          perm_dat
        }

        perm_model <- rbioClass_svm(x = perm_x, y = object$inputY,
                                    center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                    tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                    verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
        perm_accu <- perm_model$tot.accuracy   # permutation model accuracy
        perm_accu_dfm <- data.frame(nperm = i, stats = perm_accu, row.names = NULL)
        perm_accu_dfm
      }
    }
  } else {  # parallel computing
    # set up cpu cluster
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores, clusterType = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # permutation test
    if (perm.method == "by_y"){  # permutate label
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind", .packages = c("foreach", "RBioFS")) %dopar% {
        perm_y <- object$inputY[sample(1:length(object$inputY))]  # sample label permutation
        perm_model <- rbioClass_svm(x = object$center.scaledX$centerX, y = factor(perm_y, levels = unique(perm_y)),
                                    center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                    tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                    verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
        perm_accu <- perm_model$tot.accuracy  # permutation model accuracy

        perm_accu_dfm <- data.frame(nperm = i, stats = perm_accu, row.names = NULL)
        perm_accu_dfm
      }
    } else {  # subset by label first, then permutate data per feature
      perm_dfm <- foreach(i = 1:nperm, .combine = "rbind", .packages = c("foreach", "RBioFS")) %dopar% {
        perm_x <- foreach(m = unique(levels(object$inputY)), .combine = "rbind") %do% {
          sub_dat <- object$center.scaledX$centerX[which(object$inputY == m), ]  # subsetting the centre-scaled X by label (Y)
          perm_dat <- foreach(n = 1:ncol(sub_dat), .combine = "cbind") %do% {
            col <- sub_dat[sample(1:nrow(sub_dat)), n, drop = FALSE]  # permutation for each feature
            col
          }
          perm_dat
        }

        perm_model <- rbioClass_svm(x = perm_x, y = object$inputY,
                                    center.scale = FALSE, svm.cross.k = object$svm.cross.k, tune.method = object$tune.method,
                                    tune.cross.k = object$tune.cross.k, tune.boot.n = object$tune.boot.n,
                                    verbose = FALSE)  # permutated data modelling. NOTE: the x data is already scaled and centred
        perm_accu <- perm_model$tot.accuracy   # permutation model accuracy
        perm_accu_dfm <- data.frame(nperm = i, stats = perm_accu, row.names = NULL)
        perm_accu_dfm
      }
    }
  }

  perm_accu_val <- perm_dfm$tot.accuracy
  p.val <- (length(which(perm_accu_val >= orig_accu)) + 1) / (nperm + 1)
  if (verbose) cat("Done!\n")

  ## output
  perm_stats <- rbind(data.frame(nperm = 0, stats = orig_accu, row.names = NULL), perm_dfm)
  out <- list(perm.method = perm.method, nperm = nperm,  perm.stats = "tot.accuracy",
                       original.stats = orig_accu, p.value = p.val, perm.results = perm_stats)
  class(out) <- "rbiosvm_perm"
  assign(paste(deparse(substitute(object)), "_perm", sep = ""), out, envir = .GlobalEnv)

  # print results
  if (verbose) cat("\n")
  if (verbose) cat("Permutation test restuls (p value):")
  if (verbose) cat(p.val)
  if (verbose) cat("\n")

  # export to the directory
  if (verbose) cat(paste0("Permutation stats results stored in csv file: ", paste(deparse(substitute(object)), "_perm.csv", sep = ""), ". \n"))
  write.csv(file = paste0(deparse(substitute(object)), ".perm.csv"), perm_stats, row.names = FALSE)

  ## plot
  if (perm.plot){
    svm_permutation_test <- out
    rbioUtil_perm_plot(perm_res = svm_permutation_test, ...)
  }
}


#' @title rbioClass_svm_predict
#'
#' @description Prediction function for SVM analysis. The function calculates the predicted value for unknown sample data using the input SVM model.
#' @param object A \code{rbiosvm} object.
#' @param newdata Input data to be classified. Make sure it is a \code{matrix} class and has the same variables as the model, i.e. same number of columns as the training data.
#' @param sampleLabel.vector A character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param center.scale.newdata If to center the newdata. When \code{TRUE}, it will also apply the same scaling option as the \code{object}. Default is \code{TRUE}.
#' @param prob.method Method to calculate classification probability. Options are \code{"logistic"}, \code{"softmax"} and \code{"Bayes"}. See details for more information. Default is \code{"logistic"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return  A \code{prediction} obejct. The items of the object are:
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
#' @details Although optional, the \code{newdata} matrix should be centered prior to testing, with the same scaling setting as the input \code{rbiosvm} object. The option \code{center.scale.newdata = FALSE} is
#' for the already centered the data matrix. This center.scale process should use training data's column mean and column standard deviation.
#'
#' The default posterior probability calculation method \code{"logistic"} is the \code{e1071} pacakge's implementation of logistic regression model.
#' See \code{\link{rbioClass_plsda_predict}} for description for "Bayes" and "softmax" method.
#'
#' If \code{sampleLabel.vector = NULL} or missing, the function uses row numbers as label.
#'
#' @import ggplot2
#' @import pls
#' @importFrom klaR NaiveBayes
#' @examples
#' \dontrun{
#'
#'
#' }
#' @export
rbioClass_svm_predcit <- function(object, newdata, sampleLabel.vector = NULL,
                                  center.scale.newdata = TRUE, prob.method = "logistic",
                                  verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c("rbiosvm"))) stop("object needs to be a \"rbiosvm\" class.")
  if (!class(newdata) %in% c("matrix", "data.frame")) stop("newdata has to be either a matrix or data.frame object.")
  if (ncol(newdata) != ncol(object$inputX)) stop("newdata needs to have the same number of variables, i.e. columns, as the object.")
  if (!prob.method %in% c("logistic",  "Bayes", "softmax")) stop("Probability method should be either \"softmax\" or \"Bayes\".")
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data while center.scale.newdata = TRUE.")
  }
  if (!missing(sampleLabel.vector) & !is.null(sampleLabel.vector)){
    if (length(sampleLabel.vector) != nrow(newdata)){
      cat("")
      sampleLabel.vector <- NULL
    }
  }

  ## center data with the option of scaling
  if (class(newdata) == "data.frame"){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (center.scale.newdata){
    if (verbose) cat("Data center.scaled using training data column mean and sd, prior to modelling.")
    centerdata <- t((t(newdata) - object$center.scaledX$meanX) / object$center.scaledX$columnSD)
    test <- centerdata
  } else {
    centerdata <- NULL
    test <- newdata
  }

  ## predict
  rownames(test) <- 1:nrow(test)
  # pred <- predict(object = object, newdata = test, type = "response", probability = TRUE)
  pred <- predict(object = object, newdata = test, probability = TRUE)
  pred_mtx <- attr(pred, "probabilities")  # probability matrix

  if (is.null(dim(pred_mtx))) {  # if only one sample
    pred_mtx <- t(as.matrix(pred_mtx))
  }
  if (!missing(sampleLabel.vector) & !is.null(sampleLabel.vector)){
    rownames(pred_mtx) <- sampleLabel.vector
  }

  ## posterior
  if (prob.method == "logistic"){
    prob <- pred_mtx  # probability matrix
  } else if (prob.method == "Bayes"){
    if (!is.null(object$center.scaledX)){
      training_mtx <- as.matrix(object$center.scaledX$centerX)
    } else {
      training_mtx <- as.matrix(object$inputX)
    }
    trainingpred <- predict(object = object, newdata = training_mtx, type = "response", probability = TRUE)
    trainingpred <- attributes(trainingpred)$probabilities
    bayes.prob <- klaR::NaiveBayes(x = trainingpred, grouping = object$inputY, usekernel = TRUE)
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

  ## output
  predicted.value <- pred
  attributes(predicted.value) <- NULL
  out <- list(classifier.class = class(object), predicted.value = predicted.value, prob.method = prob.method, probability.summary = prob_dfm,
              raw.newdata = newdata, center.scale = center.scale.newdata, center.scaled.newdata = centerdata)
  class(out) <- "prediction"
  assign(paste(deparse(substitute(object)), "_svm_predict", sep = ""), out, envir = .GlobalEnv)
}
