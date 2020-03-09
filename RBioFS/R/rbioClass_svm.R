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
#' \code{model.type}: the SVM model type, "classification" or "regression".
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
#' The function also supports regression study, in which case, the performance metric is \code{RMSE}.
#'
#' @importFrom e1071 svm tune tune.control
#' @importFrom matrixStats colSds
#' @examples
#' \dontrun{
#' svm_model <- rbioClass_svm(x = training_set[, -1], y = training_set[, 1], kernel = "radial", center = TRUE, scale = FALSE)
#' }
#' @export
rbioClass_svm <- function(x, y, center.scale = TRUE,
                          kernel = c("radial", "linear", "polynomial", "sigmoid"), svm.cross.k = 10,
                          tune.method = c("cross", "boot", "fix"), tune.cross.k = 10, tune.boot.n = 10, ...,
                          verbose = TRUE){
  ## check arguments
  if (is.factor(y)){
    if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
    model_type <- "classification"
    y <- factor(y, levels = unique(y))
  } else {
    model_type <- "regression"
  }
  if (!class(x) %in% c("data.frame", "matrix") & !is.null(dim(x))) stop("x needs to be a matrix, data.frame or vector.")
  if (class(x) == "data.frame" | is.vector(x)){
    if (verbose) cat("x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }
  # if (class(y) != "factor"){
  #   if (verbose) cat("y is converted to factor. \n")
  #   y <- factor(y, levels = unique(y))
  # }
  # if (!kernel %in% c("radial", "linear", "polynomial", "sigmoid")) stop("kernel needs to be exactly one of \"radial\", \"linear\", \"polynomial\", or \"sigmoid\".")
  kernel <- match.arg(tolower(kernel), c("radial", "linear", "polynomial", "sigmoid"))
  tune.method <- match.arg(tolower(tune.method), c("cross", "boot", "fix"))


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

  ## weight evaluation according to the model type
  if (model_type == "classification"){
    if (length(unique(table(y))) != 1){  # test if the sample is balanced
      if(any(table(y) == 0)) {
        warning("Not all groups present in the unbalanced training data, discard missing group and set weight for other groups as 1")
        wgt <- rep(1, length(unique(y)))
        names(wgt) <- unique(y)
      } else {
        wgt <- length(y) / table(y)
      }
    } else {  # balanced sample, each class weight is 1
      wgt <- unique(table(y)) / table(y)
    }
  } else{
    wgt <- NULL
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
           cross = svm.cross.k, ...)
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
  m$model.type <- model_type
  class(m) <- c("rbiosvm", "svm")
  return(m)
}


#' @export
print.rbiosvm <- function(x, ...){
  cat("Parameters:\n")
  cat("   SVM-Type: ", c("C-classification", "nu-classification",
                         "one-classification", "eps-regression", "nu-regression")[x$type +
                                                                                    1], "\n")
  cat(" SVM-Kernel: ", c("linear", "polynomial", "radial",
                         "sigmoid")[x$kernel + 1], "\n")
  if (x$type == 0 || x$type == 3 || x$type == 4)
    cat("       cost: ", x$cost, "\n")
  if (x$kernel == 1)
    cat("     degree: ", x$degree, "\n")
  cat("      gamma: ", x$gamma, "\n")
  if (x$kernel == 1 || x$kernel == 3)
    cat("     coef.0: ", x$coef0, "\n")
  if (x$type == 1 || x$type == 2 || x$type == 4)
    cat("         nu: ", x$nu, "\n")
  if (x$type == 3) {
    cat("    epsilon: ", x$epsilon, "\n")
    if (x$compprob)
      cat("Sigma: ", x$sigma, "\n")
  }
  cat("\nNumber of Support Vectors: ", x$tot.nSV)
  cat("\n")
}


#' @title rbioClass_svm_ncv_fs
#'
#' @description Nested cross-validation assessment for SVM classification. It evaluates the overall performace of SVM modelling given the training data.
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param univariate.fs If to use limma-based univariate reduction. Default is \code{FALSE}.
#' @param uni.log2trans Only set if \code{univariate.fs = TRUE}, if to log2 transform data before univariate reduction.
#' @param uni.contrast Only set if \code{univariate.fs = TRUE} and for a classificaiton study, the contrast for the univariate analysis. Default is \code{NULL}.
#' @param uni.alpha Only set if \code{univariate.fs = TRUE}, the p value alpha. Default is \code{0.05}.
#' @param uni.fdr Only set if \code{univariate.fs = TRUE}, if to use FDR for the p value. Default is \code{FALSE}.
#' @param center.scale Logical, wether center and scale the data, i.e. subtracting mean (col.mean) and deviding by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param cross.k Fold of nested cross validation, i.e. outer loop. Default is \code{10}.
#' @param cross.best.model.method The method to select the best cv models for feature selection. Options are \code{"median"} and \code{"none"}. Default is \code{"median"}.
#' @param tune.method Parameter tuning method, i.e. innter loop. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{rbioClass_svm}.
#' @param fs.method Feature selection method. Only \code{"rf"} (i.e. random forest) is supported so far. Default is \code{"rf"}.
#' @param rf.ifs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the initial feature selection step of recursive feature selection. Default is \code{501}.
#' @param rf.sfs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the sequential forward selection step of recursive feature selection. Default is \code{501}.
#' @param fs.count.cutoff A integer for feature vote cutoff. Default is outer loop cross-valiation fold \code{cross.k}.
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return Returns a SVM model object, with classes "svm" and "rbiosvm".
#'
#' Additional items for \code{rbiosvm_nestedcv}:
#'
#' \code{univariate.fs}
#'
#' \code{cv.fold}: number of (outer) cross-validation fold
#'
#' \code{randomized.sample.index}: randomized sample order for (outer) cross-validation
#'
#' \code{model.type}: the SVM model type, "classification" or "regression".
#'
#' \code{tot.nested.accuracy.summary}: total (i.e. mean) nested cross-validation accuracy, if \code{model_type = "classification"}
#'
#' \code{tot.nested.RMSE.summary}: total (i.e. mean) nested cross-validation RMSE, if \code{model_type = "regression"}
#'
#' \code{tot.nested.rsq.summary}: total (i.e. mean) nested cross-validation R2, if \code{model_type = "regression"}
#'
#' \code{nested.accuracy}: accuracy for each cross-validation iteration, if \code{model_type = "classification"}
#'
#' \code{nested.RMSE}: RMSE for each cross-validation iteration, if \code{model_type = "regression"}
#'
#' \code{nested.rsq}: R2 for each cross-validation iteration, if \code{model_type = "regression"}
#'
#' \code{nested.fs.count}: fs count for all the cv models
#'
#' \code{nested.cv.models}: all the cv models with test data. NOTE: the test data is center-scaled with training data if \code{center.scale = TRUE}.
#'
#' \code{best.nested.method}
#'
#' \code{best.nested.index}: index of the best models in the total model list
#'
#' \code{best.nested.acc.summary}
#'
#' \code{best.nested.rmse.summary}
#'
#' \code{best.nested.rsq.summary}
#'
#' \code{best.nested.accu}
#'
#' \code{best.nested.rmse}
#'
#' \code{best.nested.rsq}
#'
#' \code{best.nested.fs.count}: fs count for only the best models
#'
#' \code{fs.method}: feature selection method
#'
#' \code{fs.count.threshold}: count threshold for generating consensus feature list
#'
#' \code{selected.features}
#'
#' \code{tune.method}: (inner) loop method for SVM grid search
#'
#' \code{tune.cross.k}: fold number for (inner) loop if cross-validation is chosen
#'
#' \code{tune.boot.n}: iternation number for (inner) loop if bootstrap is chosen
#'
#' \code{run.time}: total run time
#'
#' @details
#'
#' For now, RBioFS implementation of two-step random forest feature selection is used to select features based on nested cross-validation.
#' Resulted features from each nested cross-validation round are voted. Features with votes equal or greater than the cutoff are reported as selected features.
#'
#' It is also a good idea to set \code{fs.count.cutoff} as \code{cross.k - 1}. Notably, \code{fs.count.cutoff = 1} is "no threshold",
#' meaning maximum number of features selected from the nested cross-validation are reported.
#'
#' When \code{cross.best.model.method = "median"}, the function only use models with accuracy/RMSE equal or better than the median valaue
#' for feature count threholding. When there is no change in perforamce across cv models, the function behaves same as \code{cross.best.model.method = "none"}
#'
#'
#' The function also supports regression study, in which case, the performance metric is \code{RMSE}.
#'
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma lmFit topTable makeContrasts eBayes
#' @importFrom splines ns
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
rbioClass_svm_ncv_fs <- function(x, y,
                                 univariate.fs = FALSE, uni.log2trans = TRUE, uni.contrast = NULL,
                                 uni.alpha = 0.05, uni.fdr = FALSE,
                                 center.scale = TRUE,
                                 kernel = c("radial", "linear", "polynomial", "sigmoid"),
                                 cross.k = 10, cross.best.model.method = c("none", "median"),
                                 tune.method = c("cross", "boot", "fix"),
                                 tune.cross.k = 10, tune.boot.n = 10, ...,
                                 fs.method = "rf", rf.ifs.ntree = 1001, rf.sfs.ntree = 1001,
                                 fs.count.cutoff = cross.k,
                                 parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                                 verbose = TRUE){
  ## initiate the run time and set seed
  start_time <- Sys.time()

  ## check arguments
  if (!fs.method %in% c("rf")) stop("So far, fs.method has to be \"rf\". More methods will be implemented")
  if (is.factor(y)) {
    if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
    model_type <- "classification"
  } else {
    model_type <- "regression"
  }
  # if (class(y) != "factor"){
  #   if (verbose) cat("y is converted to factor. \n")
  #   y <- factor(y, levels = unique(y))
  # }
  if (model_type == "classification" && univariate.fs) {
    if (is.null(uni.contrast)) {
      stop("When univariate.fs = TRUE, uni.contrast needs to be set for classification study.")
    } else if (!is.character(uni.contrast)) {
      stop("uni.contrast needs to be a character string.")
    }
  }

  if (cross.k > nrow(x)) stop("Cross-validation fold setting cross.k exceeded limit. Hint: max at total sample number.\n")
  if (class(x) == "data.frame"){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }

  if (fs.count.cutoff %% 1 != 0 | fs.count.cutoff < 1 | fs.count.cutoff > cross.k) stop("fs.count.cutoff should be an integer between 1 and cross.k.")
  # if (!kernel %in% c("radial", "linear", "polynomial", "sigmoid")) stop("kernel needs to be exactly one of \"radial\", \"linear\", \"polynomial\", or \"sigmoid\".")
  kernel <- match.arg(tolower(kernel), c("radial", "linear", "polynomial", "sigmoid"))
  cross.best.model.method <- match.arg(tolower(cross.best.model.method), c("none", "median"))
  tune.method <- match.arg(tolower(tune.method), c("cross", "boot", "fix"))
  if (parallelComputing){
    clusterType <- match.arg(clusterType, c("PSOCK", "FORK"))
  }

  ## nested cv
  # processing training data
  dfm <- data.frame(y, x, check.names = FALSE)
  random_sample_idx <- sample(nrow(dfm))
  dfm_randomized <- dfm[random_sample_idx, ]  # randomize samples
  fold <- cut(seq(1:nrow(dfm_randomized)), breaks = cross.k, labels = FALSE)  # create a factor object with fold indices

  # nested CV function
  nestedcv_func <- function(i) {
    if (verbose) cat(paste0("Nested CV iteration: ", i, "|", cross.k, "..."))
    cv_training <- dfm_randomized[which(fold != i, arr.ind = TRUE), ]
    cv_training_x <- cv_training[, -1]
    cv_training_y <- cv_training[, 1]

    # optional uni
    # outcome should be a list of features to subset the cv_training first
    if (univariate.fs){
      if (model_type == "classification") {
        cv_design <- model.matrix(~ 0 + cv_training_y)
        colnames(cv_design) <- levels(cv_training_y)
      } else {
        nsy <- splines::ns(cv_training_y)
        cv_design <- model.matrix(~ nsy)
      }

      if (uni.log2trans) {
        E <- apply(t(cv_training_x), c(1, 2), FUN = function(x)log2(x + 2))
      } else {
        E <- t(cv_training_x)
      }

      cv_pair <- data.frame(ProbeName = seq(ncol(cv_training) - 1), pair = colnames(cv_training)[-1])
      cv_sample <- paste0(colnames(E), "_", cv_training_y)
      cv_idx <- data.frame(sampleid = colnames(E), group = cv_training_y, sample = cv_sample)
      cv_rawlist <- list(E = E, genes = cv_pair, targets = cv_idx)
      cv_normdata <- uni_PreProc(rawlist = cv_rawlist, offset = 2, normMethod = "quantile", bgMethod = "none")

      cv_fit <- lmFit(cv_normdata$E, design = cv_design, weights = cv_normdata$ArrayWeight)
      if (model_type == "classification"){
        contra_string <- unlist(strsplit(uni.contrast, split = ","))
        contra_string <- gsub(" ", "", contra_string, fixed = TRUE)  # remove all the white space
        # NOTE: below: use do.call to unpack arguments for makeContrasts()
        cv_contra <- do.call(makeContrasts, c(as.list(contra_string), levels = as.list(parse(text = "cv_design"))))
        cv_fit <- contrasts.fit(cv_fit, contrasts = cv_contra)
      }
      cv_fit <- eBayes(cv_fit)
      cv_fit_dfm <- suppressMessages(topTable(cv_fit, number = Inf))
      cv_fit_dfm$feature <- rownames(cv_fit_dfm)

      if (uni.fdr){ # FDR
        sig_idx <- which(cv_fit_dfm$adj.P.Val <= uni.alpha)
        if (length(sig_idx) < 1){
          pcutoff <- uni.alpha
        } else {
          pcutoff <- max(cv_fit_dfm$P.Value[sig_idx])
        }
      } else { # NON-FDR
        pcutoff <- uni.alpha
      }
      # output
      uni_sig_fs <- as.character(cv_fit_dfm[cv_fit_dfm$P.Value < pcutoff, "feature"])

      # update the cv training data
      cv_training <- cv_training[, c("y", uni_sig_fs)]
      cv_training_x <- cv_training[, -1]
      cv_training_y <- cv_training[, 1]
    } else {
      uni_sig_fs <- NULL
    }

    # fs
    if (center.scale){
      fs_training_x <- center_scale(cv_training_x, scale = TRUE)$centerX  # without y
    } else {
      fs_training_x <- cv_training_x
    }
    rbioFS_rf_initialFS(objTitle = paste0("svm_nested_iter_", i), x = fs_training_x, y = cv_training_y, nTimes = 50,
                        nTree = rf.ifs.ntree, parallelComputing = parallelComputing, clusterType = clusterType, plot = FALSE)
    # fs <- svm_nested_initial_FS$feature_initial_FS
    rbioFS_rf_SFS(objTitle = paste0("svm_nested_iter_", i),
                  x = eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$training_initial_FS, y = cv_training_y, nTimes = 50,
                  nTree = rf.sfs.ntree, parallelComputing = parallelComputing, clusterType = clusterType, plot = FALSE)

    if (length(eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features) > 1){
      fs <- eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features
    } else {
      fs <- eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$feature_initial_FS
    }

    # cv svm
    cv_m <- rbioClass_svm(x = fs_training_x[, fs], y = cv_training_y, center.scale = center.scale,
                          svm.cross.k = 0, tune.method = tune.method,
                          tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE, ...)

    # processing test data
    fs_test <- dfm_randomized[which(fold == i, arr.ind = TRUE), ][, c("y", fs)]  # preseve y and selected fetures
    if (center.scale){ # using training data mean and sd
      centered_newdata <- t((t(fs_test[, -1]) - cv_m$center.scaledX$meanX) / cv_m$center.scaledX$columnSD)
      fs_test[, -1] <- centered_newdata
    } else {
      centered_newdata <- NULL
    }
    pred <- predict(cv_m, newdata = fs_test[, -1])
    if (model_type == "classification"){
      accu <- sum(diag(table(pred, fs_test$y))) / length(fs_test$y)  # accuracy = total TP / total (TP: true positive)
      tmp_out <- list(univariate.fs = univariate.fs, uni.sig.fs = uni_sig_fs, selected.features = fs,
                      cv_svm_model = cv_m, nested.cv.accuracy = accu, cv_test_data = fs_test)
    } else {
      error <- pred - fs_test$y
      rmse <- sqrt(mean(error^2))
      rsq <- cor(pred, fs_test$y)
      tmp_out <- list(univariate.fs = univariate.fs, uni.sig.fs = uni_sig_fs, selected.features = fs,
                      cv_svm_model = cv_m, nested.cv.rmse = rmse, nested.cv.rsq = rsq, cv_test_data = fs_test)
    }
    if (verbose) cat("Done!\n")
    # foreach output
    return(tmp_out)
  }

  # computing
  if (verbose){
    cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
    cat(paste0("Data center.scale: ", ifelse(center.scale, " ON\n", " OFF\n")))
    cat(paste0("Univariate reduction: ", ifelse(univariate.fs, " ON\n", " OFF\n")))
    cat("\n")
    cat("Nested cross-validation with feature selection (speed depending on hardware configuration)...\n")
  }
  nested.cv.list <- vector(mode = "list", length = cross.k)
  nested.cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS")) %do% nestedcv_func(i)
  names(nested.cv.list) <- paste0("cv_fold_", c(1:cross.k))

  # below: cv.model.idx: best models index
  if (model_type == "classification"){
    nested.accu <- foreach(i = 1:cross.k, .combine = "c") %do% {
      nested.cv.list[[i]]$nested.cv.accuracy
    }
    nested.rmse <- NULL
    nested.rsq <- NULL
    cv.median <- median(nested.accu)
    cv.model.idx <- which(nested.accu >= cv.median)  # the higher the better
    best.nested.accu <- nested.accu[cv.model.idx]
    best.nested.rmse <- NULL
    best.nested.rsq <- NULL
  } else {
    nested.accu <- NULL
    nested.rmse<- foreach(i = 1:cross.k, .combine = "c") %do% {
      nested.cv.list[[i]]$nested.cv.rmse
    }
    nested.rsq <- foreach(i = 1:cross.k, .combine = "c") %do% {
      nested.cv.list[[i]]$nested.cv.rsq
    }
    cv.median <- median(nested.rmse)
    cv.model.idx <- which(nested.rmse <= cv.median)  # the lower the better

    best.nested.accu <- NULL
    best.nested.rmse <- nested.rmse[cv.model.idx]
    best.nested.rsq <- nested.rsq[cv.model.idx]
  }

  if (cross.best.model.method == "median"){
    final.cv.list <- nested.cv.list[cv.model.idx]
  } else {
    final.cv.list <- nested.cv.list
    cv.model.idx <- seq(cross.k)
  }

  nested.fs <- foreach(i = 1:length(nested.cv.list), .combine = "c") %do% {
    nested.cv.list[[i]]$selected.features
  }
  fs.count <- sort(table(nested.fs), decreasing = TRUE)

  # we only use best.fs.count for thresholding becuase it is going to be the same as nest.fs when cross.best.model.method = "none"
  best.nested.fs <- foreach(i = 1:length(final.cv.list), .combine = "c") %do% {
    final.cv.list[[i]]$selected.features
  }
  best.nested.fs.count <- sort(table(best.nested.fs), decreasing = TRUE)

  if (cross.best.model.method == "median"){  # adjust the threshold to 1 if only one model existed
    if (!any(fs.count >= fs.count.cutoff)) {
      fs.count.cutoff <-  1
    }
  }

  # end time
  end_time <- Sys.time()

  ## output
  if (model_type == "classification"){
    tot.nested.acc.summary <- c(mean(nested.accu), sd(nested.accu), sd(nested.accu)/sqrt(cross.k))
    names(tot.nested.acc.summary) <- c("tot.nested.accuracy", "sd", "sem")
    best.nested.acc.summary <- c(mean(best.nested.accu), sd(best.nested.accu), sd(best.nested.accu)/sqrt(length(cv.model.idx)))
    names(best.nested.acc.summary) <- c("best.nested.acc.summary", "sd", "sem")
    tot.nested.rmse.summary <- NULL
    tot.nested.rsq.summary <- NULL
    best.nested.rmse.summary <- NULL
    best.nested.rsq.summary <- NULL
  } else {
    tot.nested.acc.summary <- NULL
    best.nested.acc.summary <- NULL

    tot.nested.rmse.summary <- c(mean(nested.rmse), sd(nested.rmse), sd(nested.rmse)/sqrt(cross.k))
    names(tot.nested.rmse.summary) <- c("tot.nested.RMSE", "sd", "sem")
    best.nested.rmse.summary <- c(mean(best.nested.rmse), sd(best.nested.rmse), sd(best.nested.rmse)/sqrt(length(cv.model.idx)))
    names(best.nested.rmse.summary) <- c("best.nested.RMSE.summary", "sd", "sem")

    tot.nested.rsq.summary <- c(mean(nested.rsq), sd(nested.rsq), sd(nested.rsq)/sqrt(cross.k))
    names(tot.nested.rsq.summary) <- c("tot.nested.rsq", "sd", "sem")
    best.nested.rsq.summary <- c(mean(best.nested.rsq), sd(best.nested.rsq), sd(best.nested.rsq)/sqrt(length(cv.model.idx)))
    names(best.nested.rsq.summary) <- c("best.nested.rsq.summary", "sd", "sem")
  }
  selected.features <- names(best.nested.fs.count[which(best.nested.fs.count >= fs.count.cutoff)])

  # display
  if (verbose) {
    cat("\n")
    cat("SVM model type: ")
    cat(model_type)
    if (model_type == "classification"){
      cat("\n\n")
      cat("Best cross-validation accuracy summary: \n")
      print(best.nested.acc.summary)
    } else {
      cat("\n\n")
      cat("Best cross-validation RMSE summary: \n")
      print(best.nested.rmse.summary)
    }
    cat("\n")
    cat("Final (best) cross-validation selected features: \n")
    cat(selected.features)
  }

  # run time
  runtime <- end_time - start_time

  # export to environment
  out <- list(univariate.fs = univariate.fs,
              cv.fold = fold,
              randomized.sample.index = random_sample_idx,
              model.type = model_type,
              tot.nested.accuracy.summary = tot.nested.acc.summary,
              tot.nested.RMSE.summary = tot.nested.rmse.summary,
              tot.nested.rsq.summary = tot.nested.rsq.summary,
              nested.accuracy = nested.accu,
              nested.RMSE = nested.rmse,
              nested.rsq = nested.rsq,
              nested.fs.count = fs.count,
              nested.cv.models = nested.cv.list,
              best.nested.method = cross.best.model.method,
              best.nested.index = cv.model.idx,
              best.nested.accuracy.summary = best.nested.acc.summary,
              best.nested.rmse.summary = best.nested.acc.summary,
              best.nested.rsq.summary = best.nested.rsq.summary,
              best.nested.accuracy = best.nested.accu,
              best.nested.rmse = best.nested.rmse,
              best.nested.rsq = best.nested.rsq,
              best.nested.fs.count = best.nested.fs.count,
              fs.method = fs.method,
              fs.count.threshold = fs.count.cutoff,
              selected.features = selected.features,
              tune.method = tune.method,
              tune.cross.k = if(tune.method == "cross") tune.cross.k else NULL,
              tune.boot.n = if(tune.method == "boot") tune.boot.n else NULL,
              run.time = paste0(signif(runtime[[1]], 4), " ", attributes(runtime)[2]))
  class(out) <- "rbiosvm_nestedcv"
  return(out)
}

#' @export
print.rbiosvm_nestedcv <- function(x, ...){
  cat("SVM model type: ")
  cat(x$model.type)
  cat("\n")
  cat("Univariate reduction: ")
  cat(x$univariate.fs)
  cat("\n")
  if (x$model.type == "classification") {
    cat("Total nested cross-validation accuracy:\n")
    print(x$tot.nested.accuracy.summary)
    cat("\n")
    cat("Best nested cross-validation accuracy:\n")
    print(x$best.nested.accuracy.summary)
  } else {
    cat("Total nested cross-validation RMSE:\n")
    print(x$tot.nested.RMSE.summary)
    cat("\n")
    cat("Best nested cross-validation RMSE:\n")
    print(x$best.nested.RMSE.summary)
  }
  cat("\n")
  cat(paste0("Consensus selected features (count threshold: ", x$fs.count.threshold,"):", "\n"))
  cat(x$selected.features)
  cat("\n\n")
  cat("Nested CV run time: ")
  cat(x$run.time)
}


#' @title rbioClass_svm_roc_auc()
#'
#' @description ROC-AUC analysis and ploting for SVM model
#' @param object A \code{rbiosvm} object.
#' @param newdata A data matrix or vector for test data. Make sure it is a \code{matrix} or \code{vector} without labels, as well as the same feature numbers as the training set.
#' @param newdata.y For regression model only, the vector for the new data's continuous outcome variable. Default is \code{NULL}
#' @param y.threshold For regression model only, a numeric vector as the theshold(s) to categorize the continuous outcome variable into groups. Default is \code{if (object$model.type == "regression") round(median(object$inputY)) else NULL}
#' @param newdata.label For classification model only, the correspoding label vector to the data. Make sure it is a \code{factor} object. Defaults is \code{NULL}.
#' @param center.scale.newdata Logical, wether center and scale the newdata with training data mean and standard deviation. Default is \code{TRUE}.
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
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
#'          When \code{newdata} is not provided, the function uses the training data from the input SVM object.
#'
#'          Although optional, the \code{newdata} matrix should use training data's column mean and column standard deviation to center.scale prior to ROC-AUC analysis.
#'          The option \code{center.scaled.newdata = FALSE} is used when the whole (training and test sets) data were center.scaled before SVM training and testing.
#'
#'          For regression models, it is ok to provide multiple thresholds for \code{y.threshold}.
#'
#'          Along with the data frame, the \code{svm_roc_auc} output includes the \code{roc} object from the \code{pROC} package, with which the stats can be done to compare ROCs.
#'
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
rbioClass_svm_roc_auc <- function(object, newdata = NULL, newdata.label = NULL,
                                  newdata.y = NULL,
                                  y.threshold = if (object$model.type == "regression") round(median(object$inputY)) else NULL,
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
  if (is.null(newdata)) {
    cat("No newdata input, proceed with training data.\n\n")
    newdata <- object$inputX
    if (object$model.type == "classification") {
      newdata.label <- object$inputY
    } else {
      newdata.y <- object$inputY
    }
  }
  if (!class(newdata) %in% c("data.frame", "matrix") & !is.null(dim(newdata))) stop("newdata needs to be a matrix, data.frame or vector.")
  if (class(newdata) == "data.frame" | is.null(dim(newdata))){
    if (verbose) cat("newdata converted to a matrix object.\n")
    newdata <- as.matrix(sapply(newdata, as.numeric))
  }
  if (ncol(newdata) != ncol(object$inputX)) stop("test data should have the same number of variables as the training data.")
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data while center.scale.newdata = TRUE.")
  }
  y <- newdata.y
  if (object$model.type == "classification"){
    if (verbose && !is.null(newdata.y)) {
      cat("newdata.y set to NULL for classification model. \n")
      newdata.y <- NULL
    }
    y <- NULL
    if (class(newdata.label) != "factor"){
      if (verbose) cat("newdata.label is converted to factor. \n")
      newdata.label <- factor(newdata.label, levels = unique(newdata.label))
    }
    reg.group <- NULL
  } else {  # for regression
    if (is.null(newdata.y) || length(newdata.y) != nrow(newdata)) stop("Please set the correct outcome value vector newdata.y for regression model.")
    if (is.null(y.threshold)) stop("Please set an appropriate value for y.threhold for regression model.")
    if (!is.numeric(y.threshold)) stop("y.threshold only takes a numeric vector.")
    y.threshold <- sort(unique(y.threshold))
    threshold.length <- length(y.threshold)
    reg.group.names <- paste0("case", seq(threshold.length + 1))
    rg <- c(0, y.threshold, ceiling(max(newdata.y)*1.1))
    newdata.label <- cut(y, rg)
    reg.group <- levels(newdata.label)
    names(reg.group) <- reg.group.names
    levels(newdata.label) <- reg.group.names
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
  if (object$model.type == "regression"){
    pred <- cut(pred, rg, labels = reg.group.names)
  }
  outcome <- newdata.label  # origial label

  # auc
  roc_auc_list <- vector(mode = "list", length = length(levels(outcome)))
  roc_auc_list[] <- foreach(i = 1:length(levels(outcome))) %do% {
    response <- outcome
    levels(response)[-i] <- "others"
    predictor <- dummy(pred)
    predictor <- as.matrix(predictor[, i], ncol = 1)
    splt <- split(predictor, response)  # split function splist array according to a factor
    controls <- splt$others
    cases <- splt[[levels(outcome)[i]]]
    perf <- tryCatch(pROC::roc(controls = controls, cases = cases, smooth = plot.smooth, ci = TRUE),
                     error = function(err){
                       cat("Curve not smoothable. Proceed without smooth.\n")
                       pROC::roc(controls = controls, cases = cases, smooth = FALSE, ci = TRUE)
                     })
    if (length(levels(outcome)) == 2){
      cat(paste0("AUC - ", levels(outcome)[i], ": ", perf$auc, "\n"))
    } else {
      cat(paste0(" AUC - ", levels(outcome)[i], " (vs Others): ", perf$auc, "\n"))
    }
    perf
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
  out <- list(model.type = object$model.type,
              y.threshold = y.threshold,
              regression.categories = reg.group,
              svm.roc_object = roc_auc_list,
              svm.roc_dataframe = roc_dfm,
              input.newdata = newdata,
              input.newdata.y = newdata.y,
              input.newdata.label = newdata.label,
              newdata.center.scaled = centered_newdata)
  class(out) <- "svm_roc_auc"
  assign(paste(deparse(substitute(object)), "_svm_roc_auc", sep = ""), out, envir = .GlobalEnv)

  ## plotting
  if (rocplot){
    if (verbose) cat(paste("Plot being saved to file: ", deparse(substitute(object)),".svm.roc.pdf...", sep = ""))  # initial message

    plt <- ggplot(data = roc_dfm, aes(x = fpr, y = tpr, group = group, colour = group)) +
      geom_line(aes(linetype = group), size = plot.lineSize) +
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
    # grid.newpage()
    ggsave(filename = paste(deparse(substitute(object)),".svm.roc.pdf", sep = ""), plot = plt,
           width = plot.Width, height = plot.Height, units = "mm",dpi = 600)
    grid.draw(plt)
    if (verbose) cat("Done!\n")
  }
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
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
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
    perm_model <- rbioClass_svm(x = object$center.scaledX$centerX, y = perm_y,
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
  by_feature_per_y_func <- function(i){
    set.seed(i)
    perm_x <- foreach(m = unique(levels(object$inputY)), .combine = "rbind") %do% {
      sub_dat <- object$center.scaledX$centerX[which(object$inputY == m), ]  # subsetting the centre-scaled X by label (Y)
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


#' @title rbioClass_svm_predict
#'
#' @description Prediction function for SVM analysis. The function calculates the predicted value for unknown sample data using the input SVM model.
#' @param object A \code{rbiosvm} object.
#' @param newdata Input data to be classified. Make sure it is a \code{matrix} class and has the same variables as the model, i.e. same number of columns as the training data.
#' @param sampleLabel.vector A character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param center.scale.newdata If to center the newdata. When \code{TRUE}, it will also apply the same scaling option as the \code{object}. Default is \code{TRUE}.
#' @param newdata.y Only requried for regression study. The default is \code{NULL}.
#' @param prob.method Method to calculate classification probability. Options are \code{"logistic"}, \code{"softmax"} and \code{"Bayes"}. See details for more information. Default is \code{"logistic"}.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will not affect error or warning messeages.
#' @return  A \code{prediction} obejct. The items of the object are:
#'
#' \code{model.type}
#'
#' \code{classifier.class}
#'
#' \code{predited.value}
#'
#' \code{tot.predict.RMSE}
#'
#' \code{prob.method}: Method to caluculate posterier probability
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
#' @details Although optional, the \code{newdata} matrix should be centered prior to testing, with the same scaling setting as the input \code{rbiosvm} object.
#'          The option \code{center.scale.newdata = FALSE} is for the already centered the data matrix. This center.scale process should use training data's
#'          column mean and column standard deviation.
#'
#'          The default posterior probability calculation method \code{"logistic"} is the \code{e1071} pacakge's implementation of logistic regression model.
#'          See \code{\link{rbioClass_plsda_predict}} for description for "Bayes" and "softmax" method.
#'
#'          If \code{sampleLabel.vector = NULL} or missing, the function uses row numbers as label.
#'
#'          When the model type is \code{"regression"}, the value of the irrelavent items is set to \code{NULL}. Likewise, when the model type is \code{"classification"},
#'          \code{newdata.y} and \code{tot.predict.RMSE} are set to \code{NULL}
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
rbioClass_svm_predcit <- function(object, newdata, center.scale.newdata = TRUE,
                                  sampleLabel.vector = NULL, newdata.y = NULL,
                                  prob.method = c("logistic",  "Bayes", "softmax"),
                                  verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c("rbiosvm"))) stop("object needs to be a \"rbiosvm\" class.")
  if (missing(newdata) || is.null(newdata)) stop("please provide newdata.")
  if (!class(newdata) %in% c("matrix", "data.frame")) stop("newdata has to be either a matrix or data.frame object.")
  if (ncol(newdata) != ncol(object$inputX)) stop("newdata needs to have the same number of variables, i.e. columns, as the object.")
  # if (!prob.method %in% c("logistic",  "Bayes", "softmax")) stop("Probability method should be either \"softmax\" or \"Bayes\".")
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data while center.scale.newdata = TRUE.")
  }
  if (!missing(sampleLabel.vector) & !is.null(sampleLabel.vector)){
    if (length(sampleLabel.vector) != nrow(newdata)){
      cat("Sample label vector not the same length as newdata. Proceed without custom sample labels.\n")
      sampleLabel.vector <- NULL
    }
  }
  if (object$model.type == "classification"){
    newdata.y <- NULL
    y <- newdata.y
    prob.method <- match.arg(prob.method, c("logistic",  "Bayes", "softmax"))
    probability <- TRUE
  } else {
    if (missing(newdata.y) || is.null(newdata.y) || length(newdata.y) != nrow(newdata)) stop("Please set the correct outcome value vector newdata.y for regression study.")
    y <- newdata.y
    prob.method <- NULL
    probability <- FALSE
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
  pred <- predict(object = object, newdata = test, probability = probability)

  if (object$model.type == "classification"){
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
    tot.pred.rmse <- NULL
  } else {  # for regression study
    prob_dfm <- NULL
    error <- pred - y
    tot.pred.rmse <- sqrt(mean(error^2))
  }

  ## output
  predicted.value <- pred
  attributes(predicted.value) <- NULL
  out <- list(model.type = object$model.type,
              classifier.class = class(object),
              predicted.value = predicted.value,
              tot.predict.RMSE = tot.pred.rmse,
              prob.method = prob.method,
              probability.summary = prob_dfm,
              raw.newdata = newdata,
              center.scale = center.scale.newdata,
              center.scaled.newdata = centerdata,
              newdata.y = y)
  class(out) <- "prediction"
  assign(paste(deparse(substitute(object)), "_svm_predict", sep = ""), out, envir = .GlobalEnv)
}



#' @export
print.prediction <- function(x, ...){
  cat("Model type: ", x$model.type, "\n")
  cat("\n")
  if (x$model.type == "classification"){
    cat("Prediction results:\n")
    print(x$probability.summary[, c(1:2, 5)])
  } else {
    cat("Total RMSE on test set: ")
    cat(x$tot.predict.RMSE)
    cat("\n")
  }
  cat("\n")
  cat(paste0("Classifier class:\n"))
  print(x$classifier.class)
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
#' The R2 is calculated as follwoing: 1 - rss/tss.
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
