#' @title rbioClass_svm
#'
#' @description Support Vector Machine (SVM) modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param center.scale Logical, whether center and scale the data, i.e. subtracting mean (col.mean) and deviding by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param svm.cross.k Fold of cross validation. Default is \code{10}.
#' @param tune.method Parameter tuning method. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{svm} function from \code{e1071} pacakge.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
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
#' If \code{svm.cross.k} is set to 0, no cross validation is used for the final modelling.
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
    # if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
    model_type <- "classification"
    y <- factor(y, levels = unique(y))
  } else {
    model_type <- "regression"
  }
  if (!any(class(x) %in% c("data.frame", "matrix")) & !is.null(dim(x))) stop("x needs to be a matrix, data.frame or vector.")
  if (any(class(x) == "data.frame") | is.vector(x)){
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
#' @description Nested cross-validation feature selection and SVM modelling: CV-rRF-FS-SVM.
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param univariate.fs If to use limma-based univariate reduction. Default is \code{FALSE}.
#' @param uni.log2trans Only set if \code{univariate.fs = TRUE}, if to log2 transform data before univariate reduction.
#' @param uni.contrast Only set if \code{univariate.fs = TRUE} and for a classification study, the contrast for the univariate analysis. Default is \code{NULL}.
#' @param uni.alpha Only set if \code{univariate.fs = TRUE}, the p value alpha. Default is \code{0.05}.
#' @param uni.fdr Only set if \code{univariate.fs = TRUE}, if to use FDR for the p value. Default is \code{FALSE}.
#' @param center.scale Logical, whether center and scale the data, i.e. subtracting mean (col.mean) and dividing by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param cross.k Fold of nested cross validation, i.e. outer loop. Default is \code{10}.
#' @param cross.best.model.method The method to select the best cv models for feature selection. Options are \code{"median"} and \code{"none"}. Default is \code{"median"}.
#' @param tune.method Parameter tuning method, i.e. inner loop. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{rbioClass_svm}.
#' @param fs.method Feature selection method. Only \code{"rf"} (i.e. random forest) is supported so far. Default is \code{"rf"}.
#' @param rf.ifs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the initial feature selection step of recursive feature selection. Default is \code{501}.
#' @param rf.sfs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the sequential forward selection step of recursive feature selection. Default is \code{501}.
#' @param fs.count.cutoff A integer for feature vote cutoff. Default is outer loop cross-validation fold \code{cross.k}.
#' @param parallelComputing Whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Returns a nested CV SVM model object: \code{rbiosvm_nestedcv}.
#'
#' Additional items for \code{rbiosvm_nestedcv}:
#'
#' \code{univariate.fs}
#'
#' \code{cv.fold}: number of (outer) cross-validation fold
#'
#' \code{randomized.sample.index}: randomized sample order for (outer) cross-validation. \code{NULL} for classification modelling (cv.fold already randomized).
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
#' \code{best.nested.method}: method to determine the best nested model.
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
#' For now, RBioFS implementation of two-step random forest feature selection (rRF-FS) is used to select features, and based on nested cross-validation with SVM.
#' Resulted features from each nested cross-validation round are voted. Features with votes equal or greater than the cutoff are reported as selected features.
#'
#' It is also a good idea to set \code{fs.count.cutoff} as \code{cross.k - 1}. Notably, \code{fs.count.cutoff = 1} is "no threshold",
#' meaning maximum number of features selected from the nested cross-validation are reported.
#'
#' When \code{cross.best.model.method = "median"}, the function only use models with accuracy/RMSE equal or better than the median value
#' for feature count threholding. When there is no change in performance across cv models, the function behaves same as \code{cross.best.model.method = "none"}
#'
#' For parallel computing, the nested CV rRF-FS part uses parallel for the rRF-FS functions \code{rbioFS_rf_initialFS} and \code{rbioFS_rf_SFS},
#' i.e. the \code{parallelComputing} argument in these two functions. This would lead to as many CPU threads as set likely being used due to the number of trees
#' often exceeding the number of threads. However, the SVM process is paralleled per CV iteration, e.g. 10-fold CV would use 10 CPU threads, which may lead to
#' insufficient CPU thread utilization. This problem would not likely be solved unless the SVM dependency \code{e1071} is updated with parallel computing,
#' or a different SVM dependency that has parallel computing is used.
#'
#' The function also supports regression study, in which case, the performance metric is \code{RMSE}.
#'
#' For classification modelling, the 10-fold CV segmentation is stratified, whereas the "normal" random resampling is used for regression modelling.
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
    # if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
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

  if (cross.k > nrow(x)) stop("Nested cross-validation fold setting cross.k exceeded limit. Hint: max at total sample number.\n")
  if (any(class(x) == "data.frame")){
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

  ## nested cv-fs
  # processing training data
  dfm <- data.frame(y, x, check.names = FALSE)
  if (model_type == "regression"){
    random_sample_idx <- sample(nrow(dfm))
    dfm_randomized <- dfm[random_sample_idx, ]  # randomize samples
    fold <- cut(1:nrow(dfm_randomized), breaks = cross.k, labels = FALSE)  # create a factor object with fold indices
  } else {
    dfm_randomized <- dfm  # no randomization at this point yet
    y_summary <- table(y)
    fold <- foreach(i = 1:length(levels(y)), .combine = "c") %do% {
      min_rep <- y_summary[i] %/% cross.k
      if (min_rep > 0){
        remain_size <- y_summary[i] %% cross.k
        v <- rep(1:cross.k, times = min_rep)
        if (remain_size > 0) v <- c(v, sample(1:cross.k, remain_size))
        sample(v)  # randomization
      } else {
        sample(1:cross.k, size = y_summary[i]) # randomization
      }
    }
    random_sample_idx <- NULL
  }

  # nested CV-FS function
  nested_cvfs_func <- function(i) {
    if (verbose) cat(paste0("Iteration: ", i, "|", cross.k, "..."))
    cv_training <- dfm_randomized[which(fold != i, arr.ind = TRUE), ]
    cv_training_x <- cv_training[, -1]
    # below: remove columns with constant value
    cv_training_x <- cv_training_x[,!apply(cv_training_x, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
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

      if (length(uni_sig_fs) == 0) {
        stop("\nNo statistically significant features found. Try runing with larger uni.alpha value, or univarate.fs = FALSE.")
      } else {
        cv_training <- cv_training[, c("y", uni_sig_fs), drop = FALSE]
        cv_training_x <- cv_training[, -1, drop = FALSE]
        cv_training_y <- cv_training[, 1]
      }
    } else {
      uni_sig_fs <- NULL
    }

    # fs
    if (center.scale){
      fs_training_x <- center_scale(cv_training_x, scale = TRUE)$centerX  # without y
    } else {
      fs_training_x <- cv_training_x
    }

    fs <- tryCatch({  # rRF-FS with error handling
      rbioFS_rf_initialFS(objTitle = paste0("svm_nested_iter_", i), x = fs_training_x, y = cv_training_y, nTimes = 50,
                          nTree = rf.ifs.ntree,
                          parallelComputing = parallelComputing, clusterType = clusterType, n_cores = n_cores,
                          plot = FALSE)
      # fs <- svm_nested_initial_FS$feature_initial_FS
      rbioFS_rf_SFS(objTitle = paste0("svm_nested_iter_", i),
                    x = eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$training_initial_FS, y = cv_training_y, nTimes = 50,
                    nTree = rf.sfs.ntree,
                    parallelComputing = parallelComputing, clusterType = clusterType, n_cores = n_cores,
                    plot = FALSE)
      if (verbose) cat("Done!\n")
      if (length(eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features) > 1){
        out <- eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features
      } else {
        out <- eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$feature_initial_FS
      }
    },
    error = function(e){
      cat(paste0("\nErrors occurred during rRF-FS, skipping rRF-FS for CV iteration: ", i, "...\n", "\tError message: ", e, "\n"))
      out <- colnames(fs_training_x)
      return(out)
    })

    # processing test data
    fs_test <- dfm_randomized[which(fold == i, arr.ind = TRUE), ][, c("y", fs)]  # preserve y and selected features

    cv_fs_out <- list(
      fs = fs,
      fs_training_x = fs_training_x,
      cv_training_y = cv_training_y,
      uni_sig_fs = uni_sig_fs,
      fs_test = fs_test
    )
    return(cv_fs_out)
  }

  # nested CV-M function: SVM modelling for nested CV iterations
  nested_cvm_func <- function(i, ...){
    fs <- nested_cvfs.list[[i]]$fs
    fs_training_x <- nested_cvfs.list[[i]]$fs_training_x
    cv_training_y <- nested_cvfs.list[[i]]$cv_training_y
    uni_sig_fs <- nested_cvfs.list$uni_sig_fs

    cv_m <- rbioClass_svm(x = fs_training_x[, fs], y = cv_training_y, center.scale = center.scale,
                          svm.cross.k = 0, tune.method = tune.method, kernel = kernel,
                          tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE, ...)

    # processing test data
    fs_test <- nested_cvfs.list[[i]]$fs_test
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

    cat(paste0("cvm iteration ", i, "...Done!\n"))
    return(tmp_out)
  }


  # computing
  if (verbose){
    cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
    cat(paste0("Data center.scale: ", ifelse(center.scale, " ON\n", " OFF\n")))
    cat(paste0("Univariate reduction: ", ifelse(univariate.fs, " ON\n", " OFF\n")))
    cat("\n")
    cat("Nested cross-validation with feature selection (speed depending on hardware configuration): \n")
    cat("Nested CV rRF-FS: \n")
  }
  nested_cvfs.list <- vector(mode = "list", length = cross.k)
  nested_cvfs.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "stop") %do% nested_cvfs_func(i)
  names(nested_cvfs.list) <- paste0("cv_fold_", c(1:cross.k))

  if (verbose) cat(paste0("SVM for nested CV rRF-FS assessment (", cross.k, " iterations)..."))
  if (parallelComputing) {
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = "FORK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    nested_cvm.list <- vector(mode = "list", length = cross.k)
    nested_cvm.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "stop") %dopar% nested_cvm_func(i, ...)
    names(nested_cvm.list) <- paste0("cv_fold_", c(1:cross.k))
  } else {
    nested_cvm.list <- vector(mode = "list", length = cross.k)
    nested_cvm.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "stop") %do% nested_cvm_func(i, ...)
    names(nested_cvm.list) <- paste0("cv_fold_", c(1:cross.k))
  }
  if (verbose) cat("Done!\n")

  # below: cv.model.idx: best models index
  nested.cv.list <- nested_cvm.list  # set up nested.cv.list instead of using new nested_cvm.list for compatibility
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

  # we only use best.fs.count for thresholding because it is going to be the same as nest.fs when cross.best.model.method = "none"
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
    cat("\n\n")
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


#' @title rbioClass_svm_ncv_fs_v2
#'
#' @description Nested cross-validation feature selection and SVM modelling: CV-rRF-FS-SVM. This is V2 with slightly different parallel computing implementation
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param univariate.fs If to use limma-based univariate reduction. Default is \code{FALSE}.
#' @param uni.log2trans Only set if \code{univariate.fs = TRUE}, if to log2 transform data before univariate reduction.
#' @param uni.contrast Only set if \code{univariate.fs = TRUE} and for a classification study, the contrast for the univariate analysis. Default is \code{NULL}.
#' @param uni.alpha Only set if \code{univariate.fs = TRUE}, the p value alpha. Default is \code{0.05}.
#' @param uni.fdr Only set if \code{univariate.fs = TRUE}, if to use FDR for the p value. Default is \code{FALSE}.
#' @param center.scale Logical, whether center and scale the data, i.e. subtracting mean (col.mean) and dividing by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param cross.k Fold of nested cross validation, i.e. outer loop. Default is \code{10}.
#' @param cross.best.model.method The method to select the best cv models for feature selection. Options are \code{"median"} and \code{"none"}. Default is \code{"median"}.
#' @param tune.method Parameter tuning method, i.e. inner loop. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{rbioClass_svm}.
#' @param fs.method Feature selection method. Only \code{"rf"} (i.e. random forest) is supported so far. Default is \code{"rf"}.
#' @param rf.ifs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the initial feature selection step of recursive feature selection. Default is \code{501}.
#' @param rf.sfs.ntree Set only when \code{fs.method = "rf"}, ntree setting for the sequential forward selection step of recursive feature selection. Default is \code{501}.
#' @param fs.count.cutoff A integer for feature vote cutoff. Default is outer loop cross-validation fold \code{cross.k}.
#' @param parallelComputing Whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose Whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Returns a nested CV SVM model object: \code{rbiosvm_nestedcv}.
#'
#' Additional items for \code{rbiosvm_nestedcv}:
#'
#' \code{univariate.fs}
#'
#' \code{cv.fold}: number of (outer) cross-validation fold
#'
#' \code{randomized.sample.index}: randomized sample order for (outer) cross-validation. \code{NULL} for classification modelling (cv.fold already randomized).
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
#' \code{best.nested.method}: method to determine the best nested model.
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
#' This "v2" implements parallel computing when running nested_cv sub-function, as opposed to the RF-FS function inside the nested_cv like the original
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
#' For classification modelling, the 10-fold CV segmentation is stratified, whereas the "normal" random resampling is used for regression modelling.
#'
#'
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma lmFit topTable makeContrasts eBayes
#' @importFrom splines ns
#' @examples
#' \dontrun{
#' svm_nestedcv <- rbioClass_svm_ncv_fs_v2(x = mydata[, -c(1:2)],
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
rbioClass_svm_ncv_fs_v2 <- function(x, y,
                                    univariate.fs = FALSE, uni.log2trans = TRUE, uni.contrast = NULL,
                                    uni.alpha = 0.05, uni.fdr = FALSE,
                                    center.scale = TRUE,
                                    kernel = c("radial", "linear", "polynomial", "sigmoid"),
                                    cross.k = 10, cross.best.model.method = c("none", "median"),
                                    tune.method = c("cross", "boot", "fix"),
                                    tune.cross.k = 10, tune.boot.n = 10,
                                    fs.method = "rf", rf.ifs.ntree = 1001, rf.sfs.ntree = 1001,
                                    fs.count.cutoff = cross.k,
                                    parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                                    verbose = TRUE){
  ## initiate the run time and set seed
  start_time <- Sys.time()

  ## check arguments
  if (!fs.method %in% c("rf")) stop("So far, fs.method has to be \"rf\". More methods will be implemented")
  if (is.factor(y)) {
    # if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
    model_type <- "classification"
  } else {
    model_type <- "regression"
  }

  if (model_type == "classification" && univariate.fs) {
    if (is.null(uni.contrast)) {
      stop("When univariate.fs = TRUE, uni.contrast needs to be set for classification study.")
    } else if (!is.character(uni.contrast)) {
      stop("uni.contrast needs to be a character string.")
    }
  }

  if (cross.k > nrow(x)) stop("Nested cross-validation fold setting cross.k exceeded limit. Hint: max at total sample number.\n")
  if (any(class(x) == "data.frame")){
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
  if (model_type == "regression"){
    random_sample_idx <- sample(nrow(dfm))
    dfm_randomized <- dfm[random_sample_idx, ]  # randomize samples
    fold <- cut(1:nrow(dfm_randomized), breaks = cross.k, labels = FALSE)  # create a factor object with fold indices
  } else {
    dfm_randomized <- dfm  # no randomization at this point yet
    y_summary <- table(y)
    fold <- foreach(i = 1:length(levels(y)), .combine = "c") %do% {
      min_rep <- y_summary[i] %/% cross.k
      if (min_rep > 0){
        remain_size <- y_summary[i] %% cross.k
        v <- rep(1:cross.k, times = min_rep)
        if (remain_size > 0) v <- c(v, sample(1:cross.k, remain_size))
        sample(v)  # randomization
      } else {
        sample(1:cross.k, size = y_summary[i]) # randomization
      }
    }
    random_sample_idx <- NULL
  }

  # nested CV function
  nestedcv_func <- function(i) {
    # if (verbose) cat(paste0("Nested CV iteration: ", i, "|", cross.k, "..."))
    cv_training <- dfm_randomized[which(fold != i, arr.ind = TRUE), ]
    cv_training_x <- cv_training[, -1]
    # below: remove columns with constant value
    cv_training_x <- cv_training_x[,!apply(cv_training_x, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
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

      if (length(uni_sig_fs) == 0) {
        stop("No statistically significant features found. Try runing with larger uni.alpha value, or univarate.fs = FALSE.")
      } else {
        cv_training <- cv_training[, c("y", uni_sig_fs), drop = FALSE]
        cv_training_x <- cv_training[, -1, drop = FALSE]
        cv_training_y <- cv_training[, 1]
      }
    } else {
      uni_sig_fs <- NULL
    }

    # fs
    if (center.scale){
      fs_training_x <- center_scale(cv_training_x, scale = TRUE)$centerX  # without y
    } else {
      fs_training_x <- cv_training_x
    }

    fs <- tryCatch({  # rRF-FS with error handling
      rbioFS_rf_initialFS(objTitle = paste0("svm_nested_iter_", i), x = fs_training_x, y = cv_training_y, nTimes = 50,
                          nTree = rf.ifs.ntree,
                          parallelComputing = FALSE, clusterType = clusterType, n_cores = n_cores,
                          plot = FALSE)
      rbioFS_rf_SFS(objTitle = paste0("svm_nested_iter_", i),
                    x = eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$training_initial_FS, y = cv_training_y, nTimes = 50,
                    nTree = rf.sfs.ntree,
                    parallelComputing = FALSE, clusterType = clusterType, n_cores = n_cores,
                    plot = FALSE)

      if (length(eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features) > 1){
        out <- eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features
      } else {
        out <- eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$feature_initial_FS
      }
    },
    error = function(e){
      cat(paste0("Errors occurred during rRF-FS, skipping rRF-FS for CV iteration: ", i, "...\n", "\tError message: ", e))
      out <- colnames(fs_training_x)
      return(out)
    })

    cv_m <- rbioClass_svm(x = fs_training_x[, fs], y = cv_training_y, center.scale = center.scale,
                          svm.cross.k = 0, tune.method = tune.method, kernel = kernel,
                          tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE)
    # processing test data
    fs_test <- dfm_randomized[which(fold == i, arr.ind = TRUE), ][, c("y", fs)]  # preserve y and selected features
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
    # if (verbose) cat("Done!\n")
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
  if (parallelComputing) {
    # set up cpu cluster
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # computing
    nested.cv.list <- vector(mode = "list", length = cross.k)
    nested.cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "pass") %dopar% nestedcv_func(i)
    names(nested.cv.list) <- paste0("cv_fold_", c(1:cross.k))
  } else {
    nested.cv.list <- vector(mode = "list", length = cross.k)
    nested.cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "pass") %do% nestedcv_func(i)
    names(nested.cv.list) <- paste0("cv_fold_", c(1:cross.k))
  }

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


#' @title rbioClass_svm_cv
#'
#' @description Cross-validation modelling and assessment for SVM modelling. It evaluates the overall performance of SVM modelling given the training data.
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param center.scale Logical, whether center and scale the data, i.e. subtracting mean (col.mean) and dividing by standard deviation (col.sd). Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param cross.k Fold of nested cross validation, i.e. outer loop. Default is \code{10}.
#' @param cross.best.model.method The method to select the best cv models for feature selection. Options are \code{"median"} and \code{"none"}. Default is \code{"median"}.
#' @param tune.method Parameter tuning method, i.e. inner loop. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{rbioClass_svm}.
#' @param parallelComputing whether to use parallel computing or not. Default is \code{TRUE}.
#' @param n_cores Only set when \code{parallelComputing = TRUE}, the number of CPU cores to use. Default is \code{detectCores() - 1}, or the total number cores minus one.
#' @param clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return Returns a CV SVM model object: \code{rbiosvm_cv}.
#'
#' Items for \code{rbiosvm_cv}:
#'
#'      \code{univariate.fs}
#'
#'      \code{cv.fold}: number of (outer) cross-validation fold
#'
#'      \code{randomized.sample.index}: randomized sample order for (outer) cross-validation. \code{NULL} for classification modelling (cv.fold already randomized).
#'
#'      \code{model.type}: the SVM model type, "classification" or "regression".
#'
#'      \code{tot.accuracy.summary}: total (i.e. mean) nested cross-validation accuracy, if \code{model_type = "classification"}
#'
#'      \code{tot.RMSE.summary}: total (i.e. mean) nested cross-validation RMSE, if \code{model_type = "regression"}
#'
#'      \code{tot.rsq.summary}: total (i.e. mean) nested cross-validation R2, if \code{model_type = "regression"}
#'
#'      \code{accuracy}: accuracy for each cross-validation iteration, if \code{model_type = "classification"}
#'
#'      \code{RMSE}: RMSE for each cross-validation iteration, if \code{model_type = "regression"}
#'
#'      \code{rsq}: R2 for each cross-validation iteration, if \code{model_type = "regression"}
#'
#'      \code{cv.models}: all the CV SVM models
#'
#'      \code{final.cv.models}: the CV SVM models selected by the \code{cross.best.model.method}.
#'                              when \code{cross.best.model.method = NULL}, \code{final.cv.models} is the same as \code{cv.models}.
#'
#'      \code{nested.cv.models}: all the cv models with test data. NOTE: the test data is center-scaled with training data if \code{center.scale = TRUE}.
#'
#'      \code{best.method}: method to determine the best model.
#'
#'      \code{best.index}: index of the best models in the total model list
#'
#'      \code{best.acc.summary}
#'
#'      \code{best.rmse.summary}
#'
#'      \code{best.rsq.summary}
#'
#'      \code{best.accu}
#'
#'      \code{best.rmse}
#'
#'      \code{best.rsq}
#'
#'      \code{tune.method}: (inner) loop method for SVM grid search
#'
#'      \code{tune.cross.k}: fold number for (inner) loop if cross-validation is chosen
#'
#'      \code{tune.boot.n}: iteration number for (inner) loop if bootstrap is chosen
#'
#'      \code{run.time}: total run time
#'
#' @details
#'
#' This function is very similar to \code{\link{rbioClass_svm_ncv_fs}}.
#' However, this function is intended for evaluating the final model, whereas the latter is for feature selection.
#'
#' When \code{cross.best.model.method = "median"}, the function only use models with accuracy/RMSE equal or better than the median valaue
#' for feature count threshold. When there is no change in performance across cv models, the function behaves same as \code{cross.best.model.method = "none"}
#'
#' The function also supports regression study, in which case, the performance metric is \code{RMSE}.
#'
#' For classification modelling, the 10-fold CV segmentation is stratified, whereas the "normal" random resampling is used for regression modelling.
#'
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom limma lmFit topTable makeContrasts eBayes
#' @importFrom splines ns
#' @examples
#' \dontrun{
#' svm_cv <- rbioClass_svm_cv(x = mydata[, -c(1:2)],
#'                                      y = factor(mydata$Conditions, levels = unique(mydata$Conditions)),
#'                                      center.scale = TRUE,
#'                                      cross.k = 5,
#'                                      tune.method = "cross",
#'                                      tune.cross.k = 5,
#'                                      parallelComputing = TRUE, clusterType = "FORK")
#'
#' }
#' @export
rbioClass_svm_cv <- function(x, y,
                             center.scale = TRUE,
                             kernel = c("radial", "linear", "polynomial", "sigmoid"),
                             cross.k = 10, cross.best.model.method = c("none", "median"),
                             tune.method = c("cross", "boot", "fix"),
                             tune.cross.k = 10, tune.boot.n = 10, ...,
                             parallelComputing = TRUE, n_cores = parallel::detectCores() - 1, clusterType = c("PSOCK", "FORK"),
                             verbose = TRUE){
  ## initiate the run time and set seed
  start_time <- Sys.time()

  ## check arguments
  if (is.factor(y)) {
    # if (nlevels(y) > 3) warning("y has more than three groups. SVM is not recommended.\n")
    model_type <- "classification"
  } else {
    model_type <- "regression"
  }

  if (cross.k > nrow(x)) stop("Cross-validation fold setting cross.k exceeded limit. Hint: max at total sample number.\n")
  if (any(class(x) == "data.frame")){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    x <- as.matrix(sapply(x, as.numeric))
  }

  kernel <- match.arg(tolower(kernel), c("radial", "linear", "polynomial", "sigmoid"))
  cross.best.model.method <- match.arg(tolower(cross.best.model.method), c("none", "median"))
  tune.method <- match.arg(tolower(tune.method), c("cross", "boot", "fix"))
  if (parallelComputing){
    clusterType <- match.arg(clusterType, c("PSOCK", "FORK"))
  }

  ## nested cv
  # processing training data
  dfm <- data.frame(y, x, check.names = FALSE)
  if (model_type == "regression"){
    random_sample_idx <- sample(nrow(dfm))
    dfm_randomized <- dfm[random_sample_idx, ]  # randomize samples
    fold <- cut(1:nrow(dfm_randomized), breaks = cross.k, labels = FALSE)  # create a factor object with fold indices
  } else {
    dfm_randomized <- dfm  # no randomization at this point yet
    y_summary <- table(y)
    fold <- foreach(i = 1:length(levels(y)), .combine = "c") %do% {
      min_rep <- y_summary[i] %/% cross.k
      if (min_rep > 0){
        remain_size <- y_summary[i] %% cross.k
        v <- rep(1:cross.k, times = min_rep)
        if (remain_size > 0) v <- c(v, sample(1:cross.k, remain_size))
        sample(v)  # randomization
      } else {
        sample(1:cross.k, size = y_summary[i]) # randomization
      }
    }
    random_sample_idx <- NULL
  }

  # nested CV function
  cv_func <- function(i, ...) {
    if (verbose) cat(paste0("CV iteration: ", i, "|", cross.k, "..."))
    cv_training <- dfm_randomized[which(fold != i, arr.ind = TRUE), ]
    cv_training_x <- cv_training[, -1]
    cv_training_y <- cv_training[, 1]

    # fs
    if (center.scale){
      cv_training_x <- center_scale(cv_training_x, scale = TRUE)$centerX  # without y
    } else {
      cv_training_x <- cv_training_x
    }

    cv_m <- rbioClass_svm(x = cv_training_x, y = cv_training_y, center.scale = center.scale,
                    svm.cross.k = 0, tune.method = tune.method, kernel = kernel,
                    tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE, ...)

    # processing test data
    cv_test <- dfm_randomized[which(fold == i, arr.ind = TRUE), ]

    if (center.scale){ # using training data mean and sd
      centered_newdata <- t((t(cv_test[, -1]) - cv_m$center.scaledX$meanX) / cv_m$center.scaledX$columnSD)
      cv_test[, -1] <- centered_newdata
    } else {
      centered_newdata <- NULL
    }
    pred <- predict(cv_m, newdata = cv_test[, -1])
    if (model_type == "classification"){
      cv_test$y <- droplevels(cv_test$y)
      accu <- sum(diag(table(pred, cv_test$y))) / length(cv_test$y)  # accuracy = total TP / total (TP: true positive)
      tmp_out <- list(cv_svm_model = cv_m, cv.accuracy = accu, cv_test_data = cv_test)
    } else {
      error <- pred - cv_test$y
      rmse <- sqrt(mean(error^2))
      rsq <- cor(pred, cv_test$y)
      tmp_out <- list(cv_svm_model = cv_m, cv.rmse = rmse, cv.rsq = rsq, cv_test_data = cv_test)
    }
    if (verbose) cat("Done!\n")
    # foreach output
    return(tmp_out)
  }

  # computing
  if (verbose){
    cat(paste0("Parallel computing:", ifelse(parallelComputing, " ON\n", " OFF\n")))
    cat(paste0("Data center.scale: ", ifelse(center.scale, " ON\n", " OFF\n")))
    cat("\n")
    cat("Cross-validation (speed depending on hardware configuration)...\n")
  }
  cv.list <- vector(mode = "list", length = cross.k)
  if (parallelComputing) {
    # set up cpu cluster
    n_cores <- n_cores
    cl <- makeCluster(n_cores, type = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # run function
    cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "pass") %dopar% cv_func(i)
  } else {
    cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "pass") %do% cv_func(i)
  }
  names(cv.list) <- paste0("cv_fold_", c(1:cross.k))

  # below: cv.model.idx: best models index
  if (model_type == "classification"){
    accu <- foreach(i = 1:cross.k, .combine = "c") %do% {
      cv.list[[i]]$cv.accuracy
    }
    rmse <- NULL
    rsq <- NULL
    cv.median <- median(accu)
    cv.model.idx <- which(accu >= cv.median)  # the higher the better
    best.accu <- accu[cv.model.idx]
    best.rmse <- NULL
    best.rsq <- NULL
  } else {
    accu <- NULL
    rmse<- foreach(i = 1:cross.k, .combine = "c") %do% {
      cv.list[[i]]$cv.rmse
    }
    rsq <- foreach(i = 1:cross.k, .combine = "c") %do% {
      cv.list[[i]]$cv.rsq
    }
    cv.median <- median(rmse)
    cv.model.idx <- which(rmse <= cv.median)  # the lower the better

    best.accu <- NULL
    best.rmse <- rmse[cv.model.idx]
    best.rsq <- rsq[cv.model.idx]
  }

  if (cross.best.model.method == "median"){
    final.cv.list <- cv.list[cv.model.idx]
  } else {
    final.cv.list <- cv.list
    cv.model.idx <- seq(cross.k)
  }

  # end time
  end_time <- Sys.time()

  ## output
  if (model_type == "classification"){
    tot.acc.summary <- c(mean(accu), sd(accu), sd(accu)/sqrt(cross.k))
    names(tot.acc.summary) <- c("tot.accuracy", "sd", "sem")
    best.acc.summary <- c(mean(best.accu), sd(best.accu), sd(best.accu)/sqrt(length(cv.model.idx)))
    names(best.acc.summary) <- c("best.acc.summary", "sd", "sem")
    tot.rmse.summary <- NULL
    tot.rsq.summary <- NULL
    best.rmse.summary <- NULL
    best.rsq.summary <- NULL
  } else {
    tot.acc.summary <- NULL
    best.acc.summary <- NULL

    tot.rmse.summary <- c(mean(rmse), sd(rmse), sd(rmse)/sqrt(cross.k))
    names(tot.rmse.summary) <- c("tot.RMSE", "sd", "sem")
    best.rmse.summary <- c(mean(best.rmse), sd(best.rmse), sd(best.rmse)/sqrt(length(cv.model.idx)))
    names(best.rmse.summary) <- c("best.RMSE.summary", "sd", "sem")

    tot.rsq.summary <- c(mean(rsq), sd(rsq), sd(rsq)/sqrt(cross.k))
    names(tot.rsq.summary) <- c("tot.rsq", "sd", "sem")
    best.rsq.summary <- c(mean(best.rsq), sd(best.rsq), sd(best.rsq)/sqrt(length(cv.model.idx)))
    names(best.rsq.summary) <- c("best.rsq.summary", "sd", "sem")
  }

  # display
  if (verbose) {
    cat("\n")
    cat("SVM model type: ")
    cat(model_type)
    if (model_type == "classification"){
      cat("\n\n")
      cat("Best cross-validation accuracy summary: \n")
      print(best.acc.summary)
    } else {
      cat("\n\n")
      cat("Best cross-validation RMSE summary: \n")
      print(best.rmse.summary)
    }
  }

  # run time
  runtime <- end_time - start_time

  # export to environment
  out <- list(cv.fold = fold,
              randomized.sample.index = random_sample_idx,
              model.type = model_type,
              tot.accuracy.summary = tot.acc.summary,
              tot.RMSE.summary = tot.rmse.summary,
              tot.rsq.summary = tot.rsq.summary,
              accuracy = accu,
              RMSE = rmse,
              rsq = rsq,
              cv.models = cv.list,
              final.cv.models = final.cv.list,
              best.method = cross.best.model.method,
              best.index = cv.model.idx,
              best.accuracy.summary = best.acc.summary,
              best.rmse.summary = best.acc.summary,
              best.rsq.summary = best.rsq.summary,
              best.accuracy = best.accu,
              best.rmse = best.rmse,
              best.rsq = best.rsq,
              tune.method = tune.method,
              tune.cross.k = if(tune.method == "cross") tune.cross.k else NULL,
              tune.boot.n = if(tune.method == "boot") tune.boot.n else NULL,
              run.time = paste0(signif(runtime[[1]], 4), " ", attributes(runtime)[2]))
  class(out) <- "rbiosvm_cv"
  return(out)
}


#' @title rbioClass_svm_predict
#'
#' @description Prediction function for SVM analysis. The function calculates the predicted value for unknown sample data using the input SVM model.
#' @param object A \code{rbiosvm} object.
#' @param newdata Input data to be classified. Make sure it is a \code{matrix} class and has the same variables as the model, i.e. same number of columns as the training data.
#' @param export.name String. Optional user defined export name prefix. Default is \code{NULL}.
#' @param sampleID.vector A character vector containing annotation (i.e. labels) for the samples. Default is \code{NULL}.
#' @param center.scale.newdata If to center the newdata. When \code{TRUE}, it will also apply the same scaling option as the \code{object}. Default is \code{TRUE}.
#' @param newdata.y Only required for regression study. The default is \code{NULL}.
#' @param prob.method Method to calculate classification probability. Options are \code{"logistic"}, \code{"softmax"} and \code{"Bayes"}. See details for more information. Default is \code{"logistic"}.
#' @param verbose whether to display messages. Default is \code{TRUE}. This will not affect error or warning messages.
#' @return  A \code{prediction} object. The items of the object are:
#'
#' \code{model.type}
#'
#' \code{classifier.class}
#'
#' \code{predited.value}
#'
#' \code{tot.predict.RMSE}
#'
#' \code{prob.method}: Method to calculate posterior probability
#'
#' \code{probability.summary}
#'
#' \code{raw.newdata}
#'
#' \code{newdata.id}
#'
#' \code{center.scale}
#'
#' \code{center.scaled.newdata}
#'
#' \code{newdata.y}
#'
#' @details Although optional, the \code{newdata} matrix should be centered prior to testing, with the same scaling setting as the input \code{rbiosvm} object.
#'
#'          The option \code{center.scale.newdata = FALSE} is for the already centered the data matrix. This center.scale process should use training data's
#'          column mean and column standard deviation.
#'
#'          The default posterior probability calculation method \code{"logistic"} is the \code{e1071} package's implementation of logistic regression model.
#'          See \code{\link{rbioClass_plsda_predict}} for description for "Bayes" and "softmax" method.
#'
#'          If \code{sampleID.vector = NULL} or missing, the function uses row numbers as label.
#'
#'          When the model type is \code{"regression"}, the value of the irrelevant items is set to \code{NULL}. Likewise, when the model type is \code{"classification"},
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
rbioClass_svm_predict <- function(object,
                                  newdata, center.scale.newdata = TRUE, export.name = NULL,
                                  sampleID.vector = NULL, newdata.y = NULL,
                                  prob.method = c("logistic",  "Bayes", "softmax"),
                                  verbose = TRUE){
  ## argument check
  if (!any(class(object) %in% c("rbiosvm"))) stop("object needs to be a \"rbiosvm\" class.")
  if (missing(newdata) || is.null(newdata)) stop("please provide newdata.")
  if (is.null(export.name)){
    export.name <- deparse(substitute(newdata))
  } else {
    export.name <- export.name
  }
  if (!any(class(newdata) %in% c("matrix", "data.frame"))) stop("newdata has to be either a matrix or data.frame object.")
  if (ncol(newdata) != ncol(object$inputX)) stop("newdata needs to have the same number of variables, i.e. columns, as the object.")
  # if (!prob.method %in% c("logistic",  "Bayes", "softmax")) stop("Probability method should be either \"softmax\" or \"Bayes\".")
  if (center.scale.newdata){
    if (is.null(object$center.scaledX)) stop("No center.scaledX found in training data while center.scale.newdata = TRUE.")
  }
  if (!missing(sampleID.vector) & !is.null(sampleID.vector)){
    if (length(sampleID.vector) != nrow(newdata)){
      cat("Sample ID vector not the same length as newdata. Proceed without custom sample labels.\n")
      sampleID.vector <- NULL
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
  if (any(class(newdata) == "data.frame")){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    # newdata <- as.matrix(sapply(newdata, as.numeric)) # testing
    newdata <- as.matrix(newdata)  # testing
  }
  if (center.scale.newdata){
    if (verbose) cat("Data center.scaled using training data column mean and sd, prior to predicting.")
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
    if (!missing(sampleID.vector) & !is.null(sampleID.vector)){
      rownames(pred_mtx) <- sampleID.vector
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
      bayes.prob$train.posterior <- predict(bayes.prob)$posterior  # calculate posterior probability for
      bayes.prob$x <- NULL

      bayespred <- predict(object = bayes.prob, newdata = pred_mtx)
      prob <- bayespred$posterior
    } else {
      group <- colnames(pred_mtx)
      # calculate probability
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
              newdata.id = sampleID.vector,
              center.scale = center.scale.newdata,
              center.scaled.newdata = centerdata,
              newdata.y = y)
  class(out) <- "prediction"
  assign(paste0(export.name, "_svm_predict"), out, envir = .GlobalEnv)
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


