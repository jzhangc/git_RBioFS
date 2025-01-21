#' @title uni_PreProc
#'
#' @description a legacy function borrowed from RBioArray. Data pre-processing function for the microarary data, which is now used for univariate reduction during svm nested CV.
#' @param rawlist Input data, either a list, \code{EList} or \code{MAList} object.
#' @param logTrans If to perform a log transformation on the data or not. Default is \code{FALSE}.
#' @param logTransMethod If \code{logTrans = TRUE}, set which method to use for the transformation, \code{"log2"} or \code{"log10"}. Default is \code{"log2"}.
#' @param logTransObjT If \code{logTrans = TRUE}, set the file name for the output \code{csv} file containing the log transformed data.
#' @param logTransParallelComputing If \code{logTrans = TRUE}, set if to use parallel computing for the transformation or not. Default is \code{FALSE}.
#' @param bgMethod Background correction method. Default is \code{"auto"}. See \code{backgroundCorrect()} function from \code{limma} package for details.
#' @param normMethod Normalization method. Default is \code{"quantile"}. See \code{normalizeBetweenArrays()} function from \code{limma} package for details.
#' @param ... arguments for \code{backgroundCorrect.matrix()} or \code{backgroundCorrect()} functions from \code{limma} package.
#' @details The function does not use design matrix for array weight calculation.
#'          Therefore, the DE analysis based on the output from this function will yield slightly different results from the \code{\link{rbioarray_filter_combine}}.
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


#' @title rbioClass_svm_ncv_fs_legacy
#'
#' @description This is a legacy version of the function, which has nested cv svm modelling only running on a single CPU thread.
#'              Nested cross-validation feature selection and SVM modelling: CV-rRF-FS-SVM.
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
#' For now, RBioFS implementation of two-step random forest feature selection is used to select features based on nested cross-validation.
#' Resulted features from each nested cross-validation round are voted. Features with votes equal or greater than the cutoff are reported as selected features.
#'
#' It is also a good idea to set \code{fs.count.cutoff} as \code{cross.k - 1}. Notably, \code{fs.count.cutoff = 1} is "no threshold",
#' meaning maximum number of features selected from the nested cross-validation are reported.
#'
#' When \code{cross.best.model.method = "median"}, the function only use models with accuracy/RMSE equal or better than the median value
#' for feature count threholding. When there is no change in performance across cv models, the function behaves same as \code{cross.best.model.method = "none"}
#'
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
#' svm_nestedcv <- rbioClass_svm_ncv_fs_legacy(x = mydata[, -c(1:2)],
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
rbioClass_svm_ncv_fs_legacy <- function(x, y,
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
    if (verbose) cat(paste0("Nested CV iteration: ", i, "|", cross.k, "..."))
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

      # update the cv training data
      # cv_training <- cv_training[, c("y", uni_sig_fs), drop = FALSE]
      # cv_training_x <- cv_training[, -1, drop = FALSE]
      # cv_training_y <- cv_training[, 1]
      if (length(uni_sig_fs) == 0) {
        stop("No statistically significant features found. Try runing with larger uni.alpha value, or univarate.fs = FALSE.")
        # cv_training_x <- cv_training[, -1]
        # cv_training_y <- cv_training[, 1]
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

      if (length(eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features) > 1){
        out <- eval(parse(text = paste0("svm_nested_iter_", i, "_SFS")))$selected_features
      } else {
        out <- eval(parse(text = paste0("svm_nested_iter_", i, "_initial_FS")))$feature_initial_FS
      }
    },
    # warning = function(w){
    #   cat(paste0("Warnings occurred during rRF-FS, skipping rRF-FS for CV iteration: ", i, "..."))
    #   out <- colnames(fs_training_x)
    #   return(out)
    # },
    error = function(e){
      cat(paste0("Errors occurred during rRF-FS, skipping rRF-FS for CV iteration: ", i, "...\n", "\tError message: ", e))
      out <- colnames(fs_training_x)
      return(out)
    })

    cv_m <- rbioClass_svm(x = fs_training_x[, fs], y = cv_training_y, center.scale = center.scale,
                          svm.cross.k = 0, tune.method = tune.method, kernel = kernel,
                          tune.cross.k = tune.cross.k, tune.boot.n = tune.boot.n, verbose = FALSE, ...)
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
  nested.cv.list[] <- foreach(i = 1:cross.k, .packages = c("foreach", "RBioFS"), .errorhandling = "stop") %do% nestedcv_func(i)
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
