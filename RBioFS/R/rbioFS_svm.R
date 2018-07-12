#' @title rbioFS_svm
#'
#' @description Support Vector Machine (SVM) modelling
#' @param x Input data matrix (e.g., independent variables, predictors, features, X, etc). Make sure it is either a matrix or a dataframe.
#' @param y Input response variable (e.g.,dependent variables, Y etc). Make sure it is \code{factor} class.
#' @param center Logical, wether center the data, i.e. subtracting mean. Default is \code{TRUE}.
#' @param scale Logical, whether to scale the data, i.e. dividing by standard deviation. Default is \code{TRUE}.
#' @param kernel SVM kernel. Options are \code{"linear", "ploynomial", "radial", "sigmoid"}. Default is \code{"radial"}, aka RBF.
#' @param tune.method Parameter tuning method. Options are \code{"cross"} (i.e. cross validation), \code{"boot"} (i.e. bootstrap), and \code{"fix"}. Default is \code{"cross"}.
#' @param tune.cross.k Set only when \code{tune.method = "cross"}, fold number for cross validation. Default is \code{10}.
#' @param tune.boot.n Set only when \code{tune.method = "boot"}, bootstrap iterations. Default is \code{10}.
#' @param ... Additional arguments for \code{svm} function from \code{e1071} pacakge.
#' @param verbose Wether to display messages. Default is \code{TRUE}. This will be affect error or warning messeages.
#' @return Returns a SVM model object, with classes "svm" and "rbiosvm".
#' @details Parameter tuning is for gamma (not applicable when \code{kernel = "linear"}) and cost.
#' @importFrom e1071 svm tune tune.control
#' @examples
#' \dontrun{
#' svm_model <- rbioFS_svm(x = training_set[, -1], y = training_set[, 1], kernel = "radial", center = TRUE, scale = FALSE)
#' }
#' @export
rbioFS_svm <- function(x, y, center = TRUE, scale = TRUE,
                       kernel = "radial",
                       tune.method = "cross", tune.cross.k = 10, tune.boot.n = 10, ...,
                       verbose = TRUE){
  ## check arguments
  if (class(x) == "data.frame"){
    if (verbose) cat("data.frame x converted to a matrix object.\n")
    x <- as.matrix(x)
  }

  ## data processing
  if (center){
    if (verbose) cat(paste0("Data centered with the scale option ", ifelse(scale, "\"ON\" ", "\"OFF\" "), "prior to modelling..."))
    centered_X <- center_scale(x, scale = scale)  # center data with the option of scaling
    if (verbose) cat("DONE!\n")
    scale_X <- "See centerX"
    X <- centered_X$centerX
  } else {
    if (scale){
      if (verbose) cat(paste0("Data scaling prior to modelling..."))
      col.mean <- colMeans(x, na.rm = TRUE)
      col.sd <- matrixStats::colSds(x, center = col.mean, na.rm = TRUE) # matrixStats::colSds
      scale_X <- t(t(x) / col.sd)  # scale without centering
      X <- scale_X
      if (verbose) cat("DONE!\n")
    } else {
      scale_X <- NULL
      X <- x
    }
    centered_X <- NULL
  }
  y <- y

  ## weight evaluation

  ## svm
  # tune parameters
  gamma_start <- ifelse(is.vector(X), 1, 1 / ncol(X))
  if (verbose) cat(paste0("Parameter tuning for cost and gamma with ", tune.method, " method (speed depending on hardware configuration)..."))
  svm_tuned <- tune(svm, kernel = kernel,
                    train.x = X, train.y = y, ranges = list(gamma = gamma_start * 2^(-5:5), cost = 2^(-10:8)),
                    tunecontrol = tune.control(sampling = tune.method, cross = tune.cross.k, nboot = tune.boot.n))
  if (verbose) cat("DONE!\n")

  # svm
  m <- svm(x = X, y = y, kernel = kernel,cost = svm_tuned$best.parameters$cost, gamma = svm_tuned$best.parameters$gamma,
           scale = FALSE,...)

  # return
  m$inputX <- x
  m$inputY <- y
  m$scaleX <- scale_X
  m$centerX <- centered_X
  class(m) <- c("svm", "rbiosvm")
  return(m)
}
