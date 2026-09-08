# Fixture builders for rbioClass_svm_cv_pr_auc tests.
#
# These builders construct minimal but realistic `rbiosvm` / `rbiosvm_cv` /
# `rbiosvm_nestedcv` objects using a real e1071mc SVM model, mirroring the
# structure produced by rbioClass_svm() and its CV workflows.

# Build a real e1071mc SVM model with probability output enabled.
build_svm_model <- function(X, y, seed = 42) {
  if (!is.null(seed)) set.seed(seed)
  m <- e1071mc::svm(x = X, y = y, kernel = "radial", probability = TRUE,
                     cost = 1, gamma = 1 / ncol(X), cross = 0)
  m$inputX <- X
  m$inputY <- y
  m$center.scaledX <- NULL
  m$model.type <- "classification"
  class(m) <- c("rbiosvm", "svm")
  m
}

# Build a single CV fold: an `rbiosvm` model plus its held-out `cv_test_data`
# (a data.frame with a `y` factor column and feature columns, no scaling info).
build_cv_fold <- function(model, Xtest, ytest) {
  td <- data.frame(Xtest, y = ytest, check.names = FALSE)
  list(cv_svm_model = model,
       cv.accuracy = NA,
       cv_test_data = td)
}

# Assemble a list of folds into an `rbiosvm_cv` or `rbiosvm_nestedcv` object.
build_cv_object <- function(folds, model_type = "classification",
                           nested = FALSE) {
  if (nested) {
    structure(list(nested.cv.models = folds, model.type = model_type),
              class = "rbiosvm_nestedcv")
  } else {
    structure(list(cv.models = folds, model.type = model_type),
              class = "rbiosvm_cv")
  }
}

# 2-class design: 30 + 30 samples, a margin on feature 1 so the SVM separates.
make_two_class_data <- function(seed = 1) {
  set.seed(seed)
  n <- 60; p <- 6
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", 1:p)))
  X[1:30, 1] <- X[1:30, 1] + 2
  X[31:60, 1] <- X[31:60, 1] - 2
  y <- factor(rep(c("A", "B"), each = 30))
  list(X = X, y = y)
}

# 3-class design: 15 + 15 + 15 samples, a margin on a distinct feature each.
make_three_class_data <- function(seed = 2) {
  set.seed(seed)
  n <- 45; p <- 5
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("v", 1:p)))
  X[1:15, 1] <- X[1:15, 1] + 2
  X[16:30, 2] <- X[16:30, 2] + 2
  X[31:45, 3] <- X[31:45, 3] + 2
  y <- factor(rep(c("A", "B", "C"), each = 15))
  list(X = X, y = y)
}
