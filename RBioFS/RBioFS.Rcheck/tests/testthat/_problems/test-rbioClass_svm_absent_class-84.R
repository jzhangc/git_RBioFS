# Extracted from test-rbioClass_svm_absent_class.R:84

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RBioFS", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
make_absent_class_fixture <- function(seed = 2) {
  d <- make_three_class_data(seed = seed)
  X <- d$X; y <- d$y
  tr <- 16L:45L   # training set: B and C only (no A)
  te <- 1L:15L    # test set: all class A
  model <- build_svm_model(X[tr, ], y[tr])
  # confirm the model's probability output has no "A" column
  pred <- predict(model, newdata = X[te, ],
                  decision.values = TRUE, probability = TRUE)
  prob <- attr(pred, "probabilities")
  list(X = X, y = y, tr = tr, te = te, model = model, prob = prob)
}

# test -------------------------------------------------------------------------
fx <- make_absent_class_fixture()
folds <- vector("list", 1)
folds[[1]] <- build_cv_fold(fx$model, fx$X[fx$te, ], fx$y[fx$te])
names(folds) <- "cv_fold_1"
obj <- build_cv_object(folds)
res <- rbioClass_svm_cv_roc_auc(object = obj,
                                  fileprefix = withr::local_tempfile(),
                                  rocplot = FALSE, plot.smooth = FALSE,
                                  verbose = FALSE)
expect_null(res)
