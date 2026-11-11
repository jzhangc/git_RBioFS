# Regression tests for the "absent-class CV fold" crash.
#
# A CV training fold that never observed a class produces an SVM model with no
# probability column for that class. The per-class loops then did
#   pred_prob[, levels(response)[i]]
# which throws "subscript out of bounds" for the missing column. The same latent
# bug existed in five sites that were fixed alongside rbioClass_svm_cv_pr_auc:
#   rbioClass_svm_roc_auc          (rbioClass_svm_assessment.R)
#   rbioClass_svm_roc_auc_inter    (single-model path + svm_cv_rocauc_helper)
#   rbioClass_svm_cv_roc_auc       (cv path)
#   svm_cv_rocauc_helper           (reached via rbioClass_svm_cv_roc_auc_v2)
#   rbioClass_svm_pr_auc
# Each test below drives the corresponding function with a model trained on B+C
# only and a test set of class A, so the "A" probability column is absent, and
# verifies the function returns gracefully instead of crashing.

# Build an SVM model trained without class A (only B and C).
# The test data is class A only -> "A" has no column in the probability output.
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

# ----- 1. rbioClass_svm_roc_auc (site 173) -----
test_that("rbioClass_svm_roc_auc handles absent-class newdata without crashing", {
  fx <- make_absent_class_fixture()
  expect_true(!("A" %in% colnames(fx$prob)))
  res <- rbioClass_svm_roc_auc(object = fx$model,
                                newdata = fx$X[fx$te, ],
                                newdata.label = fx$y[fx$te],
                                center.scale.newdata = FALSE,
                                fileprefix = withr::local_tempfile(),
                                rocplot = FALSE, verbose = FALSE)
  # no analyzable class -> NULL (no crash)
  expect_null(res)
})

# ----- 2. rbioClass_svm_roc_auc_inter single-model path (site 423) -----
test_that("rbioClass_svm_roc_auc_inter handles absent-class newdata without crashing", {
  fx <- make_absent_class_fixture()
  expect_true(!("A" %in% colnames(fx$prob)))
  res <- rbioClass_svm_roc_auc_inter(object = fx$model,
                                     newdata = fx$X[fx$te, ],
                                     newdata.label = fx$y[fx$te],
                                     center.scale.newdata = FALSE,
                                     fileprefix = withr::local_tempfile(),
                                     rocplot = FALSE, verbose = FALSE)
  expect_null(res)
})

# ----- 3. svm_cv_rocauc_helper via rbioClass_svm_cv_roc_auc_v2 (site 670) -----
test_that("rbioClass_svm_cv_roc_auc_v2 handles absent-class fold without crashing", {
  fx <- make_absent_class_fixture()
  folds <- vector("list", 1)
  folds[[1]] <- build_cv_fold(fx$model, fx$X[fx$te, ], fx$y[fx$te])
  names(folds) <- "cv_fold_1"
  obj <- build_cv_object(folds)
  res <- rbioClass_svm_cv_roc_auc_v2(object = obj,
                                     fileprefix = withr::local_tempfile(),
                                     rocplot = FALSE, verbose = FALSE)
  # the single fold had no analyzable class -> NULL (no crash)
  expect_null(res)
})

# ----- 4. rbioClass_svm_cv_roc_auc cv path (site 852) -----
test_that("rbioClass_svm_cv_roc_auc handles absent-class fold without crashing", {
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
})

# ----- 5. rbioClass_svm_pr_auc (site 2153) -----
test_that("rbioClass_svm_pr_auc handles absent-class newdata without crashing", {
  fx <- make_absent_class_fixture()
  expect_true(!("A" %in% colnames(fx$prob)))
  res <- rbioClass_svm_pr_auc(object = fx$model,
                              newdata = fx$X[fx$te, ],
                              newdata.label = fx$y[fx$te],
                              center.scale.newdata = FALSE,
                              fileprefix = withr::local_tempfile(),
                              prplot = FALSE, verbose = FALSE)
  # returns a svm_pr_auc object; the absent class entry is NULL (no crash)
  expect_s3_class(res, "svm_pr_auc")
  expect_null(res$svm.pr_object["A"][[1]])
  # F1 (argmax-based) is still computed
  expect_setequal(names(res$svm.f1), c("per.class", "macro", "micro"))
})
