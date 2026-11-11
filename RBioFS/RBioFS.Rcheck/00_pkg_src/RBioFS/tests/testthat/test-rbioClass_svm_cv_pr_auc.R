# Tests for rbioClass_svm_cv_pr_auc()
#
# Exercises the function against real e1071mc SVM models built in
# helper-cv_pr_auc.R. Two-class and multi-class CV scenarios, nested vs.
# single CV, regression rejection, simpleError fold removal, and the
# absen-class edge case (a CV training fold that never observed a class).

# ----- fixture: 2-class stratified CV (3 folds) -----
make_two_class_cv <- function() {
  d <- make_two_class_data(seed = 1)
  X <- d$X; y <- d$y
  folds <- vector("list", 3)
  # stratified split: each fold's test set contains both classes
  tr1 <- c(rep(1:30, 1)[1:10], rep(31:60, 1)[1:10])  # 10 A + 10 B
  tr1 <- c(1:10, 31:40)
  te1 <- c(11:20, 41:50)
  tr2 <- c(11:20, 41:50)
  te2 <- c(1:10, 31:40)
  tr3 <- c(1:10, 11:20, 31:40, 41:50)
  te3 <- c(21:30, 51:60)
  folds[[1]] <- build_cv_fold(build_svm_model(X[tr1, ], y[tr1]), X[te1, ], y[te1])
  folds[[2]] <- build_cv_fold(build_svm_model(X[tr2, ], y[tr2]), X[te2, ], y[te2])
  folds[[3]] <- build_cv_fold(build_svm_model(X[tr3, ], y[tr3]), X[te3, ], y[te3])
  names(folds) <- paste0("cv_fold_", 1:3)
  list(obj = build_cv_object(folds), folds = folds)
}

# ----- fixture: 3-class stratified CV (2 folds) -----
make_three_class_cv <- function() {
  d <- make_three_class_data(seed = 2)
  X <- d$X; y <- d$y
  tr  <- c(1:3, 16:18, 31:33)
  te  <- c(4:15, 19:30, 34:45)
  folds <- vector("list", 2)
  folds[[1]] <- build_cv_fold(build_svm_model(X[tr, ], y[tr]),  X[te, ], y[te])
  folds[[2]] <- build_cv_fold(build_svm_model(X[te, ], y[te]),  X[tr, ], y[tr])
  names(folds) <- paste0("cv_fold_", 1:2)
  list(obj = build_cv_object(folds), folds = folds, X = X, y = y)
}

# ----- 1. 2-class: basic structure + AUC range + F1 sanity -----
test_that("2-class cv: returns a per-fold list with expected keys and valid AUCs", {
  fx <- make_two_class_cv()
  res <- rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  expect_length(res, 3)
  expect_named(res, paste0("cv_fold_", 1:3))
  expect_setequal(names(res[[1]]), c("svm.pr_object", "svm.pr_dataframe",
                                     "svm.f1", "svm.pr_auc"))
  # per-class AUCs in [0, 1]
  auc <- res[[1]]$svm.pr_auc
  expect_named(auc, c("A", "B"))
  expect_true(all(auc >= 0 & auc <= 1))
  # F1 structure
  expect_setequal(names(res[[1]]$svm.f1), c("per.class", "macro", "micro"))
  expect_named(res[[1]]$svm.f1$per.class, c("class", "F1"))
  expect_true(all(res[[1]]$svm.f1$per.class$F1 >= 0 &
                  res[[1]]$svm.f1$per.class$F1 <= 1))
  expect_true(res[[1]]$svm.f1$macro$F1 >= 0 && res[[1]]$svm.f1$macro$F1 <= 1)
  expect_true(res[[1]]$svm.f1$micro$F1 >= 0 && res[[1]]$svm.f1$micro$F1 <= 1)
})

# ----- 2. 2-class: pr_dataframe columns + group labels -----
test_that("2-class cv: pr_dataframe has correct columns and per-class groups", {
  fx <- make_two_class_cv()
  res <- rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  df <- res[[1]]$svm.pr_dataframe
  expect_setequal(names(df), c("precision", "recall", "threshold", "group"))
  expect_setequal(unique(df$group), c("A", "B"))
  # precision & recall in [0,1]
  expect_true(all(df$precision >= 0 & df$precision <= 1))
  expect_true(all(df$recall >= 0 & df$recall <= 1))
})

# ----- 3. 3-class: vs Others group labels + 3-class AUCs/F1 -----
test_that("3-class cv: pr_dataframe uses 'vs Others' group labels and has 3 classes", {
  fx <- make_three_class_cv()
  res <- rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  expect_length(res, 2)
  df <- res[[1]]$svm.pr_dataframe
  expect_setequal(unique(df$group), c("A (vs Others)", "B (vs Others)", "C (vs Others)"))
  expect_named(res[[1]]$svm.pr_auc, c("A", "B", "C"))
  expect_true(all(res[[1]]$svm.pr_auc >= 0 & res[[1]]$svm.pr_auc <= 1))
  # 3-class-specific F1 has 3 rows
  expect_true(nrow(res[[1]]$svm.f1$per.class) == 3)
  expect_setequal(res[[1]]$svm.f1$per.class$class, c("A", "B", "C"))
})

# ----- 4. nested CV object -----
test_that("nested cv: rbiosvm_nestedcv is accepted and returns per-fold list", {
  fx <- make_three_class_cv()
  nest <- build_cv_object(fx$folds, nested = TRUE)
  expect_s3_class(nest, "rbiosvm_nestedcv")
  res <- rbioClass_svm_cv_pr_auc(object = nest, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  expect_length(res, 2)
  expect_setequal(names(res[[1]]), c("svm.pr_object", "svm.pr_dataframe",
                                     "svm.f1", "svm.pr_auc"))
})

# ----- 5. regression model rejected -----
test_that("regression model is rejected", {
  fx <- make_two_class_cv()
  reg <- build_cv_object(fx$folds, model_type = "regression")
  expect_error(
    rbioClass_svm_cv_pr_auc(object = reg, fileprefix = withr::local_tempfile(),
                            prplot = FALSE, verbose = FALSE),
    "classification"
  )
})

# ----- 6. wrong class rejected -----
test_that("non-rbiosvm_cv / rbiosvm_nestedcv object is rejected", {
  bad <- structure(list(cv.models = list(), model.type = "classification"),
                   class = "not_a_valid_class")
  expect_error(
    rbioClass_svm_cv_pr_auc(object = bad, fileprefix = withr::local_tempfile(),
                            prplot = FALSE, verbose = FALSE),
    "classes"
  )
})

# ----- 7. simpleError fold removed -----
test_that("fold returned as simpleError is removed from the result", {
  fx <- make_three_class_cv()
  folds <- fx$folds
  folds[[2]] <- structure(list(), class = "simpleError")
  obj <- build_cv_object(folds)
  res <- rbioClass_svm_cv_pr_auc(object = obj, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  expect_length(res, 1)
  expect_named(res, "cv_fold_1")
})

# ----- 8. absent-class edge case (fixed bug) -----
# A CV fold whose training set never saw class A produces an SVM model with
# probabilities only for B and C; the test fold's single class (A) has no
# matching probability column. Previously this crashed with "subscript out of
# bounds"; the fixed function skips the class gracefully and drops the fold.
test_that("absent-class fold is skipped gracefully (no crash)", {
  d <- make_three_class_data(seed = 2)
  X <- d$X; y <- d$y
  tr  <- 16:45   # training set: B and C only (no A)
  te  <- 1:15    # test set: all class A
  model <- build_svm_model(X[tr, ], y[tr])
  # the model's probability output only has B, C columns
  pred <- predict(model, newdata = X[te, ], decision.values = TRUE, probability = TRUE)
  prob <- attr(pred, "probabilities")
  expect_true(!("A" %in% colnames(prob)))

  folds <- vector("list", 1)
  folds[[1]] <- build_cv_fold(model, X[te, ], y[te])
  names(folds) <- "cv_fold_1"
  obj <- build_cv_object(folds)
  # should not throw
  res <- rbioClass_svm_cv_pr_auc(object = obj, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  # the single fold had no analyzable class -> removed -> empty result
  expect_length(res, 0)
})

# ----- 9. prplot = TRUE writes PDFs per class (2-class) -----
test_that("prplot writes one PDF per class", {
  fx <- make_two_class_cv()
  withr::with_tempdir({
    res <- rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = "run",
                                   prplot = TRUE, verbose = FALSE)
    pdfs <- list.files(pattern = "^run\\.cv_pr\\..*\\.pdf$")
    expect_setequal(pdfs, c("run.cv_pr.A.pdf", "run.cv_pr.B.pdf"))
  })
})

# ----- 10. prplot = FALSE writes no PDFs -----
test_that("prplot = FALSE produces no PDF output", {
  fx <- make_two_class_cv()
  withr::with_tempdir({
    res <- rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = "run",
                                   prplot = FALSE, verbose = FALSE)
    pdfs <- list.files(pattern = "^run\\.cv_pr\\..*\\.pdf$")
    expect_length(pdfs, 0)
  })
})

# ----- 11. F1-micro equals overall accuracy for this data set -----
test_that("F1-micro equals overall accuracy (diag/sum of confusion matrix)", {
  fx <- make_two_class_cv()
  res <- rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = withr::local_tempfile(),
                                 prplot = FALSE, verbose = FALSE)
  # recompute accuracy from the fold's cv_test_data and the model
  fold <- fx$folds[[1]]
  cv_test <- fold$cv_test_data
  cv_test_x <- cv_test[, !names(cv_test) %in% "y"]
  pred_lab <- predict(fold$cv_svm_model, newdata = cv_test_x)
  acc <- sum(as.character(pred_lab) == as.character(cv_test$y)) / length(cv_test$y)
  expect_equal(res[[1]]$svm.f1$micro$F1, acc, tolerance = 1e-9)
})

# ----- 12. result exported to .GlobalEnv -----
test_that("result is exported to the global environment", {
  fx <- make_two_class_cv()
  # When fileprefix is supplied the result is stored in .GlobalEnv under
  # '<fileprefix>_svm_cv_pr_auc' (for rbiosvm_cv) respectively
  # '<fileprefix>_svm_nestedcv_pr_auc' (for rbiosvm_nestedcv).
  rbioClass_svm_cv_pr_auc(object = fx$obj, fileprefix = "exp_prefix",
                            prplot = FALSE, verbose = FALSE)
  expect_true(exists("exp_prefix_svm_cv_pr_auc", envir = .GlobalEnv))
  on.exit(rm(list = "exp_prefix_svm_cv_pr_auc", envir = .GlobalEnv))
})
