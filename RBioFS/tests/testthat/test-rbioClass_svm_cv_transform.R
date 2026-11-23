# Regression tests for the non-leakage CV pipeline fixes:
#
#  (T1) the CV test fold must be transformed with the parameters fitted on
#       that fold's TRAINING rows. Previously the training data was transformed
#       once in the FS/CV-prep step and then transformed AGAIN inside
#       rbioClass_svm() (center.scale = TRUE passthrough); the stored transform
#       parameters were therefore near-identity (mean ~ 0, SD ~ 1) and the test
#       fold was left effectively UNtransformed while the model was trained on
#       transformed data (train/test scale mismatch). Fixed by fitting the
#       transform once, passing center.scale = FALSE to the inner model, and
#       re-applying the stored transform object to the test fold.
#  (T2) best.rmse.summary / best.nested.rmse.summary must hold the RMSE summary,
#       not the accuracy summary (previously NULL for regression).
#  (T3) sampleIds in rbioClass_svm_cv / _ncv_fs protects only the INNER CV;
#       a warning must be emitted (row-level outer folds can split groups),
#       and a length mismatch must error early.
#  (T4) uni.sig.fs in the nested CV output must be the per-fold univariate
#       feature list (previously always NULL, read from the wrong list element).
#  (T5) v3 group-aware OUTER folds: rows sharing a sampleId stay in one fold.
#  (T6) performance sanity: zcore_minmax is not slower than base::scale
#       (it is built on compiled matrixStats + sweep).

# ---- shared fixture: features with large offset / scale so that the
# ---- transform actually changes the data (a no-op transform would pass T1
# ---- trivially) ----
make_cv_fixture <- function(seed = 11, n = 60, p = 5) {
  set.seed(seed)
  x <- cbind(
    f_offset = rnorm(n, mean = 50, sd = 5),
    f_wide   = rnorm(n, mean = 0, sd = 100),
    f_shift  = rnorm(n, mean = -20, sd = 1),
    f_plain  = rnorm(n),
    f_plain2 = rnorm(n))
  if (p < 5) x <- x[, seq_len(p), drop = FALSE]
  y <- factor(c(rep("A", n %/% 2), rep("B", n - n %/% 2)))
  list(x = x, y = y)
}

# T1a: rbioClass_svm_cv, z-score path
test_that("T1a: rbioClass_svm_cv test fold uses the TRAIN-fold transform (z-score)", {
  d <- make_cv_fixture()
  set.seed(2026)
  obj <- rbioClass_svm_cv(x = d$x, y = d$y, cross.k = 3,
                          tune.method = "cross", tune.cross.k = 3,
                          parallelComputing = FALSE, verbose = FALSE)
  fold <- obj$cv.fold  # original row order (classification: no shuffle)
  for (i in 1:3) {
    tr <- which(fold != i); te <- which(fold == i)
    params <- zcore_minmax(d$x[tr, ], zscore = TRUE, minmax = FALSE, verbose = FALSE)
    expected <- RBioFS:::.apply_zcore_minmax(d$x[te, ], params = params, verbose = FALSE)
    observed <- as.matrix(obj$cv.models[[i]]$cv_test_data[, -1, drop = FALSE])
    expect_equal(as.numeric(expected), as.numeric(observed), tolerance = 1e-10,
                 info = paste("fold", i, "test transform mismatch"))
    # the stored test data must NOT be the raw test rows (the old bug)
    expect_false(isTRUE(all.equal(expected, d$x[te, ], tolerance = 1e-6)))
  }
})

# T1b: rbioClass_svm_cv, min-max path
test_that("T1b: rbioClass_svm_cv test fold uses the TRAIN-fold transform (min-max)", {
  d <- make_cv_fixture(seed = 12)
  set.seed(2026)
  obj <- rbioClass_svm_cv(x = d$x, y = d$y, cross.k = 3,
                          zscore = FALSE, minmax = TRUE,
                          tune.method = "cross", tune.cross.k = 3,
                          parallelComputing = FALSE, verbose = FALSE)
  fold <- obj$cv.fold
  for (i in 1:3) {
    tr <- which(fold != i); te <- which(fold == i)
    params <- zcore_minmax(d$x[tr, ], zscore = FALSE, minmax = TRUE, verbose = FALSE)
    expected <- RBioFS:::.apply_zcore_minmax(d$x[te, ], params = params, verbose = FALSE)
    observed <- as.matrix(obj$cv.models[[i]]$cv_test_data[, -1, drop = FALSE])
    expect_equal(as.numeric(expected), as.numeric(observed), tolerance = 1e-10)
  }
})

# T1c: rbioClass_svm_ncv_fs_v2 (rRF-FS in the loop) test-fold transform
test_that("T1c: ncv_fs_v2 test fold uses the TRAIN-fold transform", {
  d <- make_cv_fixture(seed = 13, n = 48)
  set.seed(2026)
  obj <- rbioClass_svm_ncv_fs_v2(x = d$x, y = d$y, cross.k = 2,
                                 rf.ifs.ntree = 50, rf.sfs.ntree = 50,
                                 tune.method = "cross", tune.cross.k = 3,
                                 parallelComputing = FALSE, verbose = FALSE)
  fold <- obj$cv.fold
  for (i in 1:2) {
    tr <- which(fold != i); te <- which(fold == i)
    fs <- obj$nested.cv.models[[i]]$selected.features
    # the transform is column-wise, so fitting on the selected columns of the
    # training rows gives identical parameters to subsetting the full fit
    params <- zcore_minmax(d$x[tr, fs, drop = FALSE], zscore = TRUE, minmax = FALSE, verbose = FALSE)
    expected <- RBioFS:::.apply_zcore_minmax(d$x[te, fs, drop = FALSE], params = params, verbose = FALSE)
    observed <- as.matrix(obj$nested.cv.models[[i]]$cv_test_data[, fs, drop = FALSE])
    expect_equal(as.numeric(expected), as.numeric(observed), tolerance = 1e-10)
  }
})

# T1d: ncv_fs (v1, two-stage list pipeline) test-fold transform
test_that("T1d: ncv_fs (v1) test fold uses the TRAIN-fold transform", {
  d <- make_cv_fixture(seed = 14, n = 48)
  set.seed(2026)
  obj <- rbioClass_svm_ncv_fs(x = d$x, y = d$y, cross.k = 2,
                              rf.ifs.ntree = 50, rf.sfs.ntree = 50,
                              tune.method = "cross", tune.cross.k = 3,
                              parallelComputing = FALSE, verbose = FALSE)
  fold <- obj$cv.fold
  for (i in 1:2) {
    tr <- which(fold != i); te <- which(fold == i)
    fs <- obj$nested.cv.models[[i]]$selected.features
    # the transform is column-wise, so fitting on the selected columns of the
    # training rows gives identical parameters to subsetting the full fit
    params <- zcore_minmax(d$x[tr, fs, drop = FALSE], zscore = TRUE, minmax = FALSE, verbose = FALSE)
    expected <- RBioFS:::.apply_zcore_minmax(d$x[te, fs, drop = FALSE], params = params, verbose = FALSE)
    observed <- as.matrix(obj$nested.cv.models[[i]]$cv_test_data[, fs, drop = FALSE])
    expect_equal(as.numeric(expected), as.numeric(observed), tolerance = 1e-10)
  }
})

# T2: regression summaries must not be NULL and must be the RMSE summary
test_that("T2: rbioClass_svm_cv regression best.rmse.summary is the RMSE summary", {
  set.seed(2026)
  x <- matrix(rnorm(60 * 4), 60, 4)
  y <- rnorm(60)
  obj <- rbioClass_svm_cv(x = x, y = y, cross.k = 3,
                          tune.method = "cross", tune.cross.k = 3,
                          parallelComputing = FALSE, verbose = FALSE)
  expect_null(obj$best.accuracy.summary)
  expect_false(is.null(obj$best.rmse.summary))
  best.idx <- obj$best.index
  expect_equal(unname(obj$best.rmse.summary),
               c(mean(obj$best.rmse), sd(obj$best.rmse),
                 sd(obj$best.rmse) / sqrt(length(best.idx))),
               tolerance = 1e-12)
  # the best RMSEs are a subset of the per-fold RMSEs (median rule)
  expect_true(all(obj$best.rmse <= median(obj$RMSE) + 1e-12))
})

# T3a: sampleIds in rbioClass_svm_cv warns about row-level outer folds
test_that("T3a: rbioClass_svm_cv with sampleIds warns that outer folds are row-level", {
  d <- make_cv_fixture(seed = 15, n = 40)
  sids <- as.character(rep(1:20, each = 2))
  expect_warning(
    obj <- rbioClass_svm_cv(x = d$x, y = d$y, cross.k = 2, sampleIds = sids,
                            tune.method = "cross", tune.cross.k = 3,
                            parallelComputing = FALSE, verbose = FALSE),
    "row-level")
  # inner CV still runs and records the ids (subset to the 20 training rows)
  expect_equal(length(obj$cv.models[[1]]$cv_svm_model$sampleIds), 20)
})

# T3b: sampleIds length mismatch errors early
test_that("T3b: rbioClass_svm_cv rejects mis-sized sampleIds", {
  d <- make_cv_fixture(seed = 16, n = 40)
  expect_error(
    rbioClass_svm_cv(x = d$x, y = d$y, cross.k = 2, sampleIds = 1:5,
                     tune.method = "cross", tune.cross.k = 3,
                     parallelComputing = FALSE, verbose = FALSE),
    "same length")
})

# T3c: sampleIds with conflicting class labels errors early (ncv_fs v1)
test_that("T3c: ncv_fs rejects same-id rows with different labels", {
  d <- make_cv_fixture(seed = 17, n = 40)
  sids <- as.character(rep(1:20, each = 2))
  y <- d$y
  y[2] <- "B"  # row 2 shares id "1" with row 1 (class A)
  expect_error(
    rbioClass_svm_ncv_fs(x = d$x, y = y, cross.k = 2, sampleIds = sids,
                         rf.ifs.ntree = 50, rf.sfs.ntree = 50,
                         tune.method = "cross", tune.cross.k = 3,
                         parallelComputing = FALSE, verbose = FALSE),
    "differing class labels")
})

# T4: uni.sig.fs is the per-fold univariate list (not NULL)
test_that("T4: ncv_fs univariate.fs stores per-fold uni.sig.fs", {
  # two features strongly associated with the class so limma finds them
  set.seed(2026)
  n <- 60
  y <- factor(c(rep("A", 30), rep("B", 30)))
  x <- cbind(sig1 = rnorm(n) + as.numeric(y == "B") * 4,
             sig2 = rnorm(n) - as.numeric(y == "B") * 4,
             noise = rnorm(n))
  # uni.log2trans = FALSE: the fixture's sig2 column takes negative values,
  # which are invalid for the default log2(x + 2) transform (produces NaNs)
  obj <- rbioClass_svm_ncv_fs(x = x, y = y, cross.k = 2,
                              univariate.fs = TRUE, uni.contrast = "B-A",
                              uni.log2trans = FALSE, uni.alpha = 0.5,
                              rf.ifs.ntree = 50, rf.sfs.ntree = 50,
                              tune.method = "cross", tune.cross.k = 3,
                              parallelComputing = FALSE, verbose = FALSE)
  for (i in 1:2) {
    expect_false(is.null(obj$nested.cv.models[[i]]$uni.sig.fs))
    expect_type(obj$nested.cv.models[[i]]$uni.sig.fs, "character")
  }
  # both informative features should be retained by the univariate filter
  all_uni <- sort(unique(unlist(lapply(obj$nested.cv.models, `[[`, "uni.sig.fs"))))
  expect_true(all(c("sig1", "sig2") %in% all_uni))
})

# T5: v3 outer folds keep groups intact
test_that("T5: ncv_fs_v3 outer folds keep sampleIds together", {
  set.seed(2026)
  grp <- rep(1:12, each = 2)           # 12 groups x 2 rows
  y <- factor(c(rep("A", 12), rep("B", 12)))  # groups 1-6 A, 7-12 B
  x <- matrix(rnorm(24 * 4), 24, 4)
  sids <- as.character(grp)
  obj <- rbioClass_svm_ncv_fs_v3(x = x, y = y, cross.k = 4, sampleIds = sids,
                                 rf.ifs.ntree = 50, rf.sfs.ntree = 50,
                                 tune.method = "cross", tune.cross.k = 3,
                                 parallelComputing = FALSE, verbose = FALSE)
  # each group maps to exactly one outer fold
  expect_true(all(vapply(unique(sids), function (g)
    length(unique(obj$cv.fold[sids == g])) == 1L, logical(1))))
  # mapping stored on the object is consistent with cv.fold
  expect_equal(unname(obj$sampleIdToFold[obj$sampleIds]), as.integer(obj$cv.fold))
  # all folds non-empty
  expect_equal(length(unique(obj$cv.fold)), 4)
})

# T6: zcore_minmax performance sanity (compiled matrixStats + sweep path)
test_that("T6: zcore_minmax is competitive with base::scale", {
  set.seed(2026)
  big <- matrix(rnorm(500 * 2000), 500, 2000)
  t_scale <- system.time(base::scale(big))["elapsed"]
  t_zm    <- system.time(zcore_minmax(big, zscore = TRUE, minmax = FALSE, verbose = FALSE))["elapsed"]
  # generous margin to stay robust on slow CI machines; the point is that the
  # matrixStats-based path is not an order of magnitude slower
  expect_lt(t_zm, max(3 * t_scale, 5))
})
