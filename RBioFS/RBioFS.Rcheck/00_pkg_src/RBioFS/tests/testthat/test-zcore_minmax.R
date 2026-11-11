# Tests for zcore_minmax()
#
# Drop-in, faster replacement for center_scale() with z-score and/or min-max
# transforms. Exercises: matrix / data.frame / vector inputs; z-score only,
# min-max only, and the combined z-score -> min-max order; re-applying the
# fitted parameters to new data; equivalence with the legacy center_scale;
# and the error / warning edge cases (non-numeric data, missing values,
# zero-SD, zero-range, both flags FALSE).

# Re-apply the stored parameters of a zcore_minmax result to new data, mirroring
# the documented (z-score first, then min-max) pipeline.
apply_fitted <- function(res, new) {
  z <- new
  if (res$zscore) {
    z <- sweep(sweep(z, 2, res$meanX, "-"), 2, res$columnSD, "/")
  }
  if (res$minmax) {
    z <- sweep(sweep(z, 2, res$colMin, "-"), 2, (res$colMax - res$colMin), "/")
  }
  z
}

# Small, varied fixture data: rows = samples, cols = features. Different raw
# scales/spreads per column so the transforms are meaningful.
make_fixture <- function(seed = 7, n = 24, p = 5) {
  set.seed(seed)
  X <- matrix(rnorm(n * p, mean = 10, sd = 3), n, p,
              dimnames = list(paste0("s", 1:n), paste0("f", 1:p)))
  X[, 3] <- X[, 3] * 10       # very different scale
  X[, 5] <- X[, 5] + 100      # very different offset
  X
}

# ----- 1. z-score only (drop-in for center_scale(scale = TRUE)) -----
test_that("z-score only: structure, names, and [0,1]-free standardized output", {
  X <- make_fixture()
  res <- zcore_minmax(X, zscore = TRUE, minmax = FALSE, verbose = FALSE)

  expect_s3_class(res, "zcore_minmax")
  expect_setequal(names(res), c("centerX", "zscore", "meanX", "columnSD",
                                "minmax", "colMin", "colMax", "n.na"))
  expect_true(res$zscore)
  expect_false(res$minmax)
  expect_null(res$colMin)
  expect_null(res$colMax)
  # standardized columns have ~mean 0 and ~SD 1 (ignore names on the comparison)
  expect_equal(unname(colMeans(res$centerX, na.rm = TRUE)),
               rep(0, ncol(X)), tolerance = 1e-8)
  expect_equal(unname(apply(res$centerX, 2, sd)), rep(1, ncol(X)), tolerance = 1e-8)
})

# ----- 2. data.frame input, all-numeric -----
test_that("data.frame input is coerced and transformed correctly", {
  X <- make_fixture()
  df <- as.data.frame(X)
  res <- zcore_minmax(df, zscore = TRUE, minmax = FALSE, verbose = FALSE)
  # identical result to the matrix path
  res_mat <- zcore_minmax(X, zscore = TRUE, minmax = FALSE, verbose = FALSE)
  expect_equal(res$centerX, res_mat$centerX)
  expect_equal(colnames(res$centerX), paste0("f", 1:5))
})

# ----- 3. vector input -> single column -----
test_that("numeric vector input is treated as a single feature column", {
  v <- c(1, 2, 3, 4, 5, 6)
  res <- zcore_minmax(v, zscore = TRUE, minmax = FALSE, verbose = FALSE)
  expect_equal(unname(dim(res$centerX)), c(6, 1))
  expect_equal(unname(colMeans(res$centerX)), 0, tolerance = 1e-8)
  expect_equal(unname(apply(res$centerX, 2, sd)), 1, tolerance = 1e-8)
})

# ----- 4. min-max only -----
test_that("min-max only: each column mapped onto [0, 1]", {
  X <- make_fixture()
  res <- zcore_minmax(X, zscore = FALSE, minmax = TRUE, verbose = FALSE)
  expect_null(res$meanX)
  expect_null(res$columnSD)
  expect_true(res$minmax)
  expect_equal(unname(apply(res$centerX, 2, min)), rep(0, ncol(X)), tolerance = 1e-10)
  expect_equal(unname(apply(res$centerX, 2, max)), rep(1, ncol(X)), tolerance = 1e-10)
})

# ----- 5. combined: z-score THEN min-max -----
test_that("combined transform applies z-score first, then min-max, yielding [0,1]", {
  X <- make_fixture()
  res <- zcore_minmax(X, zscore = TRUE, minmax = TRUE, verbose = FALSE)
  expect_true(res$zscore)
  expect_true(res$minmax)
  # min-max is the last step, so the output is in [0, 1]
  expect_true(all(res$centerX >= -1e-10 & res$centerX <= 1 + 1e-10))
  expect_equal(unname(apply(res$centerX, 2, min)), rep(0, ncol(X)), tolerance = 1e-10)
  expect_equal(unname(apply(res$centerX, 2, max)), rep(1, ncol(X)), tolerance = 1e-10)
  # the stored min/max bounds are the min/max of the z-scored (intermediate) data,
  # i.e. not of the raw input -> confirms z-score ran first
  Xz <- sweep(sweep(X, 2, colMeans(X, na.rm = TRUE), "-"),
              2, apply(X, 2, sd), "/")
  expect_equal(res$colMin, apply(Xz, 2, min), tolerance = 1e-8)
  expect_equal(res$colMax, apply(Xz, 2, max), tolerance = 1e-8)
})

# ----- 6. re-applying fitted params reproduces the transform -----
test_that("fitted parameters re-applied to the same data reproduce centerX", {
  X <- make_fixture(seed = 11)
  res <- zcore_minmax(X, zscore = TRUE, minmax = TRUE, verbose = FALSE)
  reapplied <- apply_fitted(res, X)
  expect_equal(reapplied, res$centerX, tolerance = 1e-10)
})

# ----- 7. equivalence with the legacy center_scale (z-score only) -----
test_that("z-score only matches the legacy center_scale(scale = TRUE) idiom", {
  X <- make_fixture(seed = 3)
  res <- zcore_minmax(X, zscore = TRUE, minmax = FALSE, verbose = FALSE)
  # legacy computation
  m <- colMeans(X, na.rm = TRUE)
  s <- matrixStats::colSds(X, center = m, na.rm = TRUE)
  legacy <- t((t(X) - m) / s)
  expect_equal(res$centerX, legacy, tolerance = 1e-10)
  # the returned params match too
  expect_equal(res$meanX, m, tolerance = 1e-12)
  expect_equal(res$columnSD, s, tolerance = 1e-12)
})

# ----- 8. error: non-numeric data.frame -----
test_that("non-numeric column in a data.frame raises an error", {
  X <- make_fixture()
  df <- as.data.frame(X)
  df$grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  df$note <- letters[1:nrow(df)]
  expect_error(zcore_minmax(df, verbose = FALSE), "non-numeric")
})

# ----- 9. error: NA with na.rm = FALSE -----
test_that("missing values error out when na.rm = FALSE", {
  X <- make_fixture()
  X[2, 3] <- NA
  expect_error(zcore_minmax(X, na.rm = FALSE, verbose = FALSE), "missing value")
})

# ----- 10. warning: NA with na.rm = TRUE, and n.na recorded -----
test_that("missing values are tolerated (with warning) when na.rm = TRUE", {
  X <- make_fixture()
  X[2, 3] <- NA
  X[5, 5] <- NA
  res <- suppressWarnings(zcore_minmax(X, zscore = TRUE, na.rm = TRUE, verbose = FALSE))
  expect_warning(zcore_minmax(X, zscore = TRUE, na.rm = TRUE, verbose = FALSE),
                 "missing value")
  expect_equal(res$n.na, 2)
})

# ----- 11. error: both transforms requested to be FALSE -----
test_that("both zscore = FALSE and minmax = FALSE is rejected", {
  X <- make_fixture()
  expect_error(zcore_minmax(X, zscore = FALSE, minmax = FALSE, verbose = FALSE),
               "at least one")
})

# ----- 12. warning: near-zero SD column -----
test_that("a constant (zero-SD) column triggers a scaling warning but still runs", {
  X <- make_fixture()
  X[, 2] <- 7   # constant feature -> SD 0
  res <- suppressWarnings(zcore_minmax(X, zscore = TRUE, minmax = FALSE, verbose = FALSE))
  expect_warning(zcore_minmax(X, zscore = TRUE, minmax = FALSE, verbose = FALSE),
                 "zero standard deviation")
  # the constant column becomes all-NaN (0/0) and is not silently dropped
  expect_true(all(is.na(res$centerX[, 2])))
})

# ----- 13. warning: near-zero range column under min-max -----
test_that("a constant (zero-range) column triggers a min-max warning", {
  X <- make_fixture()
  X[, 4] <- 12   # constant feature -> min == max
  expect_warning(zcore_minmax(X, zscore = FALSE, minmax = TRUE, verbose = FALSE),
                 "zero range")
})

# ----- 14. print method -----
test_that("print method produces output and returns the object invisibly", {
  X <- make_fixture()
  res <- zcore_minmax(X, zscore = TRUE, minmax = TRUE, verbose = FALSE)
  captured <- capture.output(print(res))
  expect_true(any(grepl("zcore_minmax result", captured)))
  expect_true(any(grepl("z-score -> min-max", captured)))
  # invisible: print returns the object
  expect_identical(print(res), res)
})

# ----- 15. error: non-matrix / non-df / non-numeric type -----
test_that("non-matrix, non-data.frame, non-vector input is rejected", {
  # a list is, in R, a (non-numeric) vector -> rejected as non-numeric
  expect_error(zcore_minmax(list(1, 2, 3), verbose = FALSE),
               "input vector must be numeric")
  # a character vector is also non-numeric
  expect_error(zcore_minmax("abc", verbose = FALSE), "input vector must be numeric")
  # a true non-vector (a closure) is not an accepted container -> generic message
  expect_error(zcore_minmax(function() 1, verbose = FALSE),
               "must be a matrix, data.frame or numeric vector")
})

# ----- 16. empty input rejected -----
test_that("empty matrix input is rejected", {
  empty <- matrix(numeric(0), 0, 0)
  expect_error(zcore_minmax(empty, verbose = FALSE), "non-empty")
})

# ----- 17. matrix column names are preserved -----
test_that("existing column names are preserved in the output", {
  X <- make_fixture()
  res <- zcore_minmax(X, zscore = TRUE, verbose = FALSE)
  expect_equal(colnames(res$centerX), paste0("f", 1:5))
})

# ----- 18. unnamed matrix gets default feature names -----
test_that("an unnamed matrix receives default feature names", {
  X <- matrix(rnorm(20 * 4), 20, 4)
  res <- zcore_minmax(X, zscore = TRUE, verbose = FALSE)
  expect_equal(colnames(res$centerX), paste0("feature_", 1:4))
})

# ================= .apply_zcore_minmax() re-application helper =================

# ----- 19. re-applying fitted params to the SAME data reproduces centerX -----
test_that(".apply_zcore_minmax on the same data reproduces the fitted transform (z/mm/both)", {
  X <- make_fixture()
  z    <- zcore_minmax(X, zscore = TRUE,  minmax = FALSE, verbose = FALSE)
  mm   <- zcore_minmax(X, zscore = FALSE, minmax = TRUE,  verbose = FALSE)
  both <- zcore_minmax(X, zscore = TRUE,  minmax = TRUE,  verbose = FALSE)

  expect_equal(RBioFS:::.apply_zcore_minmax(X, params = z,    verbose = FALSE),
               z$centerX, tolerance = 1e-9)
  expect_equal(RBioFS:::.apply_zcore_minmax(X, params = mm,   verbose = FALSE),
               mm$centerX, tolerance = 1e-9)
  expect_equal(RBioFS:::.apply_zcore_minmax(X, params = both, verbose = FALSE),
               both$centerX, tolerance = 1e-9)
})

# ----- 20. NULL / identity: NULL params returns the input unchanged -----
test_that(".apply_zcore_minmax with NULL params returns input unchanged (center.scale = FALSE path)", {
  X <- make_fixture()
  out <- RBioFS:::.apply_zcore_minmax(X, params = NULL, verbose = FALSE)
  expect_equal(out, X)
})

# ----- 21. legacy object (minmax flagged but no colMin/colMax): z-score only + warning -----
test_that("legacy fitted object (minmax flag TRUE, no bounds) falls back to z-score only and warns", {
  X <- make_fixture()
  z <- zcore_minmax(X, zscore = TRUE, minmax = FALSE, verbose = FALSE)
  # mimic an older object: z-score params present, minmax flag set, but no min/max bounds
  legacy <- list(zscore = TRUE, minmax = TRUE,
                 meanX = z$meanX, columnSD = z$columnSD,
                 colMin = NULL, colMax = NULL)

  out <- suppressWarnings(RBioFS:::.apply_zcore_minmax(X, params = legacy, verbose = FALSE))
  expect_equal(out, z$centerX, tolerance = 1e-9)  # z-score-only result
  expect_warning(
    RBioFS:::.apply_zcore_minmax(X, params = legacy, verbose = FALSE),
    "colMin/colMax"
  )
})

# ----- 22. the combined z-score + min-max is order-stable: == min-max-only (affine) -----
test_that("z-score then min-max equals min-max only (affine invariance), both in [0,1]", {
  X <- make_fixture()
  mm   <- zcore_minmax(X, zscore = FALSE, minmax = TRUE, verbose = FALSE)
  both <- zcore_minmax(X, zscore = TRUE,  minmax = TRUE, verbose = FALSE)
  expect_equal(both$centerX, mm$centerX, tolerance = 1e-9)
  expect_true(all(both$centerX >= -1e-9 & both$centerX <= 1 + 1e-9))
})

# ----- 23. 1-D numeric vector input is coerced to a single-column matrix -----
test_that(".apply_zcore_minmax accepts a 1-D numeric vector (single-feature fold)", {
  v <- rnorm(24, mean = 5, sd = 2)          # single feature as a plain vector
  z <- zcore_minmax(v, zscore = TRUE, minmax = FALSE, verbose = FALSE)  # fits a 1-col
  out <- RBioFS:::.apply_zcore_minmax(v, params = z, verbose = FALSE)   # re-apply, input is a vector
  expect_true(is.matrix(out))
  expect_equal(unname(out), unname(z$centerX), tolerance = 1e-9)
  # data.frame input is also coerced
  df <- data.frame(x = v, stringsAsFactors = FALSE)
  out2 <- RBioFS:::.apply_zcore_minmax(df, params = z, verbose = FALSE)
  expect_equal(unname(out2), unname(z$centerX), tolerance = 1e-9)
})
