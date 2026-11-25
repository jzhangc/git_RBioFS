# Regression tests for the explicit uni.log2trans = TRUE NaN guard.
#
# log2(x + 2) is undefined for input values below -2. Before the guard, such
# values produced NaN entries in the expression matrix, NA P.Values in the
# limma topTable, and then (a) silently dropped affected features, or (b) if
# every feature was affected, the generic "No statistically significant
# features found" stop that hid the real cause. The guard now fails fast with
# a message that points at the input data.

# 2-class design where EVERY feature value is below -2, so every
# log2(x + 2) entry is NaN and (pre-guard) all limma P.Values would be NA.
make_all_below_neg2_data <- function(seed = 11) {
  set.seed(seed)
  n <- 24; p <- 4
  x <- matrix(runif(n * p, min = -10, max = -3), n, p,
              dimnames = list(NULL, paste0("f", 1:p)))
  y <- factor(rep(c("A", "B"), each = n / 2))
  list(x = x, y = y)
}

test_that("uni.log2trans = TRUE errors clearly when ALL features are < -2 (all-NA P.Value path)", {
  d <- make_all_below_neg2_data()
  expect_true(all(d$x < -2))
  expect_error(
    rbioClass_svm_ncv_fs(x = d$x, y = d$y, cross.k = 2,
                         univariate.fs = TRUE, uni.log2trans = TRUE, uni.contrast = "B-A",
                         tune.method = "fix",
                         parallelComputing = FALSE, verbose = FALSE),
    regexp = "Check the input data for issues"
  )
  # the message must name the offending transform, not just "no features"
  expect_error(
    rbioClass_svm_ncv_fs(x = d$x, y = d$y, cross.k = 2,
                         univariate.fs = TRUE, uni.log2trans = TRUE, uni.contrast = "B-A",
                         tune.method = "fix",
                         parallelComputing = FALSE, verbose = FALSE),
    regexp = "uni.log2trans"
  )
})

test_that("uni.log2trans = TRUE fails fast when only SOME features are < -2", {
  d <- make_all_below_neg2_data()
  # one feature is fine, three are below -2
  d$x[, 1] <- runif(nrow(d$x), min = 0, max = 10)
  expect_error(
    rbioClass_svm_ncv_fs(x = d$x, y = d$y, cross.k = 2,
                         univariate.fs = TRUE, uni.log2trans = TRUE, uni.contrast = "B-A",
                         tune.method = "fix",
                         parallelComputing = FALSE, verbose = FALSE),
    regexp = "non-finite values .* for 3 feature"
  )
})

test_that("uni.log2trans = FALSE is unaffected by sub--2 values (control)", {
  withr::local_dir(withr::local_tempdir())  # rRF-FS writes iteration artifacts to the CWD
  d <- make_all_below_neg2_data()
  # Engineered class difference that survives uni_PreProc's quantile
  # normalization (which erases raw-scale shifts with matching ranks):
  #   class A: f1 rank 1 (50-60, above everything)
  #   class B: f2 rank 1 (-2.5..-1), f1 rank 2 (-5..-3), f3/f4 rank 3-4
  # so f1 maps to different quantile targets per class -> significant.
  d$x[1:12, 1] <- runif(12, 50, 60)
  d$x[13:24, 1] <- runif(12, -5, -3)
  d$x[13:24, 2] <- runif(12, -2.5, -1)
  d$x[13:24, 3:4] <- runif(24, -10, -6)
  res <- rbioClass_svm_ncv_fs(x = d$x, y = d$y, cross.k = 2,
                              univariate.fs = TRUE, uni.log2trans = FALSE, uni.contrast = "B-A",
                              tune.method = "fix",
                              parallelComputing = FALSE, verbose = FALSE)
  expect_type(res, "list")
  expect_length(res$nested.cv.models, 2)
  expect_true("f1" %in% res$nested.cv.models[[1]]$uni.sig.fs)
})
