# Targeted tests for the non-leakage-inner-CV feature (stratify / sampleIds)
# wired through RBioFS's SVM functions.  These cover the 2.4 acceptance
# criteria; the historical suite only exercises the PR-AUC path.
test_that("2.4(a) stratify=TRUE on imbalanced 2-class data does not crash and records stratify", {
  set.seed(2026)
  # very imbalanced 2-class design: 8 A, 42 B
  y <- factor(c(rep("A", 8), rep("B", 42)))
  x <- as.matrix(cbind(
    A = c(rnorm(8, 1, 1), rnorm(42, -1, 1)),
    B = rnorm(50),
    C = rnorm(50)))
  res <- rbioClass_svm(x = x, y = y, kernel = "radial",
                        svm.cross.k = 0, tune.method = "cross",
                        tune.cross.k = 5, stratify = TRUE, verbose = FALSE)
  expect_s3_class(res, "rbiosvm")
  expect_true(isTRUE(res$stratify))
  expect_true(is.numeric(res$cost))
  expect_true(is.numeric(res$gamma))
})

test_that("2.4(b) group-aware sampleIds keeps rows of one id together across inner-CV folds", {
  set.seed(2027)
  # 12 groups, 2 rows each = 24 rows; groups 1-6 class A, groups 7-12 class B
  grp <- rep(1:12, each = 2)
  y   <- factor(rep(c("A", "B"), each = 12))   # 12 rows A (groups 1-6), 12 rows B (groups 7-12)
  x   <- matrix(rnorm(24 * 4), 24, 4)
  sids <- as.character(grp)
  # the e1071mc helper must keep each id in a single fold
  fid <- e1071mc:::.svm_mc_assign_folds(24, 5, y = y, sampleIds = sids,
                                         stratify = TRUE, regression = FALSE, seed = 2027)
  no_split <- all(vapply(unique(sids), function (g)
                       length(unique(fid[sids == g])) == 1L, logical(1)))
  expect_true(no_split)
  # rbioClass_svm must run with sampleIds and record it
  res <- rbioClass_svm(x = x, y = y, kernel = "radial", svm.cross.k = 0,
                        tune.method = "cross", tune.cross.k = 5,
                        sampleIds = sids, verbose = FALSE)
  expect_s3_class(res, "rbiosvm")
  expect_equal(length(res$sampleIds), 24)
  expect_true(is.numeric(res$cost))
})

test_that("2.4(c) minmax=TRUE tuned params equal a manual scale=FALSE inner tune (scale-fix regression gate)", {
  y <- factor(c(rep("A", 30), rep("B", 30)))
  # Manual reproduction mirrors rbioClass_svm's exact RNG sequence:
  #   set.seed -> rnorm(x) -> zcore_minmax (deterministic, no RNG) -> tune (sample perm)
  # so the inner CV folds are identical and the scale=FALSE inner tune must match.
  set.seed(2028)
  x <- matrix(rnorm(60 * 6), 60, 6)
  centered <- zcore_minmax(x, zscore = FALSE, minmax = TRUE, verbose = FALSE)
  X <- centered$centerX
  gamma_start <- 1 / ncol(X)
  wgt <- if (length(unique(table(y))) != 1L) (length(y) / table(y)) else (unique(table(y)) / table(y))
  manual <- e1071mc::tune(e1071mc::svm, kernel = "radial", scale = FALSE,
                           train.x = X, train.y = y,
                           ranges = list(gamma = gamma_start * 2^(-5:5), cost = 2^(-10:8)),
                           class.weights = wgt,
                           tunecontrol = e1071mc::tune.control(
                             sampling = "cross", cross = 5, stratify = TRUE))
  bp_manual <- manual$best.parameters   # a 1-row data.frame (gamma, cost)
  # rbioClass_svm with the same seed -> same transform (deterministic) + same inner tune
  set.seed(2028)
  x2 <- matrix(rnorm(60 * 6), 60, 6)
  res <- rbioClass_svm(x = x2, y = y, kernel = "radial", center.scale = TRUE,
                        zscore = FALSE, minmax = TRUE, svm.cross.k = 0,
                        tune.method = "cross", tune.cross.k = 5,
                        stratify = TRUE, verbose = FALSE)
  # scale-fix gate: the model's cost/gamma come straight from the scale=FALSE inner tune
  expect_equal(unname(bp_manual$cost), unname(res$cost))
  expect_equal(unname(bp_manual$gamma), suppressWarnings(as.numeric(res$gamma)))
})

test_that("2.4(d) n_cores=1 (serial tune) vs n_cores=2 (parallel tune_mc) give identical cost/gamma", {
  y <- factor(c(rep("A", 30), rep("B", 30)))
  set.seed(2029); x <- matrix(rnorm(60 * 6), 60, 6)
  args <- list(x = x, y = y, kernel = "radial", svm.cross.k = 0,
               tune.method = "cross", tune.cross.k = 5,
               stratify = TRUE, verbose = FALSE)
  set.seed(2029); m1 <- do.call(rbioClass_svm, c(list(n_cores = 1L), args))
  set.seed(2029); m2 <- do.call(rbioClass_svm, c(list(n_cores = 2L), args))
  expect_true(isTRUE(all.equal(m1$cost, m2$cost)))
  expect_true(isTRUE(all.equal(suppressWarnings(m1$gamma), suppressWarnings(m2$gamma))))
})
