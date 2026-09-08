library(testthat)
library(RBioFS)

# Standard testthat runner invoked by `R CMD check` (and by:
#   Rscript tests/testthat.R)
#
# The function under test (rbioClass_svm_cv_pr_auc) is loaded from the
# installed RBioFS package.
test_check("RBioFS")
