# AGENTS.md

R package `RBioFS`: wrapper for recursive random-forest feature selection (RF-FS) plus SVM/PLS/PCA classification & regression workflows for biological/biomedical data.

## Layout

- Package source is in the **`RBioFS/` subdir**, not the repo root. All `R/`, `man/`, `DESCRIPTION`, `NAMESPACE` live under `RBioFS/`.
- `.Rbuildignore` (both root and `RBioFS/`) excludes `.Rproj`, `.Rproj.user`, `.DS_Store` from the built tarball.

## Dependencies & install (high-signal, easy to get wrong)

- **GitHub-only forks must be listed in `Remotes:`** in `RBioFS/DESCRIPTION` — CRAN/Bioconductor installers will otherwise fail to find them:
  - `jzhangc/e1071mc` (a multi-core fork of CRAN `e1071`; e1071mc is **not** on CRAN/Bioconductor).
  - `jzhangc/git_R_STATS_KBS/package/rbioplot` (fork of `RBioplot`).
- Most other deps are Bioconductor (`limma`, `edgeR`, `klaR`) — setting `options(repos = BiocManager::repositories())` and running `BiocManager::install()` is **required** before installing.

Install (use the `nightly` branch unless a stable/beta build is intended):

```r
# pak (recommended); build encoded in URL, no `ref` arg
options(repos = BiocManager::repositories())
pak::pkg_install("github::jzhangc/git_RBioFS/RBioFS@nightly")

# devtools; branch via `ref` (HEAD = stable)
devtools::install_github("jzhangc/git_RBioFS/RBioFS",
    repos = BiocManager::repositories(), ref = "nightly")
```

Build branches: **`HEAD` = stable**, `beta` = pre-stable, `nightly` = active dev source.

## Documentation / build

- `man/` and `NAMESPACE` are **roxygen2-generated**. Edit roxygen comments in `RBioFS/R/*.R`, then regenerate:
  ```r
  devtools::document()   # or: Rscript -e 'roxygen2::roxygenise()'
  ```
  Do **not** hand-edit `RBioFS/man/*.Rd` or `RBioFS/NAMESPACE`.
- Build/install from source:
  ```
  R CMD build RBioFS/       # produces RBioFS_<version>.tar.gz (gitignored)
  R CMD INSTALL RBioFS/      # or install the tarball
  ```
- `Config/roxygen2/version: 8.1.0` — keep roxygen2 >= that when regenerating.

## Validation

- **No test suite, no linter, no CI.** `RBioFS/tests/` does not exist; there is no `.lintr` or lint/format config.
- Validate changes with `R CMD check RBioFS_<version>.tar.gz` and/or manual checks. Do not claim tests pass when none run.

## Conventions

- Function name prefixes map to subpackages: `rbioClass_svm*` (SVM), `rbioFS_rf_*` / `rbioFS_*` (RF-FS & FS), `rbioClass_plsda*` / `rbioReg_plsr*` (PLS-DA/PLSR), `rbioFS_PCA*` (PCA), `rbioUtil_*` (helpers). Shiny apps end in `_app`.
- Multi-core support is baked into the existing SVM functions via an `n_cores` argument, not a separate `_mc` function: `rbioClass_svm` (and the ncv/fs family) use `e1071mc::svm_mc()` / `e1071mc::tune_mc()` when `n_cores > 1`, falling back to `e1071mc::svm()` / `tune()` when `n_cores <= 1` (or `NULL`).
- `rbioClass_svm_ncv_fs` has an explicit `v2` and `v3` variant, plus a `rbioClass_svm_ncv_fs_legacy` variant (all four are exported). `v3` adds sample-ID-aware fold assignment (samples sharing a `sampleIds` ID are kept in the same fold, with group-level stratification for classification and variance-minimising assignment for regression); it errors if samples with the same `sampleIds` carry different class labels.
- Group-fold helpers `.assign_group_folds_stratified` (classification) and `.assign_group_folds_simple` (regression) live in `RBioFS/R/rbioUtil.R` and drive `rbioClass_svm_ncv_fs_v3`.
- S3 output classes: `rbiosvm`, `rbiosvm_nestedcv`, `rbiomvr`, `rf_ifs`, `rf_sfs`, `rbiofs_pca`, `svm_roc_auc`, `prediction`, etc. Many functions set `model.type` = `"classification"` or `"regression"` on the returned object; some functions reject models of the other type.
- `verbose = TRUE` by default on most functions; it controls messages, not errors/warnings.
- `center.scale` (a.k.a. Z-score / min-max) is a pervasive preprocessing flag across SVM/PLS/PCA functions; mind consistency between train and predict.
- `zzz.R` `.onAttach` prints a startup message with citation info; keep it intact.

## R syntax styles

- DO NOT use tidyverse syntax unless otherwise requested

## Gotchas

- `AGENTS.md` and `prompts.md` are **gitignored** (see root `.gitignore`) — they will not be committed.
- `e1071` is **no longer** a dependency — SVM code relies solely on `e1071mc` (a multi-core fork/superset of CRAN `e1071`). All `@importFrom` and `:::` calls point at `e1071mc`. Watch for namespace/`tune`/`svm` clashes if `e1071` reappears when touching SVM code.
- `prompts.md` holds the active dev task list (e1071mc migration, the `n_cores` multi-core SVM path, `rbioClass_svm_ncv_fs_v3` group-fold CV); consult it before changing SVM/FS code.
