# Session Summary — Non-leakage-capable inner/final CV for `e1071mc` + `RBioFS`

## What was asked
Verify Part 1 (`e1071mc` fork) and complete Part 2 (`RBioFS`), wiring `stratify`/`sampleIds`
through the SVM functions to make inner-tuning CV and final-model CV non-leakage-capable.
`stratify` defaults `TRUE`; **serial/parallel parity** is the #1 acceptance gate.

## What I found
Both parts were **already substantially implemented in the working tree**. My work was
*verification + completion + one fix*, plus writing the new targeted tests.

## Part 1 — `e1071mc` (verified)
- Shared helper `.svm_mc_assign_folds()` (`R/e1071mc_cv_fold.R`) and all call sites
  (`tune.R`, `tune_mc.R`, `svm_multicore.R`) were correct. **Not exported** (internal) — correct.
  `parallel` already imported.
- `R CMD INSTALL` clean at v`1.7-18`.
- **Verification (a–e), all PASS:**
  - **(a) backward-compat:** I discovered the spec's "Else" formula `((perm-1)%%k)+1` was a
    **mis-statement** — the true HEAD baseline is the contiguous `cut(1:n, breaks=k)`
    (confirmed in `git show HEAD:R/tune.R:78`). The implementation's default path reproduces
    HEAD **exactly**: fold IDs + `best.performance` identical across 3 seeds.
  - **(b) `stratify=TRUE`:** imbalanced 8:A/42:B → every fold has both classes (1/2/7); the
    legacy `cut` path *does* produce empty-class folds (the bug being fixed).
  - **(c) `sampleIds`:** no group split across folds; all folds populated.
  - **(d) serial vs parallel:** identical `best.parameters` for default, `stratify=TRUE`, and
    `sampleIds` modes.
  - **(e)** `tune_mc$performances` shape unchanged (`gamma,cost,error,dispersion`; 9 rows for
    a 3×3 grid).
- **Fix applied:** `inst/NEWS.Rd` heading was `1.7-17-2-20261020` → changed to **`1.7-18`**
  to match DESCRIPTION/README/the RBioFS dependency.
- `R CMD check`: only ERROR is the **pre-existing** `svm_mc` example bug `y <- Species`
  (HEAD `man/svm_mc.Rd:153-154`, untouched); **tests pass** (`clustering.R OK`).

## Part 2 — `RBioFS` (completed + verified)
- 2.1–2.3 already in working tree and correct:
  - 2.1 `rbioClass_svm`: flags after `minmax`, inner `tune`/`tune_mc` use **`scale=FALSE`**
    (minmax fix) + thread flags; final `svm`/`svm_mc` get them; `m$stratify`/`m$sampleIds`
    stored; all validations present.
  - 2.2 callers (`ncv_fs`, `_v2`, `_v3`, `cv`, `legacy.R`) thread flags; **v3 correctly reuses
    its own `sampleIds`** for group-aware inner CV; `rbioClass_svm_perm` correctly untouched.
  - 2.3 `e1071mc (>= 1.7-18)`, `Remotes:` present, `Version: 0.8.3-20261031`.
- **2.4 done:**
  - `devtools::document()` regenerated man cleanly; `R CMD build` both packages clean.
  - **New test file** `tests/testthat/test-rbioClass_svm_nonleakage_cv.R` — **4 tests /
    10 expectations, all PASS** (a: imbalanced `stratify=TRUE` runs + `$stratify==TRUE`;
    b: group IDs unsplit; c: minmax `scale=FALSE` regression gate; d: `n_cores=1` vs `2`
    identical cost/gamma).
  - Full suite: **110 PASS / 2 FAIL** — the 2 are the **pre-existing**
    `rbioClass_svm_cv_roc_auc`/`_v2` absent-class tests (per AGENTS.md, unrelated).
    Examples passed.

## Cross-cutting — all PASS
- **Parity (3 modes):** serial `tune` ≡ parallel `tune_mc` (default, `stratify`, `sampleIds`).
- **Backward-compat:** new default ≡ HEAD `cut` folds + performance (3 seeds).
- RBioFS `(d)`: `n_cores` 1 vs 2 identical cost/gamma.

## Caveats (flagged)
1. **Version conflict (spec-internal).** Spec 1.7 wanted `1.7-17-2-20261020`, but spec 2.3
   required the dep `e1071mc (>= 1.7-18)`. In R ordering `1.7-17-2-20261020 < 1.7-18`, so
   20261020 can't satisfy 2.3. The working tree + installed package already used
   **`1.7-18`** consistently — I aligned to it (this deviates from spec 1.7's literal string;
   **flagging**).
2. **Final-model `svm_mc` CV is not 1-vs-≥2-core parity-safe** when flags are set:
   `n_cores=1` uses the C cross engine in `svm.default` (ignores flags); `n_cores≥2` uses
   R `.svm_mc_cross_validate` (honors them). **Pre-existing architectural asymmetry
   (C vs R engine), not introduced here**, and not the #1 target (the inner-tuning
   `tune`/`tune_mc` path is fully parity-clean). Nested-FS wrappers use `svm.cross.k=0`,
   so this branch doesn't fire there.
3. `mclapply` is UNIX-only (existing limitation; serial path is the Windows/PSOCK fallback).
4. `stratify=TRUE` default is a **behaviour change** (tuned params / CV metrics /
   selected features may shift) — documented in README + NEWS.
5. Group-aware CV can inflate RMSE / lower accuracy vs sample-level — correct non-leakage
   estimate, not a bug.

## Status
No commits made (working-tree changes only). Two stray machine-specific files left untracked:
`e1071mc/src/Makevars` (local build flags, not part of this deliverable) and the built tarballs.
