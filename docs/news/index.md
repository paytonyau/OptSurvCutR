# Changelog

## OptSurvCutR v0.9.9 (2026-06)

### Major Engine Enhancements

- **Integer Index-Space Mapping:** Replaced the continuous
  floating-point search space with a discrete, bounded integer lattice
  mapped directly to sorted unique data row indices. This fundamentally
  shifts the engine math from an infinite decimal space to a finite
  spectrum of actual observations, eliminating micro-decimal overfitting
  and stochastic seed drift.
- **Hybrid Execution Routing:** Upgraded threshold discovery to
  dynamically select the most efficient search engine. The package now
  defaults to a deterministic, exhaustive **`method = "systematic"`**
  grid sweep for lower dimensional profiles (1 or 2 cuts) to guarantee
  global maximum discovery with zero performance lag. For higher
  dimensions ($`K \ge 3`$), it seamlessly transitions to the
  evolutionary genetic algorithm (`rgenoud` backend).
- **Automated Proportional Hazards Diagnostics:** Integrated a native
  2-tier Schoenfeld residuals validation layer inside the core S3
  [`summary()`](https://rdrr.io/r/base/summary.html) methods. The
  package automatically flags models violating the Cox proportional
  hazards assumption, distinguishing between *Tier 1 (Proportional,
  Stable)* and *Tier 2 (Time-Varying, Dynamic Risk)* cohorts.
- **Strict Structural Headroom Guardrails:** Upgraded `nmin` cell-floor
  checking to calculate total available degrees of freedom prior to loop
  execution. When discrete data density collapses or drops below
  necessary cell floors, the engine throws clean UX alert messages and
  skips impossible configurations to prevent fatal R session crashes.

------------------------------------------------------------------------

#### 📊 Algorithmic Scaling Protocol

The table below outlines how `OptSurvCutR` automatically routes
execution and scales its internal tuning parameters (`pop.size` and
`max.generations`) based on target complexity ($`K`$ Cuts) and data
structure:

| Target Complexity ($`K`$ Cuts) | Execution Method | Default Population Size (`pop.size`) | Default Search Lifespan (`max.generations`) | Optimization Mechanics |
|----|----|----|----|----|
| **$`K = 1`$ Cut** | `systematic` | *N/A (Exhaustive Grid)* | *N/A (Exhaustive Grid)* | Full 1D Coordinate Vector Sweep |
| **$`K = 2`$ Cuts** | `systematic` | *N/A (Exhaustive Grid)* | *N/A (Exhaustive Grid)* | Full 2D Cross-Lattice Matrix Sweep |
| **$`K = 3`$ Cuts** | `genetic` | 100 | 50 | 3D Index Hyper-Lattice Traversal |
| **$`K = 4`$ Cuts** | `genetic` | 120 | 55 | 4D Hyper-Lattice Spatial Search |
| **$`K = 5`$ Cuts** | `genetic` | 150 | 60 | 5D Multi-Epitope Surface Clustering |
| **$`K = 6`$ Cuts** | `genetic` | 180 | 65 | 6D High-Dimensional Coordinate Scan |
| **$`K = 7`$ Cuts** | `genetic` | 200 | 70 | 7D Ultra-Deep Hyper-Volume Optimization |
| **$`K \ge 8`$ Cuts** | `genetic` | 250 | 80 | Complex Deep Lattice Cluster Optimization |
| **Low-Density Data (Discrete)** | `systematic` | *N/A (Auto-collapsed)* | *N/A (Auto-collapsed)* | Rigid Quantile Step Filtering / Cell Floor Defense |
| **Validation / Bootstrap Loop** | *Context Snapped* | 10 (Streamlined default) | 2 (Streamlined default) | Accelerated Resampling Stability Assessment |

## OptSurvCutR 0.9.8 (3/6/2026)

### New Core Features & Architecture

- **Optimal Cut Number Discovery (`find_cutpoint_number`)**: Introduced
  a new core function to mathematically determine the optimal number of
  cut-points (0 to 4) before running the genetic search. It supports
  model selection via Information Criteria (`AIC`, `AICc`, and `BIC`) to
  balance accuracy against complexity, preventing overfitting to
  sample-specific noise.
- **4-Tier Validation Stability System**: Completely overhauled the
  output of
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md).
  The package now automatically evaluates both the **Precision**
  (Maximum CI Width) and **Validity** (Interval Overlap) of bootstrapped
  thresholds, categorising them into four distinct tiers: OPTIMAL (Tier
  1), DISTINCT (Tier 2), CAUTION (Tier 3), and UNSTABLE (Tier 4).

#### Stability Tier Definitions

| Stability Tier | Criteria | Clinical Interpretation |
|:---|:---|:---|
| **Tier 1 (OPTIMAL)** | Narrow CI Width (\< 30%) AND No Overlap | High reliability; thresholds are highly consistent across different samples. |
| **Tier 2 (DISTINCT)** | Zero Interval Overlap (Regardless of CI Width) | The threshold value may vary, but the subpopulations remain fundamentally separate. |
| **Tier 3 (CAUTION)** | Moderate CI Width (30%–60%) WITH Overlap | Moderate instability; thresholds are likely real but sensitive to outliers. |
| **Tier 4 (UNSTABLE)** | High CI Width (\> 60%) WITH Overlap | High instability; the model is likely overfitting to sample-specific noise. |

- **Modernised Console UI (`cli` Integration)**: Replaced standard
  [`cat()`](https://rdrr.io/r/base/cat.html) and
  [`print()`](https://rdrr.io/r/base/print.html) outputs with the `cli`
  package. The package now features a highly readable, colour-coded,
  text-wrapping console UI for summaries, diagnostics, and tiered
  alerts, safely bypassing base R’s 1,000-character warning limits.

### Bare-Metal Low-Level Core Optimisation

- **Rcpp Integer-Assignment Factory (`cpp_get_group_assignments`)**:
  Migrated the deepest internal bottlenecks from high-level R loops to a
  native, compiled C++ backend via `Rcpp`. Instead of relying on R’s
  single-threaded vector copying, a high-speed C++ loop now evaluates
  continuous variables against candidate thresholds and constructs group
  partitions directly at the hardware memory layer.
- **Elimination of `findInterval` and `as.factor` Bottlenecks**:
  Replaced the native R group-allocation pipeline within the exhaustive
  grid search (`.get_stat`) and the evolutionary objective function
  (`.obj`). The engine now bypasses R’s high-level object allocation
  layers entirely, cutting memory overhead significantly during heavy
  permutation loops and genetic generation iterations.
- **Vectorised [`tabulate()`](https://rdrr.io/r/base/tabulate.html)
  Constraint Evaluation**: Replaced the expensive R
  [`table()`](https://rdrr.io/r/base/table.html) and
  `any(nlevels() < num_cuts)` validation checks with fast C++ matrix
  array metrics and vectorised
  [`tabulate()`](https://rdrr.io/r/base/tabulate.html) calls. Subgroup
  sizes are now verified instantly using raw memory pointers before any
  survival formula is compiled, saving thousands of microsecond
  allocations per evaluation cycle.
- **Hardened S3 Attribute Anchoring**: Implemented a protected S3
  structural binding using the `structure(..., class = "factor")`
  paradigm to map the C++ integer vectors into formula-compliant
  categorical fields. This prevents R’s internal object-replacement
  methods from stripping away structural attributes when dropped into
  the data pipeline, ensuring seamless, crash-free interoperability with
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html) and
  [`survival::survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html).

### Algorithmic Optimisation & Tie-Handling (Sparse Data Data-Paths)

- **Defensive Boundary Safety Buffering (The “Tie-Wedge”)**: Introduced
  an internal 2-patient mathematical tolerance window (`nmin - 2`) into
  both the systematic and genetic validation layers. This architectural
  guardrail resolves a structural vulnerability in global optimisation
  routines when evaluating heavily tied or highly zero-skewed continuous
  vectors. When a large percentage of observations share identical
  values (at the lower limit of detection or zero baselines), rigid
  group-size partitioning constraints become mathematically impossible
  to fulfil. The introduction of this tolerance window prevents strict
  gridlock rejections (`-Inf` / `NA`), allowing the core optimisation
  loops to safely navigate around data ties and locate valid optimal
  solutions.
- **Fail-Safe Compilation Guardrails**: Added an environment tracking
  verification step to identify compilation or loading failures in the
  compiled C++ shared libraries at runtime. If an outdated server or
  obscure OS distribution blocks the binary namespace, the package
  prints a clear warning and automatically triggers a graceful
  degradation path (`use_cpp = FALSE`), falling back to the transparent,
  native R loop structures without throwing a fatal crash.

### Major Performance Improvements & Memory Protections

- **C-Level Interval Math**: Replaced
  [`cut()`](https://rdrr.io/r/base/cut.html) with
  `findInterval(..., left.open = TRUE)` within the fallback mathematical
  engines. This eliminates expensive string manipulations, making
  systematic grid searches faster while perfectly preserving the
  mathematical boundaries of the groups.
- **Genetic Algorithm Memoization**: Implemented a hash-based evaluation
  cache (`eval_cache`) inside the genetic algorithm wrapper. The
  algorithm now remembers previously evaluated cut-points and bypasses
  the [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  model entirely for redundant guesses, drastically cutting computation
  time on large datasets with many generations.
- **Pre-Allocation for Genetic Search**: Shifted text manipulation and
  formula generation outside the genetic algorithm loop. The algorithm
  now uses a pre-allocated `data.frame` template, eliminating thousands
  of redundant memory allocations.
- **OS-Optimized Bootstrapping**: Upgraded the parallel processing
  backend in
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md).
  On Unix-based systems (Mac/Linux), the package now dynamically
  switches to **FORK clusters** (shared memory) instead of PSOCK. This
  drops data-transfer overhead to near-zero and drastically speeds up
  bootstrap validation.
- **Environment Memory-Isolation via Explicit Binding**: Overhauled the
  internal execution layers to pass explicitly extracted atomic vector
  streams instead of raw, nested evaluations. This isolates the parent R
  session environments from leaking into parallel background threads,
  cutting down cumulative memory inflation during high-core execution.

### Advanced Visualisation & Reporting (New Features)

- **Unified S3 Plotting Router**: Overhauled the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  ([`plot.find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/plot.find_cutpoint.md)).
  It now supports full `...` argument passthrough to `survminer`
  functions for deep customization, and includes a `return_data = TRUE`
  “escape hatch” to extract the raw, stratified plotting data.
- **Schoenfeld Diagnostic Plots & 2-Tier Alert**: Added
  `type = "diagnostic"` to automatically evaluate the proportional
  hazards assumption via
  [`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
  and plot the residuals. The
  [`summary()`](https://rdrr.io/r/base/summary.html) function now
  includes a 2-Tier diagnostic alert to warn users if predictive power
  shifts significantly over time, gracefully handling singular matrix
  edge cases.
- **Dashboard View**: Added `type = "all"` to generate a comprehensive,
  stacked composite plot (powered by `patchwork`) showing both the
  predictor distribution and the resulting survival outcome curve in a
  single clinical snapshot.
- **Interactive Web Widgets**: Introduced `optsurv_interactive()`, a
  wrapper function that converts any static OptSurvCutR plot into an
  interactive HTML widget via `plotly` (ideal for Vignettes and
  RMarkdown).
- **Clinical Aesthetics**: Implemented
  [`theme_optsurv()`](https://paytonyau.github.io/OptSurvCutR/reference/theme_optsurv.md)
  and defaulted to colorblind-safe palettes (e.g., “nejm”) to enforce a
  unified, publication-ready aesthetic across all outputs.

### Bug Fixes & Edge Cases (`find_cutpoint`)

- **Dynamic Log-Likelihood Extraction**: Fixed a bug where the `p_value`
  criterion would fail (or skip tests) when no covariates were provided.
  The null model log-likelihood is now extracted dynamically regardless
  of model length.
- **Degenerate Data Handling**: Added a safeguard to check for
  zero-variance in the `time` column, ensuring the function safely
  returns `NA` rather than crashing on pathological datasets.
- **Dropped Coefficients**: Added a check to safely return `-Inf`
  instead of `NA` if
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  completely drops a coefficient (e.g., due to extreme collinearity)
  during a `hazard_ratio` search.

### Structural Adjustments & Edge-Case Vulnerabilities (`find_cutpoint_number`)

- **Unadjusted Profile Likelihood Fallback Engine**: Resolved a critical
  evolutionary calculation bottleneck inside `.obj()` when running
  unadjusted genetic searches (`method = "genetic"`). By adding an
  implicit zero-covariate structural detector, the engine now
  dynamically bypasses manual parameter `init` beta matrix arrays when
  no adjusters are present. This prevents
  `survival::coxph(..., iter.max = 0)` from crashing on bare categorical
  split factors, ensuring unadjusted model selection tables calculate 1
  and 2 cuts perfectly.
- **Harden Nmin Evaluation Sequence**: Re-ordered the parameter
  evaluation pipeline to calculate absolute patient sample splits
  (`nmin_abs`) prior to executing total data headroom matrix checks.
  This successfully stops proportional constraints (e.g., `nmin = 0.25`)
  from silently bypassing sample capacity filters.
- **Boundary Tie Headroom Margin**: Built a 2-patient safety headroom
  filter into the information criteria data-capacity logic. This
  automatically intercepts brittle or impossible search spaces (such as
  requesting 3 cuts on 580 patients with a strict 25% allocation limit)
  and down-caps or exits gracefully before triggering genetic algorithm
  failures.

### Stability & Parallel Processing (`validate_cutpoint`)

- **CRAN-Safe Core Detector**: Built an environment check
  (`Sys.getenv("_R_CHECK_LIMIT_CORES_")`) into the core architecture.
  The package now detects when it is running on CRAN’s testing servers
  and automatically throttles itself to 2 cores to prevent rejection,
  while allowing real users to utilise maximum CPU power.
- **Parallel Variable-Name Scope Isolation**: Resolved a critical
  parallel clustering crash
  (`argument "predictor" is missing, with no default`) that occurred
  inside `.run_permutations()` during high-throughput bootstrap
  resampling. The signature now maps explicit variable text strings
  alongside the split memory streams, ensuring independent worker
  processes (`%dopar%`) can parse formulas cleanly across separate CPU
  nodes.
- **Data-Stream Alignment Fix for Scoped Models**: Fixed a dataset scope
  conflict (`argument "data" is missing, with no default`) inside the
  parallel `foreach` environment. By binding the structural data
  components directly within the function parameter signature rather
  than relying on global environmental evaluations, parallel worker
  environments successfully execute deep iterations across any operating
  system without variable scope dropouts.
- **Renamed Predictor Collision Handling**: Fixed a tracking bug inside
  `summary.find_cutpoint_number_result` where the function printed
  `"Unknown"` for the main predictor variable under specific
  covariate-adjusted models. The S3 summary framework now tracks
  original column metadata across all renaming phases to provide
  pristine diagnostic dashboards for clinical research reporting.
- **Parallel Closure Fix (`...` arguments)**: Fixed a critical bug where
  running `validate_cutpoint` with `n_cores > 1` would instantly fail if
  extra arguments (like `max.generations` or `pop.size`) were passed to
  the `foreach` loop.
- **Console Spam Prevention**: Bootstrapping with replacement naturally
  creates tied survival times, which can trigger expected model
  convergence warnings. Parallel model fits are now wrapped in
  [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) to prevent
  these benign warnings from flooding the user’s console and breaking
  the progress bar.
- **Namespace Safety**: Refactored the parallel execution block to
  reliably call the exported
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  wrapper. This guarantees that parallel worker nodes on any operating
  system can successfully locate the internal evaluation functions
  without namespace crashes.
- **S3 Summary Parameter Printing Fix**: Resolved an unescaped string
  and operator evaluation syntax error within the
  [`cli::cli_bullets`](https://cli.r-lib.org/reference/cli_bullets.html)
  parameter tracking block inside `summary.find_cutpoint_number_result`.
  All named bullet outputs are now tightly mapped to explicitly bounded,
  single-quoted parameters (`"*" =`), guaranteeing clean console
  rendering.

## OptSurvCutR 0.2.1 (2026-04-20)

### CRAN Compliance & Quality of Life Improvements

This patch release addresses CRAN reviewer feedback and polishes the
package’s console behavior and visual branding.

- **Enhanced Console Control:** Replaced all informational
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`cat()`](https://rdrr.io/r/base/cat.html) statements with
  [`message()`](https://rdrr.io/r/base/message.html) in the core
  functions
  ([`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md),
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md),
  and
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)).
  Users can now easily silence progress text by wrapping functions in
  [`suppressMessages()`](https://rdrr.io/r/base/message.html).
- **Cleaned Return Behaviors:** Removed forced
  [`print()`](https://rdrr.io/r/base/print.html) calls from the end of
  the main S3 calculation functions. Results now return silently when
  assigned to a variable, while preserving formatted output when called
  directly.
- **DESCRIPTION File Formatting:** Removed single quotes around function
  and package names in the `DESCRIPTION` file to satisfy CRAN automated
  parsers.
- **Official Branding:** Added the official `OptSurvCutR` hex logo! The
  logo is now bundled within the package (`man/figures/logo.png`) and
  fully integrated into the GitHub README and `pkgdown` website
  configuration.

## OptSurvCutR 0.2.0 (2026-04-09)

### DOCUMENTATION & STANDARDS

- Initial CRAN submission.
- Updated the package strucure (minor) to meet the submission
  requirement.

## OptSurvCutR 0.1.9.3 (2026-03-22)

### DOCUMENTATION & STANDARDS

**Graphical Abstract:** Added the new graphical abstract to the `image/`
folder for inclusion in the README.md. ![OptSurvCutR Graphical
Abstract](image/graphical_abstract.jpg)

- **Instruction Page (README.md):** Updated the main package repository
  instructions.

### COMMUNITY & OUTREACH

**Presentations:** Presented at FOSDEM 2026 in the Bioinformatics &
Computational Biology DevRoom in Brussels, Belgium. [View lightning talk
details](https://fosdem.org/2026/schedule/event/9JQJ9M-bioinformatics_lighthning_talks/)

**Presentations:** Presented at R!SK 2026 (February 18–19, 2026). This
100% online conference focuses on evaluating, measuring, and mitigating
risk across diverse industries including healthcare, finance, and
insurance. The event featured deep content sessions and live Q&A
interactions. [View presentation
abstract](https://rconsortium.github.io/Risk_website/Abstracts.html#payton-yau)

## OptSurvCutR 0.1.9.2 (2025-11-25)

### ROPENSCI REVIEW UPDATES (ITERATION 2)

#### BREAKING CHANGES

- **Renamed Functions:** `plot_optimization_curve()`
  $`\to`$[`plot_optimisation_curve()`](https://paytonyau.github.io/OptSurvCutR/reference/plot_optimisation_curve.md)
  (British English) and `plot_schoenfeld()`
  $`\to`$[`plot_cutpoint_residuals()`](https://paytonyau.github.io/OptSurvCutR/reference/plot_cutpoint_residuals.md)
  (Namespace conflict).
- **Renamed Argument:** `maxiter` is now `max.generations` in genetic
  algorithm functions to match `rgenoud`.

#### IMPROVEMENTS

- **Explicit Parameters:** Exposed `pop.size` and `max.generations` as
  explicit arguments in
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  and
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md).
- **Consistent Logic:** All functions now use
  [`floor()`](https://rdrr.io/r/base/Round.html) for `nmin` calculation
  to ensure identical sample size handling.
- **Better Validation:**
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
  now accepts `nmin` as a proportion and issues an informative message
  when using the default 90% buffer (Standard G2.10).
- **Internal Cleanup:** Removed export of helper operator `%||%`.

#### DOCUMENTATION & STANDARDS

- **British English:** Standardised all functions and documentation
  (e.g., “optimise”, “summarise”).
- **Standards Compliance:** Audited and updated all `@srrstats` tags;
  moved non-applicable standards to `R/OptSurvCutR-package.R`.
- **Testing:** Added tests for new genetic parameters and updated
  snapshot tests for S3 methods.

------------------------------------------------------------------------

## OptSurvCutR 0.1.9.1 (2025-11-20)

### POST-ROPENSCI INITIAL SUBMISSION

- Renamed `plot_diagnostics()` → `plot_schoenfeld()` to avoid name clash
  with another package.
- Stabilised and accepted all snapshot tests for pathological
  `find_cutpoint_number_result` S3 methods (now use real criterion
  `"BIC"`).
- Cleaned DESCRIPTION: removed unused Suggests (`coxphf`, `patchwork`,
  `pROC`, `tibble`, `srr`) and Remotes (`ropensci-review-tools/srr`).
- Test suite now passes with zero failures/warnings on R CMD check.

------------------------------------------------------------------------

## OptSurvCutR 0.1.9

### REFACTORING & MAINTENANCE

- **Simplified Parallel Interface:** The `use_parallel` argument was
  removed from
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md).
  Parallel execution is now controlled **exclusively** by the `n_cores`
  argument (`n_cores = 1` for sequential, `n_cores > 1` for parallel).
- **Consolidated Validation Logic:** A new internal helper
  `.validate_event_column()` was created in `R/utils-helpers.R` to
  centralise the logic for checking that an event column is numeric and
  contains only 0s and 1s. Both
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  and the internal helper `.validate_data_conditions()` (also in
  `R/utils-helpers.R`) were updated to call this function, eliminating
  duplicated code.
- **Simplified Parallel Code:** In
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md),
  the brittle `functions_to_export` block was removed. The function now
  correctly relies on the `.packages = "OptSurvCutR"` argument in the
  `foreach` loop, making it easier to maintain.

### BUG FIXES

- **Fixed Parallel Reproducibility:** A bug in
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
  that prevented true reproducibility for parallel runs (`n_cores > 1`)
  was fixed. The function now checks for and registers the
  [doRNG](https://github.com/emilioluissaenzguillen/doRNG) package when
  a `seed` is provided, ensuring results are identical regardless of the
  number of cores used. The incorrect `set.seed(i)` call inside the
  `foreach` loop was removed.
- **Improved S3 Method Robustness:** The
  [`summary.find_cutpoint_number_result()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  S3 method was updated to robustly handle `NULL` values in the
  `object$parameters` list (e.g., `method = NULL`). This prevents
  potential errors if a result object is created manually or improperly
  and aligns its behaviour with the
  [`print()`](https://rdrr.io/r/base/print.html) method.

### DOCUMENTATION

- **Clarified Log-Rank Test:** The documentation for
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  was updated to clarify that when `covariates` are provided, the
  `"logrank"` criterion is automatically generalized to the more
  appropriate **Cox score test**.

### IMPROVEMENTS

- **Test Coverage:** Increased test coverage to **88.5%** with new unit
  tests for parallel reproducibility, S3 method edge cases, and internal
  validation helpers.
- Planned: A ROpenSci submission is planned for review.
- Planned: A JOSS submission is planned post-rOpenSci review.

------------------------------------------------------------------------

## OptSurvCutR 0.1.8

### CRITICAL BUG FIXES

- **Parallel Validation:** Fixed a critical bug where
  `validate_cutpoint(use_parallel = TRUE)` would fail. This was caused
  by helper functions not being exported to the parallel workers.
- **Event Column Validation:** Fixed a major bug where non-numeric
  `outcome_event` data (e.g., `"0:LIVING"`) caused silent failures. The
  functions now “fail-fast” with an error, validating the column is
  numeric with only 0s and 1s.
- **Genetic Algorithm Reporting:** Fixed a bug where
  `method = "genetic"` would report a failure (`-Inf`) as a successful
  result. The functions now correctly check for the failure signal and
  return `NA`.

------------------------------------------------------------------------

### NEW FEATURES

- **G3.1a Compliance — Schoenfeld Residual Diagnostics**: Added
  `plot_diagnostics()` to generate publication-ready Schoenfeld residual
  plots, assessing the proportional hazards (PH) assumption.

------------------------------------------------------------------------

### IMPROVEMENTS

- **Code Quality (Refactoring):** Refactored
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  to resolve “high cyclomatic complexity” `NOTE`s by moving validation
  logic to internal helper functions.
- **Test Coverage:** Increased test coverage by adding new unit tests
  for error conditions, edge cases, and validation logic.
- **rOpenSci/`pkgcheck` Compliance:**
  - Added examples to all exported functions.
  - Added `@srrstats` tags to satisfy multi-directory requirements.
- **`R CMD check` Compliance:** Addressed all `ERROR`s, `WARNING`s, and
  `NOTE`s:
  - Added `inst/WORDLIST`, `_PACKAGE` documentation, and fixed empty
    `Rd` sections.
  - Added `codemeta.json` to `.Rbuildignore`.
- **Code Style:**
  - Reformatted code to follow Tidyverse Style Guide principles (e.g.,
    `<-` for assignment, `"` for strings, consistent spacing around
    operators).
  - Strictly enforced an 80-character maximum line length for all code,
    comments, and documentation.
  - Made Roxygen documentation, inline comments, and user-facing
    messages more concise.
  - Replaced [`sapply()`](https://rdrr.io/r/base/lapply.html) with
    type-safe [`vapply()`](https://rdrr.io/r/base/lapply.html) in S3
    methods to prevent potential bugs.
- **Internal Data:** The `crc_virome` dataset’s `status` column was
  corrected to be numeric 0/1 to match documentation and new validation
  rules.

------------------------------------------------------------------------

## OptSurvCutR 0.1.7

### NEW FEATURES

- **Covariate Adjustment Added**: Both
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  (for `method = "systematic"` and `method = "genetic"`) and
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  now support covariate adjustment via the `covariates` argument. This
  allows finding optimal cut-points and determining the optimal number
  of groups while accounting for potential confounders, providing a more
  robust assessment of a biomarker’s independent prognostic value.

### IMPROVEMENTS

- **Major Performance Optimisation**: Significantly optimised
  performance in
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  for `criterion = "hazard_ratio"` and `criterion = "p_value"` (both
  systematic and genetic methods). This was achieved by removing
  computationally expensive
  [`summary()`](https://rdrr.io/r/base/summary.html) calls, extracting
  coefficients and statistics directly from model objects, and using a
  fast, manual Likelihood Ratio Test for p-value calculation.
- **Simplified Parallelism**: Removed internal parallel processing
  (`use_parallel` argument) from `find_cutpoint` and
  `find_cutpoint_number` systematic search. This prevents potential
  issues with nested parallel calls and relies on standard external
  parallelisation approaches, such as the explicit parallel loop within
  `validate_cutpoint`.
- **Test Suite Overhaul**: Revamped the entire `testthat` suite for
  improved reliability and robustness. Replaced brittle
  `expect_snapshot()` tests with more stable checks like
  `expect_output()`, `expect_s3_class()`, and specific value
  comparisons. Corrected logic for error and warning expectations and
  improved mocking for dependency checks. Code coverage increased
  significantly (for example, to ~86%).
- **User-Friendly `rgenoud` Check**: Added clear, informative error
  messages using `cli` in
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  and
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  that trigger immediately if `method = "genetic"` is requested but the
  suggested `rgenoud` package is not installed, guiding the user on how
  to install it.
- **Dependency Management**: Moved optional dependencies `broom` (for
  plotting) and `withr` (for testing) from `Imports` to `Suggests` in
  the `DESCRIPTION` file, making the core package installation lighter.
  Removed unused `magrittr` import.
- **Simplified Evidence Labels**: Renamed evidence labels in
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  results for brevity (for example, “Substantial support” -\>
  “Substantial”).
- **Code Maintainability**: Centralised all
  [`utils::globalVariables`](https://rdrr.io/r/utils/globalVariables.html)
  definitions into `globals.R` to resolve R CMD check `NOTEs` and
  improve clarity. Removed redundant code.

### BUG FIXES

- **quiet = TRUE Message Fix**: Fixed a critical bug where
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  and
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  would fail silently (no console message) when `quiet = TRUE` was set.
  Failure messages are now always printed to the console via internal
  helpers, regardless of the `quiet` setting, improving user feedback on
  errors.
- **Genetic Algorithm Edge Cases**: Added more robust input validation
  and handling within the internal `.run_genetic_search()` function to
  gracefully manage edge cases like insufficient data variability or
  non-finite predictor ranges, preventing downstream errors and cryptic
  warnings.
- **Genetic Algorithm Monitor**: Fixed the internal monitoring function
  used by the genetic algorithm
  ([`rgenoud::genoud`](https://rdrr.io/pkg/rgenoud/man/genoud.html)) to
  correctly respect the `print.level` argument (controlled indirectly
  via user functions), ensuring progress updates are displayed or
  suppressed as intended.
- **NAMESPACE Fix (stats::)**: Resolved `NAMESPACE` errors and related
  test failures by adding explicit `stats::` calls where needed (for
  example, for
  [`stats::quantile`](https://rdrr.io/r/stats/quantile.html),
  [`stats::sd`](https://rdrr.io/r/stats/sd.html),
  [`stats::pchisq`](https://rdrr.io/r/stats/Chisquare.html)) and
  ensuring correct regeneration of the `NAMESPACE` file via
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).
- **foreach NOTE Fix**: Resolved R CMD check `NOTE` regarding “no
  visible binding for global”.

### IMPROVEMENTS

- Planned: A ROpenSci submission is planned for review.
- Planned: A JOSS submission is planned post-rOpenSci review.

------------------------------------------------------------------------

## OptSurvCutR 0.1.6

### NEW FEATURES

- Added a vignette demonstrating the use of
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  and
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
  with TCGA virome data (for example, Alphapapillomavirus as a
  predictor), guiding users through cut-point optimisation and stability
  assessment for survival analysis
  \[<https://github.com/paytonyau/OptSurvCutR/commit/aa41ca3cb3ff7fdff4cbf6cf8d5de5e4494d3500>\].
- Introduced comprehensive unit tests using `testthat`, covering core
  functions
  ([`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md),
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md),
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md))
  and edge cases like missing data or small sample sizes, with code
  coverage reporting via `covr` to ensure reliability (\>80% coverage)
  \[<https://github.com/paytonyau/OptSurvCutR/commit/6579072b087b880b448e101d8eef37f8c4fa5550>\].

### IMPROVEMENTS

- Optimised the genetic algorithm in
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  by implementing adaptive `pop.size` (for example, 50 for
  `num_cuts = 1`) and `max.generations` (for example, 75), reducing
  runtime by 20–50% for survival datasets while maintaining accuracy for
  optimal cut-point selection
  \[<https://github.com/paytonyau/OptSurvCutR/commit/e0cec20c2a72e48b39e1833742e8d7d829621f39>\].
- Enhanced error messages in
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
  to provide specific feedback on bootstrap validation failures, such as
  insufficient sample sizes or non-converging `coxph` models, improving
  user debugging experience
  \[<https://github.com/paytonyau/OptSurvCutR/commit/92ef362cec36ee9ae976f12a61b65db51cc79d94>\].
- Added a `pkgdown` GitHub Action to automatically build a package
  website, improving documentation accessibility, and updated README
  with badges for build status and code coverage to signal package
  reliability
  \[<https://github.com/paytonyau/OptSurvCutR/commit/38712b8aa3cb919556bdf0e9cba6ca27fda10a60>\].
- Updated `DESCRIPTION` with corrected URLs, dependency versions, and
  regenerated Rd files for consistent documentation across all functions
  \[<https://github.com/paytonyau/OptSurvCutR/commit/d85d4a14a1a10fb7f68537ea64ef4007c3be465a>\].

### BUG FIXES

- Fixed NA handling in
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  to robustly process survival datasets with missing predictor values,
  preventing errors in `coxph` or `survdiff` model fitting for
  real-world data like TCGA virome datasets
  \[<https://github.com/paytonyau/OptSurvCutR/commit/059a288363c9b7b30272301b709378cb58d76d2b>\].

------------------------------------------------------------------------

## OptSurvCutR 0.1.5

### IMPROVEMENTS

- Updated core functions (for example,
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md),
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md))
  with improved numerical stability and accuracy for survival model
  fitting, particularly for genetic algorithm convergence in
  high-dimensional predictors
  \[<https://github.com/paytonyau/OptSurvCutR/commit/9338d0479c34069933562cb6b360428aad9dd6fc>\].

### BUG FIXES

- Fixed bugs in script handling and input validation, improving
  reliability for edge cases like small datasets or constant predictors
  in survival analysis
  \[<https://github.com/paytonyau/OptSurvCutR/commit/1de3e1048dcaa33fa8cbb6eab0fa8d89ec5d134c>\].
- Reverted prior bug fixes to prevent potential regressions, ensuring
  stable behaviour in
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
  during bootstrap validation runs
  \[<https://github.com/paytonyau/OptSurvCutR/commit/d65fa328634a480a27bde6f64f8d615205e12225>\].

------------------------------------------------------------------------

## OptSurvCutR 0.1.0

### NEW FEATURES

- Initial release of `OptSurvCutR` for optimising cut-points in survival
  analysis.
- Added
  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  to identify optimal cut-points for continuous predictors using
  systematic or genetic algorithms (via `rgenoud`) with log-rank,
  p-value, or hazard-ratio criteria.
- Added
  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md)
  to select the optimal number of cut-points using AIC, AICc, or BIC.
- Added
  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
  for bootstrap-based stability assessment of cut-points.
- Supports survival models via
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html) and
  [`survival::survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html).
- Includes example usage with simulated churn data, adaptable to TCGA
  virome datasets (for example, Alphapapillomavirus).
