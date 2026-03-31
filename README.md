## Project Overview

`doubleShrinker` is an R package accompanying "Empirical Bayes Double Shrinkage for Combining
Biased and Unbiased Causal Estimates" by Rosenman, Dominici, and Miratrix. 

The package implements Empirical Bayes double shrinkers for combining biased and unbiased estimates of the same quantity. The package also provides competitor estimators and tools for computing Empirical Bayes  confidence intervals. Simulation scripts from the manuscript, comparing the double shrinkers' performance against competitors, live in `scripts/`.


## Installation
  
```r
# install.packages("remotes")
remotes::install_github("etrrosenman/doubleShrinkerCode")
```

## Running Simulations

Simulation scripts are in `scripts/`. Use `library(doubleShrinker)` after first installing the package. 

```r
# From the project root
Rscript scripts/MSE_Sims.R    # MSE comparison simulation
Rscript scripts/CI_Sims.R     # Confidence interval coverage simulation
```

Additional packages needed by the scripts: `mvtnorm`, `rootSolve`, `here`, `parallel`

## Package Structure

### `R/` source files
- **`estimators.R`** — Core double shrinker: `shrinker()` and the four `eb.*` hyperparameter estimators
- **`competitors.R`** — Competitor estimators: `precistion.wtd()`, `kappa.1()`, `kappa.2()`, `delta.1()`, `delta.2()`
- **`data_helpers.R`** — Data preparation: `subgroupDefs_noSplit()`, `getParamEstimates()`, `stratifiedBootstrapSample()`
- **`ci.R`** — Confidence intervals: `ebCiFunc()`, `coverage()`
- **`utils.R`** — Internal helpers (not exported): `URE`, `UREgrad`, `ebLik`, `ebLikGrad`

### `data/`
- **`synthetic_data.os.RData`** — Synthetic observational dataset
- **`synthetic_data.rct.RData`** — Synthetic experimental dataset

### `scripts/`
- **`MSE_Sims.R`** — Bootstrap simulation comparing estimators by MSE, bias, and variance
- **`CI_Sims.R`** — Bootstrap simulation evaluating CI coverage and length

## Key Functions

**Data preparation:**
- `getParamEstimates(rctSample, obsSample, outcome, propScoreAdjust)` — Estimates per-stratum RCT and OS treatment effects and variances; supports IPW propensity score adjustment
- `stratifiedBootstrapSample(sampleIds, delta)` — Fast stratified bootstrap using pre-stored IDs and `dqsample`

**Hyperparameter estimators (all call `shrinker` internally):**
- `eb.mm1()` — Method of moments 1 (moment conditions on RCT estimates)
- `eb.mm2()` — Method of moments 2 (alternative moment conditions)
- `eb.mle()` — Maximum likelihood via `constrOptim` + `L-BFGS-B`
- `eb.ure()` — Unbiased risk estimation

All four accept `lowerBound=TRUE` and `returnCIComps=TRUE` (returns list with `shrinker`, `lambda`, `a`, `c`, `etaSqLwr`, `gammaSqLwr`).

**Competitor estimators:**
- `precistion.wtd()` — Precision-weighted combination
- `kappa.1()`, `kappa.2()` — Rosenman et al. (2023) scalar/diagonal shrinkage
- `delta.1()`, `delta.2()` — Green et al. (2005) James-Stein-style shrinkage

**Confidence intervals:**
- `ebCiFunc(ebData, ed, alpha)` — EB CI half-widths using `cva()` from the `ebci` package
- `coverage()` — Checks whether true values fall within constructed intervals

## Data Conventions

- Outcome variable: `"Outcome.CHD"` (coronary heart disease, binary)
- Subgrouping variables: `"CVD"`, `"AGER"`, `"LANGLEYSCAT"`
- Treatment indicator: `Test` (1 = treated, 0 = control)
- The `ed` data frame passed to estimators has columns: `rctEst`, `rctVar`, `obsEst`, `obsVar`, `stratum`

