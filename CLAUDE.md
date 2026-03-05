# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R research codebase implementing "double shrinker" empirical Bayes estimators that combine RCT (Randomized Controlled Trial) and observational study data for subgroup causal effect estimation. The work evaluates these estimators against competitors via simulation studies using real data as the population.

## Running the Code

Run the simulation scripts from the project root (the `.Rproj` file establishes the working directory via `here()`):

```r
# In R or RStudio
source("MSE_Sims.R")   # MSE comparison simulation
source("CI_Sims.R")    # Confidence interval coverage simulation
```

Or from the command line:
```bash
Rscript MSE_Sims.R
Rscript CI_Sims.R
```

Required R packages: `mvtnorm`, `plyr`, `dqrng`, `rootSolve`, `here`, `parallel`, `ebci`, `tibble`

## Architecture

### File Structure
- **`Helper Functions.R`** — All core functions; sourced by both simulation scripts
- **`MSE_Sims.R`** — Bootstrap simulation comparing estimators by MSE, bias, and variance
- **`CI_Sims.R`** — Bootstrap simulation evaluating confidence interval coverage and length
- **`synthetic_data.os.RData`** / **`synthetic_data.rct.RData`** — Synthetic population datasets (observational study and RCT)
- **`data.os.RData`** / **`data.rct.RData`** — Real datasets

### Data Flow

Both simulation scripts follow the same pattern:
1. Load synthetic data as a population
2. Loop over subgroup variable combinations (7 total, from 1-variable to 3-variable crosses)
3. For each subgroup configuration, run `numSamples=1000` stratified bootstrap iterations in parallel (`mclapply`)
4. Each iteration: resample from population → estimate parameters → apply all estimators → record results

### Key Functions in `Helper Functions.R`

**Data preparation:**
- `subgroupDefs_noSplit()` — Tags units with subgroup labels across both datasets; computes "true" treatment effects from RCT data
- `getParamEstimates()` — Estimates per-stratum RCT and observational treatment effects and variances (supports IPW propensity score adjustment for the OS)
- `stratifiedBootstrapSample()` — Fast stratified bootstrap using pre-stored IDs and `dqsample`

**Core double shrinker:**
- `shrinker(gammaSq, etaSq, ed)` — The fundamental estimator. `gammaSq` is the variance of the bias (RCT - obs discrepancy), `etaSq` is the signal variance (true treatment effect variance). Combines RCT and OS estimates via weights `lambda` and `a`.

**Hyperparameter estimation methods (all call `shrinker` internally):**
- `eb.mm1()` — Method of moments estimator 1 (uses cross terms of rctEst and obsEst)
- `eb.mm2()` — Method of moments estimator 2 (alternative moment conditions)
- `eb.mle()` — Maximum likelihood via `constrOptim` + `L-BFGS-B` with `ebLik`/`ebLikGrad`
- `eb.ure()` — Unbiased risk estimation via minimizing `URE`/`UREgrad`

All four accept `lowerBound=TRUE` to prevent near-zero hyperparameter estimates, and `returnCIComps=TRUE` to return components needed for CI construction.

**Competitor estimators:**
- `precistion.wtd()` — Precision-weighted combination
- `kappa.1()`, `kappa.2()` — Rosenman et al. (2023) scalar/diagonal shrinkage
- `delta.1()`, `delta.2()` — Green et al. (2005) James-Stein-style shrinkage

**Confidence intervals:**
- `ebCiFunc()` — Constructs empirical Bayes CI half-widths using the `cva()` function from the `ebci` package
- `coverage()` — Checks whether true values fall within constructed intervals

### Outcome Variable

The outcome is `"Outcome.CHD"` (coronary heart disease). Subgrouping variables are `"CVD"`, `"AGER"`, and `"LANGLEYSCAT"`. The `Test` column indicates treatment assignment (1 = treated, 0 = control).
