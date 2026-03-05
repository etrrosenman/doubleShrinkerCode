set.seed(123)

# ── helper to build minimal mock datasets ─────────────────────────────────────

make_rct <- function(n_per_cell = 40, outcome = "Y") {
  d <- data.frame(
    subgroup = factor(rep(c("A", "B"), each = n_per_cell * 2)),
    Test     = rep(c(0L, 1L, 0L, 1L), each = n_per_cell),
    stringsAsFactors = FALSE
  )
  # outcome rates differ by stratum and arm
  rates <- c(A_ctl = 0.20, A_trt = 0.30, B_ctl = 0.10, B_trt = 0.15)
  d[[outcome]] <- c(
    rbinom(n_per_cell, 1, rates["A_ctl"]),
    rbinom(n_per_cell, 1, rates["A_trt"]),
    rbinom(n_per_cell, 1, rates["B_ctl"]),
    rbinom(n_per_cell, 1, rates["B_trt"])
  )
  d
}

make_obs <- function(n_per_cell = 200, outcome = "Y") {
  d <- data.frame(
    subgroup  = factor(rep(c("A", "B"), each = n_per_cell * 2)),
    Test      = rep(c(0L, 1L, 0L, 1L), each = n_per_cell),
    propScore = rep(0.5, n_per_cell * 4),
    stringsAsFactors = FALSE
  )
  rates <- c(A_ctl = 0.22, A_trt = 0.28, B_ctl = 0.12, B_trt = 0.18)
  d[[outcome]] <- c(
    rbinom(n_per_cell, 1, rates["A_ctl"]),
    rbinom(n_per_cell, 1, rates["A_trt"]),
    rbinom(n_per_cell, 1, rates["B_ctl"]),
    rbinom(n_per_cell, 1, rates["B_trt"])
  )
  d
}

# ── getParamEstimates() ───────────────────────────────────────────────────────

test_that("getParamEstimates returns a data frame with correct dimensions", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)   # two strata: A and B
  expect_named(out, c("stratum", "rctEst", "rctVar", "obsEst", "obsVar"))
})

test_that("getParamEstimates variances are strictly positive", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y")
  expect_true(all(out$rctVar > 0))
  expect_true(all(out$obsVar > 0))
})

test_that("getParamEstimates estimates are finite", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y")
  expect_true(all(is.finite(out$rctEst)))
  expect_true(all(is.finite(out$obsEst)))
})

test_that("getParamEstimates without propensity adjustment gives finite results", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y", propScoreAdjust = FALSE)
  expect_true(all(is.finite(out$rctEst)))
  expect_true(all(is.finite(out$obsEst)))
  expect_true(all(out$rctVar > 0))
  expect_true(all(out$obsVar > 0))
})

test_that("getParamEstimates stratum order matches sorted subgroup levels", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y")
  expect_equal(as.character(out$stratum), c("A", "B"))
})

# ── subgroupDefs_noSplit() ────────────────────────────────────────────────────

make_base_data <- function(n = 200) {
  data.frame(
    X1      = sample(c(0L, 1L), n, replace = TRUE),
    X2      = sample(c(0L, 1L), n, replace = TRUE),
    Test    = sample(c(0L, 1L), n, replace = TRUE),
    Y       = rbinom(n, 1, 0.2),
    stringsAsFactors = FALSE
  )
}

test_that("subgroupDefs_noSplit adds subgroup column for single variable", {
  d <- make_base_data()
  out <- subgroupDefs_noSplit("X1", d, d, outcome = "Y")
  expect_true("subgroup" %in% names(out$data.os))
  expect_true("subgroup" %in% names(out$data.rct))
  expect_equal(nlevels(out$data.os$subgroup), 2L)
})

test_that("subgroupDefs_noSplit crosses multiple variables", {
  d <- make_base_data()
  out <- subgroupDefs_noSplit(c("X1", "X2"), d, d, outcome = "Y")
  # 2 x 2 = up to 4 strata
  expect_lte(nlevels(out$data.os$subgroup), 4L)
  expect_gte(nlevels(out$data.os$subgroup), 1L)
})

test_that("subgroupDefs_noSplit returns named subgroupTaus.true", {
  d <- make_base_data()
  out <- subgroupDefs_noSplit("X1", d, d, outcome = "Y")
  expect_true(is.numeric(out$subgroupTaus.true))
  expect_false(is.null(names(out$subgroupTaus.true)))
  expect_equal(length(out$subgroupTaus.true), nlevels(out$data.rct$subgroup))
})

# ── coverage() ────────────────────────────────────────────────────────────────

test_that("coverage returns TRUE when true value is inside interval", {
  expect_true(coverage(ptEst = 0.1, interval = 0.05, trueVals = 0.12))
  expect_true(coverage(ptEst = 0.1, interval = 0.05, trueVals = 0.06))
})

test_that("coverage returns FALSE when true value is outside interval", {
  expect_false(coverage(ptEst = 0.1, interval = 0.05, trueVals = 0.20))
  expect_false(coverage(ptEst = 0.1, interval = 0.05, trueVals = 0.00))
})

test_that("coverage works vectorised", {
  ptEst    <- c(0.05, 0.10)
  interval <- c(0.02, 0.02)
  trueVals <- c(0.06, 0.20)   # first covered, second not
  out <- coverage(ptEst, interval, trueVals)
  expect_equal(out, c(TRUE, FALSE))
})
