set.seed(123)

# ── helper to build minimal mock datasets ─────────────────────────────────────

make_rct <- function(n_per_cell = 40, outcome = "Y", trt_col = "Test") {
  d <- data.frame(
    subgroup = factor(rep(c("A", "B"), each = n_per_cell * 2)),
    stringsAsFactors = FALSE
  )
  d[[trt_col]] <- rep(c(0L, 1L, 0L, 1L), each = n_per_cell)
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

make_obs <- function(n_per_cell = 200, outcome = "Y", trt_col = "Test") {
  d <- data.frame(
    subgroup  = factor(rep(c("A", "B"), each = n_per_cell * 2)),
    propScore = rep(0.5, n_per_cell * 4),
    stringsAsFactors = FALSE
  )
  d[[trt_col]] <- rep(c(0L, 1L, 0L, 1L), each = n_per_cell)
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
  out <- getParamEstimates(rct, obs, outcome = "Y", treatment = "Test")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)   # two strata: A and B
  expect_named(out, c("stratum", "rctEst", "rctVar", "obsEst", "obsVar"))
})

test_that("getParamEstimates variances are strictly positive", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y", treatment = "Test")
  expect_true(all(out$rctVar > 0))
  expect_true(all(out$obsVar > 0))
})

test_that("getParamEstimates estimates are finite", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y", treatment = "Test")
  expect_true(all(is.finite(out$rctEst)))
  expect_true(all(is.finite(out$obsEst)))
})

test_that("getParamEstimates without propensity adjustment gives finite results", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y", treatment = "Test",
                           propScoreAdjust = FALSE)
  expect_true(all(is.finite(out$rctEst)))
  expect_true(all(is.finite(out$obsEst)))
  expect_true(all(out$rctVar > 0))
  expect_true(all(out$obsVar > 0))
})

test_that("getParamEstimates stratum order matches sorted subgroup levels", {
  rct <- make_rct()
  obs <- make_obs()
  out <- getParamEstimates(rct, obs, outcome = "Y", treatment = "Test")
  expect_equal(as.character(out$stratum), c("A", "B"))
})

test_that("getParamEstimates accepts a custom treatment column name", {
  rct <- make_rct(trt_col = "assigned")
  obs <- make_obs(trt_col = "assigned")
  out <- getParamEstimates(rct, obs, outcome = "Y", treatment = "assigned")
  expect_equal(nrow(out), 2L)
  expect_true(all(is.finite(out$rctEst)))
  expect_true(all(is.finite(out$obsEst)))
})

test_that("getParamEstimates results are identical regardless of treatment column name", {
  set.seed(1)
  rct_default <- make_rct(trt_col = "Test")
  obs_default <- make_obs(trt_col = "Test")

  rct_custom <- rct_default
  obs_custom <- obs_default
  names(rct_custom)[names(rct_custom) == "Test"] <- "assigned"
  names(obs_custom)[names(obs_custom) == "Test"] <- "assigned"

  out1 <- getParamEstimates(rct_default, obs_default, outcome = "Y", treatment = "Test")
  out2 <- getParamEstimates(rct_custom,  obs_custom,  outcome = "Y", treatment = "assigned")
  expect_equal(out1$rctEst, out2$rctEst)
  expect_equal(out1$obsEst, out2$obsEst)
})

# ── subgroupDefs_noSplit() ────────────────────────────────────────────────────

make_base_data <- function(n = 200, trt_col = "Test") {
  d <- data.frame(
    X1 = sample(c(0L, 1L), n, replace = TRUE),
    X2 = sample(c(0L, 1L), n, replace = TRUE),
    Y  = rbinom(n, 1, 0.2),
    stringsAsFactors = FALSE
  )
  d[[trt_col]] <- sample(c(0L, 1L), n, replace = TRUE)
  d
}

test_that("subgroupDefs_noSplit adds subgroup column for single variable", {
  d <- make_base_data()
  out <- subgroupDefs_noSplit("X1", d, d, outcome = "Y", treatment = "Test")
  expect_true("subgroup" %in% names(out$data.os))
  expect_true("subgroup" %in% names(out$data.rct))
  expect_equal(nlevels(out$data.os$subgroup), 2L)
})

test_that("subgroupDefs_noSplit crosses multiple variables", {
  d <- make_base_data()
  out <- subgroupDefs_noSplit(c("X1", "X2"), d, d, outcome = "Y", treatment = "Test")
  # 2 x 2 = up to 4 strata
  expect_lte(nlevels(out$data.os$subgroup), 4L)
  expect_gte(nlevels(out$data.os$subgroup), 1L)
})

test_that("subgroupDefs_noSplit returns named subgroupTaus.true", {
  d <- make_base_data()
  out <- subgroupDefs_noSplit("X1", d, d, outcome = "Y", treatment = "Test")
  expect_true(is.numeric(out$subgroupTaus.true))
  expect_false(is.null(names(out$subgroupTaus.true)))
  expect_equal(length(out$subgroupTaus.true), nlevels(out$data.rct$subgroup))
})

test_that("subgroupDefs_noSplit accepts a custom treatment column name", {
  d <- make_base_data(trt_col = "assigned")
  out <- subgroupDefs_noSplit("X1", d, d, outcome = "Y", treatment = "assigned")
  expect_true(all(is.finite(out$subgroupTaus.true)))
  expect_equal(length(out$subgroupTaus.true), nlevels(out$data.rct$subgroup))
})

test_that("subgroupDefs_noSplit true effects are the same regardless of treatment column name", {
  set.seed(2)
  d_default <- make_base_data(trt_col = "Test")

  d_custom <- d_default
  names(d_custom)[names(d_custom) == "Test"] <- "assigned"

  out1 <- subgroupDefs_noSplit("X1", d_default, d_default,
                               outcome = "Y", treatment = "Test")
  out2 <- subgroupDefs_noSplit("X1", d_custom,  d_custom,
                               outcome = "Y", treatment = "assigned")
  expect_equal(out1$subgroupTaus.true, out2$subgroupTaus.true)
})
