# Shared fixture
ed3 <- data.frame(
  stratum = c("A", "B", "C"),
  rctEst  = c(0.05, 0.10, 0.03),
  rctVar  = c(0.002, 0.003, 0.001),
  obsEst  = c(0.08, 0.07, 0.06),
  obsVar  = c(0.001, 0.002, 0.0015)
)
k3 <- nrow(ed3)

# ── precistion.wtd() ──────────────────────────────────────────────────────────

test_that("precistion.wtd returns a numeric vector of length k", {
  out <- precistion.wtd(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("precistion.wtd output is convex combination of rctEst and obsEst", {
  out <- precistion.wtd(ed3, k3)
  lo <- pmin(ed3$rctEst, ed3$obsEst)
  hi <- pmax(ed3$rctEst, ed3$obsEst)
  expect_true(all(out >= lo - 1e-10 & out <= hi + 1e-10))
})

# ── kappa.1() ─────────────────────────────────────────────────────────────────

test_that("kappa.1 returns a numeric vector of length k", {
  out <- kappa.1(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("kappa.1 equals rctEst when rctEst == obsEst (no bias, lambda=1)", {
  ed_nobias <- ed3
  ed_nobias$obsEst <- ed_nobias$rctEst
  out <- kappa.1(ed_nobias, k3)
  expect_equal(out, ed_nobias$rctEst)
})

# ── kappa.2() ─────────────────────────────────────────────────────────────────

test_that("kappa.2 returns a numeric vector of length k", {
  out <- kappa.2(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

# ── delta.1() ─────────────────────────────────────────────────────────────────

test_that("delta.1 returns a numeric vector of length k (k >= 3)", {
  out <- delta.1(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("delta.1 lambda is in [0, 1], so output is a convex combination", {
  out <- delta.1(ed3, k3)
  lo <- pmin(ed3$rctEst, ed3$obsEst)
  hi <- pmax(ed3$rctEst, ed3$obsEst)
  expect_true(all(out >= lo - 1e-10 & out <= hi + 1e-10))
})

# ── delta.2() ─────────────────────────────────────────────────────────────────

test_that("delta.2 returns a numeric vector of length k", {
  out <- delta.2(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

# ── cross-estimator sanity ────────────────────────────────────────────────────

test_that("all competitor estimators are finite on the same input", {
  results <- list(
    pw     = precistion.wtd(ed3, k3),
    kappa1 = kappa.1(ed3, k3),
    kappa2 = kappa.2(ed3, k3),
    delta1 = delta.1(ed3, k3),
    delta2 = delta.2(ed3, k3)
  )
  for (nm in names(results)) {
    expect_true(all(is.finite(results[[nm]])), label = paste(nm, "is finite"))
  }
})

test_that("kappa.1 and kappa.2 agree when rctVar is uniform (scalar = diagonal)", {
  ed_uniform <- data.frame(
    rctEst = c(0.05, 0.10, 0.03),
    rctVar = rep(0.002, 3),
    obsEst = c(0.08, 0.07, 0.06),
    obsVar = c(0.001, 0.002, 0.0015)
  )
  k <- nrow(ed_uniform)
  out1 <- kappa.1(ed_uniform, k)
  out2 <- as.vector(kappa.2(ed_uniform, k))
  expect_equal(out1, out2, tolerance = 1e-10)
})
