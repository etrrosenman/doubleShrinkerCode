# Shared fixture: a 3-stratum estimator data frame similar to what
# getParamEstimates() returns in the simulation scripts.
ed3 <- data.frame(
  stratum = c("A", "B", "C"),
  rctEst  = c(0.05, 0.10, 0.03),
  rctVar  = c(0.002, 0.003, 0.001),
  obsEst  = c(0.08, 0.07, 0.06),
  obsVar  = c(0.001, 0.002, 0.0015)
)
k3 <- nrow(ed3)

# ── shrinker() ────────────────────────────────────────────────────────────────

test_that("shrinker returns a numeric vector of length k", {
  out <- shrinker(gammaSq = 0.01, etaSq = 0.05, ed = ed3)
  expect_type(out, "double")
  expect_length(out, k3)
})

test_that("shrinker returns all-zero vector when etaSq = 0", {
  out <- shrinker(gammaSq = 0.01, etaSq = 0, ed = ed3)
  expect_equal(out, rep(0, k3))
})

test_that("shrinker output is finite for sensible hyperparameters", {
  out <- shrinker(gammaSq = 0.005, etaSq = 0.02, ed = ed3)
  expect_true(all(is.finite(out)))
})

test_that("shrinker shrinks estimates toward zero (|output| <= |input|)", {
  # When etaSq is small relative to rctVar, strong shrinkage toward 0 expected.
  out <- shrinker(gammaSq = 0.001, etaSq = 0.0001, ed = ed3)
  expect_true(all(abs(out) <= abs(ed3$rctEst) + 1e-10))
})

test_that("shrinker output increases in etaSq (larger signal => less shrinkage)", {
  out_small <- shrinker(gammaSq = 0.005, etaSq = 0.001, ed = ed3)
  out_large <- shrinker(gammaSq = 0.005, etaSq = 1.0,   ed = ed3)
  # All estimates should be closer to the weighted combination for larger etaSq
  expect_true(all(abs(out_large) >= abs(out_small) - 1e-10))
})

# ── eb.mm1() ──────────────────────────────────────────────────────────────────

test_that("eb.mm1 returns a numeric vector of length k", {
  out <- eb.mm1(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("eb.mm1 returnCIComps=TRUE returns a list with required elements", {
  out <- eb.mm1(ed3, k3, lowerBound = TRUE, returnCIComps = TRUE)
  expect_type(out, "list")
  expect_named(out, c("shrinker", "lambda", "a", "c", "etaSqLwr", "gammaSqLwr"))
  expect_length(out$shrinker, k3)
  expect_length(out$lambda,   k3)
  expect_length(out$a,        k3)
  expect_length(out$c,        k3)
})

test_that("eb.mm1 lambda is in [0, 1]", {
  out <- eb.mm1(ed3, k3, returnCIComps = TRUE)
  expect_true(all(out$lambda >= 0 & out$lambda <= 1))
})

test_that("eb.mm1 a is in [0, 1]", {
  out <- eb.mm1(ed3, k3, returnCIComps = TRUE)
  expect_true(all(out$a >= 0 & out$a <= 1))
})

test_that("eb.mm1 lowerBound shrinks more than no bound", {
  out_no  <- eb.mm1(ed3, k3, lowerBound = FALSE)
  out_lb  <- eb.mm1(ed3, k3, lowerBound = TRUE)
  
  abs(out_lb) - abs(out_no)
  
  expect_true(all(abs(out_lb) <= abs(out_no) + 1e-10))
})

# ── eb.mm2() ──────────────────────────────────────────────────────────────────

test_that("eb.mm2 returns a numeric vector of length k", {
  out <- eb.mm2(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("eb.mm2 returnCIComps=TRUE returns a list with required elements", {
  out <- eb.mm2(ed3, k3, lowerBound = TRUE, returnCIComps = TRUE)
  expect_named(out, c("shrinker", "lambda", "a", "c", "etaSqLwr", "gammaSqLwr"))
})

# ── eb.mle() ──────────────────────────────────────────────────────────────────

test_that("eb.mle returns a numeric vector of length k", {
  out <- eb.mle(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("eb.mle returnCIComps=TRUE returns a list with required elements", {
  out <- eb.mle(ed3, k3, lowerBound = TRUE, returnCIComps = TRUE)
  expect_named(out, c("shrinker", "lambda", "a", "c", "etaSqLwr", "gammaSqLwr"))
})

# ── eb.ure() ──────────────────────────────────────────────────────────────────

test_that("eb.ure returns a numeric vector of length k", {
  out <- eb.ure(ed3, k3)
  expect_type(out, "double")
  expect_length(out, k3)
  expect_true(all(is.finite(out)))
})

test_that("eb.ure returnCIComps=TRUE returns a list with required elements", {
  out <- eb.ure(ed3, k3, lowerBound = TRUE, returnCIComps = TRUE)
  expect_named(out, c("shrinker", "lambda", "a", "c", "etaSqLwr", "gammaSqLwr"))
})

# ── consistency across estimators ─────────────────────────────────────────────

test_that("all four eb estimators produce finite estimates on the same input", {
  results <- list(
    mm1 = eb.mm1(ed3, k3, lowerBound = TRUE),
    mm2 = eb.mm2(ed3, k3, lowerBound = TRUE),
    mle = eb.mle(ed3, k3, lowerBound = TRUE),
    ure = eb.ure(ed3, k3, lowerBound = TRUE)
  )
  for (nm in names(results)) {
    expect_true(all(is.finite(results[[nm]])), label = paste(nm, "is finite"))
    expect_length(results[[nm]], k3 )
  }
})
