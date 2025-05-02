# Test and develop functions to create treatment effects in simulations
context("Treatment Effect Simulations")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  load_all() ## use  this during debugging
}
setDTthreads(1)

test_that("We can create effects as expected", {
  setkey(idat, bF)
  set.seed(12345)
  ##  Should be strictly null here
  idat$y1test_null <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 0, prop_blocks_0 = 1, covariate = "vb1")
  expect_equal(idat$y1test_null, idat$y0)
  ## Here the effects do differ but the shift is a zero sized shift
  idat$y1test_small <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 0, prop_blocks_0 = 0, covariate = "vb1")
  expect_lt(min(idat[, abs(mean(y1test_small - y0)), by = bF]$V1), .1)
  expect_gt(min(idat[, mean(y1test_small - y0), by = bF]$V1), -.5)
  ## Make a large shift in all of the blocks
  idat$y1test_4 <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 4, prop_blocks_0 = 0, covariate = "vb1")
  expect_lt(min(idat[, mean(y1test_4 - y0), by = bF]$V1), 5 * sd(idat$y0))
  expect_gt(min(idat[, mean(y1test_4 - y0), by = bF]$V1), 3 * sd(idat$y0))
  ## Make only some blocks have strictly no effects
  idat$y1test_Zeros <- create_effects(
    idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm_outliers,
    tau_size = 4, prop_blocks_0 = .5, covariate = "vb1"
  )
  idat[, truetau_Zeros := (y1test_Zeros - y0)]
  ## idat[,mean(truetau_Zeros),by=bF]
  expect_equal(sum(idat[, mean(truetau_Zeros), by = bF]$V1 == 0), 5)
  expect_gt(min(idat[truetau_Zeros != 0, mean(truetau_Zeros), by = bF]$V1), 3.5 * sd(idat$y0))
  expect_lt(min(idat[truetau_Zeros != 0, mean(truetau_Zeros), by = bF]$V1), 4.5 * sd(idat$y0))
  ## Now with covariates.
  idat$y1test_cov_groups <- create_effects(
    idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm_covariate,
    tau_size = 4, prop_blocks_0 = 0, covariate = "vb1", by_block = FALSE
  )
  idat[, .(mean(y1test_cov_groups - y0), mean(vb1)), by = bF]
  expect_gt(mean(idat[vb1 > median(vb1), mean(y1test_cov_groups - y0)]), mean(idat[vb1 <= median(vb1), mean(y1test_cov_groups - y0)]))
  idat$y1test_cov_cont <- create_effects(
    idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm_covariate_cont,
    tau_size = 4, prop_blocks_0 = 0, covariate = "vb1", by_block = FALSE
  )
  idat[, .(mean(y1test_cov_cont - y0), mean(vb1)), by = bF]
  expect_gt(mean(idat[vb1 > median(vb1), mean(y1test_cov_cont - y0)]), mean(idat[vb1 <= median(vb1), mean(y1test_cov_cont - y0)]))
  idat$y1test_cov_levels <- create_effects(
    idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm_covariate_levels,
    tau_size = 4, prop_blocks_0 = 0, covariate = "vbnest", by_block = FALSE
  )
  idat[, .(mean(y1test_cov_levels - y0), unique(vbnest)), by = bF]
  effect_by_levels <- idat[, .(mean(y1test_cov_levels - y0)), by = vbnest]
  expect_false(isTRUE(all.equal(effect_by_levels$V1[1], effect_by_levels$V1[2])))
  expect_false(isTRUE(all.equal(effect_by_levels$V1[1], effect_by_levels$V1[3])))
  expect_false(isTRUE(all.equal(effect_by_levels$V1[1], effect_by_levels$V1[4])))
  expect_false(isTRUE(all.equal(effect_by_levels$V1[2], effect_by_levels$V1[3])))
  expect_false(isTRUE(all.equal(effect_by_levels$V1[2], effect_by_levels$V1[4])))
  expect_false(isTRUE(all.equal(effect_by_levels$V1[3], effect_by_levels$V1[4])))
  idat[, Y_levels := Z * y1test_cov_levels + (1 - Z) * y0]
  expect_lt(anova(lm(Y_levels ~ vbnest, data = idat))$`Pr(>F)`[1], .05)
})
