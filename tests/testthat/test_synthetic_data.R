## Tests for generate_synthetic_experiment and run_synthetic_comparison

test_that("generate_synthetic_experiment returns correct structure", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0.5)

  expect_true(is.list(synth))
  expect_true("idat" %in% names(synth))
  expect_true("bdat" %in% names(synth))
  expect_true(is.data.table(synth$idat))
  expect_true(is.data.table(synth$bdat))
})

test_that("generate_synthetic_experiment produces correct dimensions", {
  set.seed(42)
  k <- 4
  max_level <- 2
  N_per_block <- 30
  synth <- generate_synthetic_experiment(k = k, max_level = max_level,
                                          N_per_block = N_per_block, t = 0.5)

  num_leaves <- k^max_level
  expect_equal(nrow(synth$bdat), num_leaves)
  expect_equal(nrow(synth$idat), num_leaves * N_per_block)
})

test_that("idat has required columns", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0.5)
  idat <- synth$idat
  expected_cols <- c("bF", "id", "y0", "trt", "trtF", "blockF", "lvls_fac", "nonnull")
  for (col in expected_cols) {
    expect_true(col %in% names(idat), info = paste("Missing column:", col))
  }
})

test_that("bdat has required columns", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0.5)
  bdat <- synth$bdat
  expected_cols <- c("bF", "blockF", "nonnull", "lvls_fac", "nb", "hwt")
  for (col in expected_cols) {
    expect_true(col %in% names(bdat), info = paste("Missing column:", col))
  }
})

test_that("treatment is balanced within each block", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0.5)
  balance <- synth$idat[, .(n_trt = sum(trt), n_ctrl = sum(1 - trt)), by = bF]
  expect_true(all(balance$n_trt == balance$n_ctrl))
})

test_that("lvls_fac works with splitSpecifiedFactorMulti", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0.5)
  bdat <- synth$bdat

  # First split should produce multiple groups
  bdat[, split1 := manytestsr::splitSpecifiedFactorMulti(bid = as.numeric(as.character(bF)), x = lvls_fac)]

  n_groups <- length(unique(bdat$split1[bdat$split1 != 0]))
  # Should get at least 2 groups (splitting works)
  expect_true(n_groups >= 2,
    info = paste("Expected >= 2 groups, got", n_groups))

  # Second split within first groups should also produce subgroups
  bdat[, split2 := manytestsr::splitSpecifiedFactorMulti(
    bid = as.numeric(as.character(bF)), x = lvls_fac
  ), by = split1]
  n_subgroups <- bdat[split2 != 0, uniqueN(interaction(split1, split2))]
  expect_true(n_subgroups >= 2,
    info = paste("Expected >= 2 subgroups, got", n_subgroups))
})

test_that("t = 0 produces all null blocks", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0)
  expect_true(all(synth$bdat$nonnull == FALSE))
})

test_that("t = 1 produces all non-null blocks", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 1)
  expect_true(all(synth$bdat$nonnull == TRUE))
})

test_that("bdat is keyed by bF", {
  set.seed(42)
  synth <- generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 20, t = 0.5)
  expect_equal(key(synth$bdat), "bF")
})

test_that("generate_synthetic_experiment validates inputs", {
  expect_error(generate_synthetic_experiment(k = 1, max_level = 2))
  expect_error(generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 3))
  expect_error(generate_synthetic_experiment(k = 3, max_level = 2, N_per_block = 21))
  expect_error(generate_synthetic_experiment(k = 3, max_level = 2, t = -0.1))
})
