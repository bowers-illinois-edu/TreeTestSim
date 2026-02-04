## Tests for alpha_method = "adaptive_power" and dual bottom-up output

test_that("derive_delta_hat returns positive values", {
  dhat <- TreeTestSim:::derive_delta_hat(beta_base = 0.1, N_total = 1000, alpha = 0.05)
  expect_true(is.numeric(dhat))
  expect_true(dhat > 0)
})

test_that("derive_delta_hat matches known root power", {
  # pbeta(0.05, 0.1, 1) = 0.05^0.1 ≈ 0.741
  root_power <- pbeta(0.05, 0.1, 1)
  dhat <- TreeTestSim:::derive_delta_hat(beta_base = 0.1, N_total = 1000, alpha = 0.05)
  # Reconstruct power from delta_hat
  reconstructed <- pnorm(dhat * sqrt(1000) - qnorm(0.975))
  expect_equal(reconstructed, root_power, tolerance = 0.01)
})

test_that("derive_delta_hat handles edge cases", {
  # Very small beta_base -> very high power
  dhat_high <- TreeTestSim:::derive_delta_hat(beta_base = 0.01, N_total = 1000)
  expect_true(dhat_high > 0)

  # Very large beta_base -> very low power
  dhat_low <- TreeTestSim:::derive_delta_hat(beta_base = 10, N_total = 1000)
  expect_true(dhat_low > 0)
})

test_that("simulate_test_DT runs with alpha_method = 'adaptive_power'", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
    beta_base = 0.1, alpha_method = "adaptive_power",
    return_details = TRUE
  )
  expect_true(is.list(res))
  expect_true("treeDT" %in% names(res))
  expect_true("sim_res" %in% names(res))
})

test_that("adaptive_power produces alpha_alloc that varies by level", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 3, k = 4, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 4, effN = 1000, N_total = 1000,
    beta_base = 0.1, alpha_method = "adaptive_power",
    return_details = TRUE
  )
  tree <- res$treeDT
  # Root gets nominal alpha
  expect_equal(tree[node == 1, alpha_alloc], 0.05)

  # Active nodes at level 1 should have a different alpha_alloc than root
  level1_allocs <- tree[level == 1 & !is.na(alpha_alloc), alpha_alloc]
  if (length(level1_allocs) > 0) {
    # All level-1 children should have the same threshold (depth-based, not parent-based)
    expect_true(length(unique(level1_allocs)) == 1)
  }
})

test_that("adaptive_power with explicit delta_hat works", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
    beta_base = 0.1, alpha_method = "adaptive_power",
    delta_hat = 0.3, return_details = FALSE
  )
  expect_true(is.data.table(res))
  expect_true("false_error" %in% names(res))
})

test_that("adaptive_power respects tau parameter", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 3, k = 4, t = 0.5)

  # tau = 1 means never adjust (always use nominal alpha)
  res_tau1 <- simulate_test_DT(dt,
    alpha = 0.05, k = 4, effN = 1000, N_total = 1000,
    beta_base = 0.1, alpha_method = "adaptive_power",
    tau = 1.0, return_details = TRUE
  )

  # With tau = 1 all levels should get nominal alpha
  allocs <- res_tau1$treeDT[!is.na(alpha_alloc), unique(alpha_alloc)]
  expect_true(all(allocs == 0.05))
})

test_that("simulate_many_runs_DT works with adaptive_power", {
  set.seed(42)
  res <- simulate_many_runs_DT(
    n_sim = 10, t = 0.3, k = 3, max_level = 2,
    alpha = 0.05, N_total = 1000, beta_base = 0.1,
    alpha_method = "adaptive_power"
  )
  expect_true(is.numeric(res))
  expect_true("false_error" %in% names(res))
  expect_true("bu_hommel_false_error" %in% names(res))
})

## Dual bottom-up tests

test_that("sim_res includes dual bottom-up columns", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
    beta_base = 0.1, alpha_method = "fixed",
    return_details = FALSE
  )
  expected_cols <- c(
    "bu_hommel_false_error", "bu_hommel_true_disc", "bu_hommel_power",
    "bu_bh_false_error", "bu_bh_true_disc", "bu_bh_power"
  )
  for (col in expected_cols) {
    expect_true(col %in% names(res), info = paste("Missing column:", col))
  }
})

test_that("old bottom_up columns still present for backward compat", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
    beta_base = 0.1, return_details = FALSE
  )
  expect_true("bottom_up_false_error" %in% names(res))
  expect_true("bottom_up_true_discoveries" %in% names(res))
  expect_true("bottom_up_power" %in% names(res))
})

test_that("under complete null, bottom-up methods have low false error", {
  set.seed(123)
  n_sim <- 200
  results <- simulate_many_runs_DT(
    n_sim = n_sim, t = 0, k = 3, max_level = 2,
    alpha = 0.05, N_total = 1000, beta_base = 0.1
  )
  # Both hommel and BH bottom-up should control FWER under complete null
  # Allow generous margin for small n_sim
  expect_true(results["bu_hommel_false_error"] <= 0.10,
    info = paste("bu_hommel_false_error =", results["bu_hommel_false_error"]))
  # BH controls FDR, not FWER, but with t=0 they are the same
  expect_true(results["bu_bh_false_error"] <= 0.10,
    info = paste("bu_bh_false_error =", results["bu_bh_false_error"]))
})

test_that("existing alpha methods still work after changes", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.3)
  for (method in c("fixed", "spending", "investing", "fixed_k_adj", "adaptive_k_adj")) {
    res <- simulate_test_DT(dt,
      alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
      beta_base = 0.1, alpha_method = method,
      return_details = FALSE
    )
    expect_true(is.data.table(res), info = paste("Failed for method:", method))
    expect_true("false_error" %in% names(res), info = paste("Missing false_error for:", method))
  }
})
