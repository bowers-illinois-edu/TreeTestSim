## Tests for the effect_size (Cohen's d) parameterization
## These encode statistical principles that the implementation must satisfy.

test_that("effect_size_to_beta matches hand-computed power", {
  # For known d, N, alpha, the beta parameter should produce the correct
  # rejection rate: pbeta(alpha, beta, 1) == expected_power
  d <- 0.5
  N <- 200
  alpha <- 0.05

  z_crit <- qnorm(1 - alpha / 2)
  expected_power <- pnorm(d * sqrt(N / 4) - z_crit)

  beta_param <- TreeTestSim:::effect_size_to_beta(d, N, alpha)
  actual_power <- pbeta(alpha, beta_param, 1)

  expect_equal(actual_power, expected_power, tolerance = 1e-6)
})

test_that("power decays with tree depth", {
  # For k=3, effect_size_to_beta should return increasing beta values

  # (weaker power) at deeper levels because sample size shrinks
  d <- 0.5
  N_total <- 1000
  k <- 3
  alpha <- 0.05
  max_level <- 3

  betas <- sapply(0:max_level, function(l) {
    TreeTestSim:::effect_size_to_beta(d, N_total / k^l, alpha)
  })

  # Beta values should be monotonically increasing (weaker signal)
  for (i in seq_along(betas)[-1]) {
    expect_gt(betas[i], betas[i - 1],
      label = paste("beta at level", i - 1, "vs level", i - 2))
  }
})

test_that("root power matches normal formula", {
  d <- 0.5
  N_total <- 1000
  alpha <- 0.05

  z_crit <- qnorm(1 - alpha / 2)
  expected_root_power <- pnorm(d * sqrt(N_total / 4) - z_crit)

  root_beta <- TreeTestSim:::effect_size_to_beta(d, N_total, alpha)
  actual_root_power <- pbeta(alpha, root_beta, 1)

  expect_equal(actual_root_power, expected_root_power, tolerance = 1e-6)
})

test_that("FWER control under complete null (t=0)", {
  set.seed(8675309)
  n_sim <- 200
  alpha <- 0.05

  results <- simulate_many_runs_DT(
    n_sim = n_sim, t = 0, k = 3, max_level = 3,
    alpha = alpha, N_total = 1000, effect_size = 0.5
  )

  sim_err <- 2 * sqrt(0.5 * (1 - 0.5) / n_sim)
  expect_lt(results["false_error"], alpha + sim_err)
})

test_that("power under complete alternative (t=1)", {
  set.seed(8675309)

  results <- simulate_many_runs_DT(
    n_sim = 50, t = 1, k = 3, max_level = 2,
    alpha = 0.05, N_total = 1000, effect_size = 0.5
  )

  # No false errors when all hypotheses are truly non-null
  expect_equal(unname(results["false_error"]), 0)
  # Should have some power to detect true effects
  expect_gt(results["power"], 0)
})

test_that("monotonicity preserved: child p_val >= parent p_val", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, N_total = 1000, effect_size = 0.5,
    return_details = TRUE
  )

  tree <- res$treeDT
  parent_ps <- tree[, .(node, parent_p = p_val)]
  children <- merge(
    tree[!is.na(parent)],
    parent_ps,
    by.x = "parent", by.y = "node",
    all.x = TRUE, sort = FALSE
  )

  # Every child with a non-NA p-value should have p >= parent's p
  expect_true(all(children$p_val >= children$parent_p, na.rm = TRUE))
})

test_that("power_decay = FALSE uses root-level beta at all depths", {
  # With power_decay = FALSE, the beta parameter should be the same at every
  # level — meaning deeper nodes are drawn from the same-strength distribution.
  # We test indirectly: with a large effect and deep tree, power_decay=FALSE
  # should produce more discoveries than power_decay=TRUE.
  set.seed(42)
  dt <- generate_tree_DT(max_level = 3, k = 3, t = 1)

  res_decay <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, N_total = 1000, effect_size = 0.3,
    power_decay = TRUE, return_details = TRUE
  )
  res_nodecay <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, N_total = 1000, effect_size = 0.3,
    power_decay = FALSE, return_details = TRUE
  )

  # Without decay, non-null p-values at deep levels are drawn from a
  # stronger distribution, so we expect weakly more discoveries
  disc_decay <- res_decay$sim_res$true_discoveries
  disc_nodecay <- res_nodecay$sim_res$true_discoveries
  expect_true(disc_nodecay >= disc_decay)
})

test_that("adaptive alpha schedule receives effect_size / 2 as delta_hat", {
  set.seed(42)
  k <- 3
  max_level <- 3
  N_total <- 1000
  alpha <- 0.05
  effect_size <- 0.5

  dt <- generate_tree_DT(max_level = max_level, k = k, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = alpha, k = k, N_total = N_total,
    effect_size = effect_size,
    alpha_method = "adaptive_power",
    return_details = TRUE
  )

  # Directly compute what the schedule should be
  expected_schedule <- manytestsr::compute_adaptive_alphas(
    k = k, delta_hat = effect_size / 2, N_total = N_total,
    tau = 0.1, max_depth = max_level + 1L, thealpha = alpha
  )

  # Level-1 children should have been assigned the schedule's depth-2 value
  level1_allocs <- res$treeDT[level == 1 & !is.na(alpha_alloc), unique(alpha_alloc)]
  if (length(level1_allocs) > 0) {
    expect_equal(unname(level1_allocs[1]), unname(expected_schedule[2]), tolerance = 1e-10)
  }
})

test_that("all parameter combos run without error", {
  set.seed(42)
  alpha_methods <- c(
    "fixed", "fixed_k_adj", "adaptive_k_adj",
    "spending", "investing", "adaptive_power"
  )
  final_adj_methods <- c("none", "fdr", "fwer")
  local_adj_methods <- c(
    "local_simes", "local_hommel_all_ps",
    "local_bh_all_ps", "local_unadj_all_ps"
  )
  power_decay_vals <- c(TRUE, FALSE)

  parms <- expand.grid(
    alpha_method = alpha_methods,
    final_adj_method = final_adj_methods,
    local_adj_method = local_adj_methods,
    power_decay = power_decay_vals,
    stringsAsFactors = FALSE
  )

  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.5)

  for (i in seq_len(nrow(parms))) {
    p <- parms[i, ]
    local_fn <- getFromNamespace(p$local_adj_method, ns = "TreeTestSim")
    expect_no_error(
      simulate_test_DT(dt,
        alpha = 0.05, k = 3, N_total = 1000, effect_size = 0.5,
        power_decay = p$power_decay, local_adj_p_fn = local_fn,
        alpha_method = p$alpha_method,
        final_global_adj = p$final_adj_method,
        return_details = FALSE
      ),
      message = paste("Failed for:", p$alpha_method, p$final_adj_method,
        p$local_adj_method, p$power_decay)
    )
  }
})

test_that("higher effect_size produces more discoveries", {
  set.seed(42)
  n_sim <- 100

  res_small <- simulate_many_runs_DT(
    n_sim = n_sim, t = 1, k = 3, max_level = 2,
    alpha = 0.05, N_total = 1000, effect_size = 0.2
  )
  res_large <- simulate_many_runs_DT(
    n_sim = n_sim, t = 1, k = 3, max_level = 2,
    alpha = 0.05, N_total = 1000, effect_size = 0.8
  )

  # Larger effect size should produce weakly more discoveries on average
  expect_gte(res_large["true_discoveries"], res_small["true_discoveries"])
})

test_that("simulate_test_DT includes root_power in sim_res", {
  set.seed(42)
  dt <- generate_tree_DT(max_level = 2, k = 3, t = 0.5)
  res <- simulate_test_DT(dt,
    alpha = 0.05, k = 3, N_total = 1000, effect_size = 0.5,
    return_details = TRUE
  )

  expect_true("root_power" %in% names(res$sim_res))

  # root_power should match the normal formula
  z_crit <- qnorm(1 - 0.05 / 2)
  expected <- pnorm(0.5 * sqrt(1000 / 4) - z_crit)
  expect_equal(res$sim_res$root_power, expected, tolerance = 1e-6)
})
