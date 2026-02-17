test_that("simulate_test_DT runs with alpha_method = 'adaptive_power_pruned'", {
  set.seed(1)
  dt <- generate_tree_DT(max_level = 3, k = 4, t = 0.3)

  res <- simulate_test_DT(
    dt,
    alpha = 0.05, k = 4, N_total = 1000,
    effect_size = 0.35,
    alpha_method = "adaptive_power_pruned",
    local_adj_p_fn = local_unadj_all_ps,
    monotonicity = FALSE,
    return_details = TRUE
  )

  expect_true(is.list(res))
  expect_true("treeDT" %in% names(res))
  expect_true("sim_res" %in% names(res))
})


test_that("adaptive_power_pruned allocates more alpha after branch pruning", {
  set.seed(1)
  dt <- generate_tree_DT(max_level = 3, k = 4, t = 0.3)

  res_static <- simulate_test_DT(
    dt,
    alpha = 0.05, k = 4, N_total = 1000,
    effect_size = 0.35,
    alpha_method = "adaptive_power",
    local_adj_p_fn = local_unadj_all_ps,
    monotonicity = FALSE,
    return_details = TRUE
  )

  set.seed(1)
  res_pruned <- simulate_test_DT(
    dt,
    alpha = 0.05, k = 4, N_total = 1000,
    effect_size = 0.35,
    alpha_method = "adaptive_power_pruned",
    local_adj_p_fn = local_unadj_all_ps,
    monotonicity = FALSE,
    return_details = TRUE
  )

  # With this seed/configuration, depth-1 has partial survival, so the
  # branch-pruned rule should increase the depth-2 threshold.
  n_level1_tested <- res_static$treeDT[level == 1 & !is.na(p_val), .N]
  expect_true(n_level1_tested > 0)
  expect_true(n_level1_tested < 4)

  static_level2_alpha <- unique(res_static$treeDT[level == 2 & !is.na(alpha_alloc), alpha_alloc])
  pruned_level2_alpha <- unique(res_pruned$treeDT[level == 2 & !is.na(alpha_alloc), alpha_alloc])

  expect_length(static_level2_alpha, 1)
  expect_length(pruned_level2_alpha, 1)
  expect_gt(pruned_level2_alpha, static_level2_alpha)
})


test_that("simulate_many_runs_DT works with adaptive_power_pruned", {
  set.seed(1)
  res <- simulate_many_runs_DT(
    n_sim = 10, t = 0.3, k = 4, max_level = 2,
    alpha = 0.05, N_total = 1000, effect_size = 0.35,
    alpha_method = "adaptive_power_pruned",
    local_adj_p_fn = local_unadj_all_ps,
    monotonicity = FALSE
  )

  expect_true(is.numeric(res))
  expect_true("false_error" %in% names(res))
})
