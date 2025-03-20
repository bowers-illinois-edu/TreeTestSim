testthat::context("Data.table-based Tree Simulation Functions")

devtools::load_all()

test_that("generate_tree_DT creates a complete tree with correct columns and levels", {
  max_level <- 2
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)

  # Expected number of nodes is sum(k^l) for l=0:max_level.
  expected_rows <- sum(k^(0:max_level))
  expect_equal(nrow(tree_dt), expected_rows)

  # The data.table should have these columns.
  expect_true(all(c("node", "level", "parent", "nonnull") %in% names(tree_dt)))

  # Check that the level assignments are correct:
  expect_equal(tree_dt[level == 0, .N], 1)
  expect_equal(tree_dt[level == 1, .N], k)
  expect_equal(tree_dt[level == 2, .N], k^2)

  # Check that the parent indices are as expected:
  # For node i > 1, parent = floor((i-2)/k) + 1.
  tree_dt[, computed_parent := ifelse(node == 1, NA_integer_, floor((node - 2) / k) + 1)]
  expect_equal(tree_dt$parent, tree_dt$computed_parent)
})

test_that("generate_tree_DT propagates nonnull flags correctly for extreme t", {
  max_level <- 3
  k <- 3

  # t = 0: all leaves false => all nodes should be false.
  dt0 <- generate_tree_DT(max_level, k, t = 0)
  expect_true(all(dt0[level == max_level, nonnull] == FALSE))
  expect_true(all(dt0$nonnull == FALSE))

  # t = 1: all leaves true => all nodes should be true.
  dt1 <- generate_tree_DT(max_level, k, t = 1)
  expect_true(all(dt1[level == max_level, nonnull] == TRUE))
  expect_true(all(dt1$nonnull == TRUE))

  # t=.5, about .5 of the leaves should be TRUE
  dt_half <- generate_tree_DT(max_level, k, t = .5)
  expect_equal(dt_half[level == max_level, sum(nonnull)], floor(.5 * nrow(dt_half[level == max_level])))
})

test_that("simulate_test_DT produces monotonic p-values", {
  # We simulate a small tree and verify that each child’s p-value is at least that of its parent.
  set.seed(1234)
  max_level <- 3
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)
  res <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1,
    local_adj_p_fn = local_simes, global_adj = "hommel", return_details = TRUE,
    alpha_method = "fixed", final_global_adj = "none"
  )

  dt_sim <- res$treeDT

  # For each node with a parent, we require:
  #    child p_val >= parent's p_val.
  # Merge dt_sim with itself to obtain parent's p-values.
  parent_p_vals <- dt_sim[, .(node, parent_p = p_val)]
  children <- dt_sim[!is.na(parent)]
  children <- merge(children, parent_p_vals, by.x = "parent", by.y = "node", all.x = TRUE, sort = FALSE)

  # Test that for each child, p_val is at least as high as parent's p_val.
  expect_true(all(children$p_val >= children$parent_p, na.rm = TRUE))

  ## Check that this holds when we allow alpha to vary
  ## TODO: I say "spending" here but I'm doing something specific that I think it not exactly
  ## the same. Also this means that I'm ignoring some of the spend_frac, etc.. arguments
  set.seed(1234)
  res_spend <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = TRUE, alpha_method = "spending"
  )

  dt_sim <- res_spend$treeDT

  # For each node with a parent, we require:
  #    child p_val >= parent's p_val.
  # Merge dt_sim with itself to obtain parent's p-values.
  parent_p_vals <- dt_sim[, .(node, parent_p = p_val)]
  children <- dt_sim[!is.na(parent)]
  children <- merge(children, parent_p_vals, by.x = "parent", by.y = "node", all.x = TRUE, sort = FALSE)
  # Test that for each child, p_val is at least as high as parent's p_val.
  expect_true(all(children$p_val >= children$parent_p, na.rm = TRUE))

  ### Now try the less conservative approach
  set.seed(1234)
  res_invest <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = TRUE, alpha_method = "investing"
  )

  dt_sim <- res_invest$treeDT

  # For each node with a parent, we require:
  #    child p_val >= parent's p_val.
  # Merge dt_sim with itself to obtain parent's p-values.
  parent_p_vals <- dt_sim[, .(node, parent_p = p_val)]
  children <- dt_sim[!is.na(parent)]
  children <- merge(children, parent_p_vals, by.x = "parent", by.y = "node", all.x = TRUE, sort = FALSE)

  # Test that for each child, p_val is at least as high as parent's p_val.
  expect_true(all(children$p_val >= children$parent_p, na.rm = TRUE))

  ## And the more conservative approach
  set.seed(1234)
  res_fixed_adj <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = TRUE, alpha_method = "fixed_k_adj"
  )

  dt_sim <- res_invest$treeDT

  # For each node with a parent, we require:
  #    child p_val >= parent's p_val.
  # Merge dt_sim with itself to obtain parent's p-values.
  parent_p_vals <- dt_sim[, .(node, parent_p = p_val)]
  children <- dt_sim[!is.na(parent)]
  children <- merge(children, parent_p_vals, by.x = "parent", by.y = "node", all.x = TRUE, sort = FALSE)

  # Test that for each child, p_val is at least as high as parent's p_val.
  expect_true(all(children$p_val >= children$parent_p, na.rm = TRUE))
})

test_that("simulate_test_DT gates branches when the local adjusted p-value exceeds alpha", {
  # Create a tree and force the children p-values to be high so that the local
  # adjusted test fails. One way to do this is to temporarily override runif()
  # so that it returns 1.

  max_level <- 3
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)

  # Override runif locally.
  local_runif <- function(n, min, max) {
    rep(max, n)
  }
  old_runif <- stats::runif
  unlockBinding("runif", as.environment("package:stats"))
  assign("runif", local_runif, envir = as.environment("package:stats"))
  on.exit({
    assign("runif", old_runif, envir = as.environment("package:stats"))
    lockBinding("runif", as.environment("package:stats"))
  })

  ## recall that the mean of a beta distribution is a/(a+1)
  set.seed(12357)
  res <- simulate_test_DT(tree_dt, alpha = 0.05, k = k, effN = 1000, N_total = 1000, beta_base = 0.1)

  dt_sim <- res$treeDT

  # For any parent that is active, if the children branch was gated the
  # children p_val should be NA. We check that for each parent where Simes
  # test fails, none of its children have been assigned a p-value.

  gated_children <- dt_sim[!is.na(parent) & is.na(p_val)]
  expect_true(nrow(gated_children) > 0)

  ## check using the hommel adjustment too
  set.seed(12357)
  res <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, effN = 1000,
    N_total = 1000, beta_base = 0.1, local_adj_p_fn = local_hommel_all_ps
  )

  dt_sim <- res$treeDT

  # For any parent that is active, if the children branch was gated the
  # children p_val should be NA. We check that for each parent where the
  # hommel adjustment is > alpha, , none of its children have been assigned a
  # p-value.

  gated_children <- dt_sim[!is.na(parent) & is.na(p_val)]
  expect_true(nrow(gated_children) > 0)
})

test_that("All arguments work. No errors", {
  ## Checking that the arguments work.
  alpha_methods <- c("fixed", "fixed_k_adj", "adaptive_k_adj", "spending", "investing")
  final_adj_methods <- c("none", "fdr", "fwer")
  local_adj_methods <- c("local_simes", "local_hommel_all_ps", "local_unadj_all_ps")
  adj_effN <- c(TRUE, FALSE)

  parms <- as.data.table(expand.grid(
    alpha_method = alpha_methods,
    final_adj_method = final_adj_methods,
    local_adj_method = local_adj_methods,
    adj_effN = adj_effN,
    stringsAsFactors = FALSE
  ))
  parms[, idx := seq_len(nrow(parms))]
  setkey(parms, "idx")

  set.seed(12345)
  res_t0_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 5, t = 0, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method, multicore = TRUE
    )
    return(tmp)
  })

  set.seed(12345)
  res_t1_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 5, t = 1, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method
    )
    return(tmp)
  })

  set.seed(12345)
  res_t_half_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 5, t = .5, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method
    )
    return(tmp)
  })

  res_t0 <- do.call("rbind", res_t0_lst)
  res_t1 <- do.call("rbind", res_t1_lst)
  res_t_half <- do.call("rbind", res_t_half_lst)

  ## basically this is a test that it ran without errors
  expect_equal(nrow(res_t0), nrow(parms))
  expect_equal(nrow(res_t1), nrow(parms))
  expect_equal(nrow(res_t_half), nrow(parms))
})


## TODO: Don't do 10,000 sims on CRAN

test_that("simulating many p-values does what we expect", {
  alpha_methods <- c("fixed", "fixed_k_adj", "adaptive_k_adj", "spending", "investing")
  final_adj_methods <- c("none", "fdr", "fwer")
  local_adj_methods <- c("local_simes", "local_hommel_all_ps", "local_unadj_all_ps")
  adj_effN <- c(TRUE, FALSE)

  parms <- as.data.table(expand.grid(
    alpha_method = alpha_methods,
    final_adj_method = final_adj_methods,
    local_adj_method = local_adj_methods,
    adj_effN = adj_effN,
    stringsAsFactors = FALSE
  ))
  parms[, idx := seq_len(nrow(parms))]
  setkey(parms, "idx")
  nrow(parms)

  n_sims <- 10000

  set.seed(123456)
  res_t0_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 10000, t = 0, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method, multicore = TRUE
    )
    x[, names(tmp) := as.list(tmp)]
    return(x)
  })

  res_t0 <- rbindlist(res_t0_lst)

  ## Should control res1 within simulation error and we are only doing 1000 sims here
  sim_err <- 2 * sqrt(.5 * (1 - .5) / 10000)
  .05 + sim_err
  summary(res_t0$false_error)
  summary(res_t0$false_error < .05 + sim_err)

  expect_lt(max(res_t0$false_error), .05 + sim_err)

  ## We don't even need to adjust power as we "split" in order to control the
  ## FWER when t=0 or t=1 (which has no errors anyway and is only shown here
  ## for completeness)

  expect_lt(max(res_t0[!(adj_effN), false_error]), .05 + sim_err)
  ## And we don't have to use any local control, or any other tactic. Monotonicity, validity (i.e. p~U()), and gating are all we need
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps", false_error]), .05 + sim_err)
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps" & !(adj_effN) & final_adj_method == "none", false_error]), .05 + sim_err)
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps" & !(adj_effN) & final_adj_method == "none" & alpha_method == "fixed", false_error]), .05 + sim_err)

  ## Check on the bottom up approach
  expect_lt(max(res_t0$bottom_up_false_error), .05 + sim_err)

  ## Power does not have meaning when all hypotheses are true
  expect_equal(unique(res_t0[["power"]]), NaN)
  expect_equal(unique(res_t0[["bottom_up_power"]]), NaN)

  ## This next is all about power. Should have no false positive rate since all
  ## hypotheses are false


  set.seed(123456)
  res_t1_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 10000, t = 1, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method, multicore = TRUE
    )
    x[, names(tmp) := as.list(tmp)]
    return(x)
  })

  res_t1 <- rbindlist(res_t1_lst)

  expect_equal(unique(res_t1$false_error), 0)
  expect_equal(unique(res_t1$bottom_up_false_error), 0)

  ## We don't even need to adjust power as we "split" in order to control the
  ## FWER when t=0 or t=1 (which has no errors anyway and is only shown here
  ## for completeness)

  ## Check on the bottom up approach
  expect_lt(max(res_t1$bottom_up_false_error), .05 + sim_err)

  ## Power does not have meaning when all hypotheses are true

  summary(res_t1)
  ## The top down method rejects more nodes and more leaves
  expect_true(all(res_t1$bottom_up_power - res_t1$power < 0))
  expect_true(all(res_t1$bottom_up_power - res_t1$leaf_power < 0))

  ### So, weak control works even when we don't split the data at each node.
  ### This is not realistic. But nice to know.

  ## Now try it with t=.1 which is .1*.27 \approx 2 leaves

  set.seed(123456)
  res_t_some_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 1000, t = .1, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method, multicore = TRUE
    )
    x[, names(tmp) := as.list(tmp)]
    return(x)
  })
  sim_err <- 2 * sqrt(.5 * (1 - .5) / 1000)
  res_t_some <- rbindlist(res_t_some_lst)

  summary(res_t_some)

  ## Which approaches control the FWER?
  parm_nms <- names(res_t_some)[1:4]

  res_t_some_ok_fwer <- res_t_some[false_error <= .05 + sim_err, ]

  ## Multiple possibilities:
  lapply(res_t_some_ok_fwer[, .SD, .SDcols = parm_nms], table)

  ## Which have the highest power?
  res_t_some_ok_fwer[order(power, decreasing = TRUE), .SD, .SDcols = c("power", "leaf_power", parm_nms)]

  ## What among those without a final adjustment method?

  res_t_some_ok_fwer[final_adj_method == "none", ][order(power, decreasing = TRUE), .SD, .SDcols = c("power", "leaf_power", parm_nms)]

  ## And without any local adjustment but with an adaptive alpha adjustment

  res_t_some_ok_fwer[final_adj_method == "none" & local_adj_method == "local_unadj_all_ps", ][order(power, decreasing = TRUE), .SD, .SDcols = c("power", "leaf_power", parm_nms)]


  ## Compare to a situation with more k and more opportunities for ungating (t=.5)

  set.seed(123456)
  (10^(3 + 1) - 1) / (10 - 1)
  sum(10^(0:3))
  sum(10^(0:4))

  res_t_half_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 1000, t = .5, k = 10, max_level = 4,
      alpha = 0.05, N_total = 1000, beta_base = 0.1,
      adj_effN = x$adj_effN,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestsSim"),
      global_adj = "hommel",
      return_details = FALSE,
      final_global_adj = x$final_adj_method,
      alpha_method = x$alpha_method, multicore = TRUE
    )
    x[, names(tmp) := as.list(tmp)]
    return(x)
  })
  sim_err <- 2 * sqrt(.5 * (1 - .5) / 1000)
  res_t_half <- rbindlist(res_t_half_lst)
  save(res_t_half, res_t_some, file = "res_t_tests.rda")
  summary(res_t_half)

  ## Which approaches control the FWER?
  parm_nms <- names(res_t_half)[1:4]

  res_t_half_ok_fwer <- res_t_half[false_error <= .05 + sim_err, ]

  ## Multiple possibilities:
  lapply(res_t_half_ok_fwer[, .SD, .SDcols = parm_nms], table)

  ## Which have the highest power?
  res_t_half_ok_fwer[order(power, decreasing = TRUE), .SD, .SDcols = c("power", "leaf_power", "bottom_up_power", parm_nms)]

  ## What among those without a final adjustment method?

  res_t_half_ok_fwer[final_adj_method == "none", ][order(power, decreasing = TRUE), .SD, .SDcols = c("power", "leaf_power", parm_nms)]

  ## And without any local adjustment but with an adaptive alpha adjustment

  res_t_half_ok_fwer[final_adj_method == "none" & local_adj_method ==
    "local_unadj_all_ps", ][order(power, decreasing = TRUE), .SD,
    .SDcols =
      c("power", "leaf_power", "num_leaves_tested", "num_leaves", parm_nms)
  ]

  res_t_some_ok_fwer[final_adj_method == "none" & local_adj_method ==
    "local_unadj_all_ps", ][order(power, decreasing = TRUE), .SD,
    .SDcols =
      c("power", "leaf_power", "num_leaves_tested", "num_leaves", parm_nms)
  ]

  ## This suggests:
  ## (1) you can avoid local_adj_methods across all children if you use the adaptive_k_adj and sample splitting.
  ## (2) with no final global adjustment, you can have power with TODO

  ## TODO: choose the best ideas for the broader sim to look at changes with k and l and t
})


##
