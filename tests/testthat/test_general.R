# devtools::load_all()

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
  # We simulate a small tree and verify that each child's p-value is at least that of its parent.
  set.seed(1234)
  max_level <- 3
  k <- 3
  t <- 0.5
  tree_dt <- generate_tree_DT(max_level, k, t)
  res <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, N_total = 1000, effect_size = 0.5,
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
  set.seed(1234)
  res_spend <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k, N_total = 1000, effect_size = 0.5,
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
    alpha = 0.05, k = k, N_total = 1000, effect_size = 0.5,
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
    alpha = 0.05, k = k, N_total = 1000, effect_size = 0.5,
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
  res <- simulate_test_DT(tree_dt, alpha = 0.05, k = k, N_total = 1000, effect_size = 0.5)

  dt_sim <- res$treeDT

  # For any parent that is active, if the children branch was gated the
  # children p_val should be NA. We check that for each parent where Simes
  # test fails, none of its children have been assigned a p-value.

  gated_children <- dt_sim[!is.na(parent) & is.na(p_val)]
  expect_true(nrow(gated_children) > 0)

  ## check using the hommel adjustment too
  set.seed(12357)
  res <- simulate_test_DT(tree_dt,
    alpha = 0.05, k = k,
    N_total = 1000, effect_size = 0.5, local_adj_p_fn = local_hommel_all_ps
  )

  dt_sim <- res$treeDT

  # For any parent that is active, if the children branch was gated the
  # children p_val should be NA. We check that for each parent where the
  # hommel adjustment is > alpha, , none of its children have been assigned a
  # p-value.

  gated_children <- dt_sim[!is.na(parent) & is.na(p_val)]
  expect_true(nrow(gated_children) > 0)
})

# We use the next set of parameters in a couple of test blocks
## Checking that the arguments work.
alpha_methods <- c("fixed", "fixed_k_adj", "adaptive_k_adj", "spending", "investing")
final_adj_methods <- c("none", "fdr", "fwer")
local_adj_methods <- c("local_simes", "local_hommel_all_ps", "local_bh_all_ps", "local_unadj_all_ps")
power_decay_vals <- c(TRUE, FALSE)

parms <- as.data.table(expand.grid(
  alpha_method = alpha_methods,
  final_adj_method = final_adj_methods,
  local_adj_method = local_adj_methods,
  power_decay = power_decay_vals,
  stringsAsFactors = FALSE
))
parms[, idx := seq_len(nrow(parms))]
setkey(parms, "idx")

test_that("All arguments work. No errors", {
  set.seed(12345)
  res_t0_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = 5, t = 0, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, effect_size = 0.5,
      power_decay = x$power_decay,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestSim"),
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
      alpha = 0.05, N_total = 1000, effect_size = 0.5,
      power_decay = x$power_decay,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestSim"),
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
      alpha = 0.05, N_total = 1000, effect_size = 0.5,
      power_decay = x$power_decay,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestSim"),
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


## See the tests.notforCRAN.R for tests that are time consuming and check the FWER etc.
## This next is just meant to be a very rough guide

test_that("simulating many p-values does what we expect", {
  skip()
  skip_on_ci()
  skip_on_cran()
  nrow(parms)

  n_sims <- 100

  set.seed(123456)
  res_t0_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = n_sims, t = 0, k = 3, max_level = 3,
      alpha = 0.05, N_total = 1000, effect_size = 0.5,
      power_decay = x$power_decay,
      local_adj_p_fn = getFromNamespace(x[["local_adj_method"]], ns = "TreeTestSim"),
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
  sim_err <- 2 * sqrt(.5 * (1 - .5) / n_sims)
  .05 + sim_err
  summary(res_t0$false_error)
  summary(res_t0$false_error < .05 + sim_err)

  expect_lt(max(res_t0$false_error), .05 + sim_err)

  ## We don't even need to adjust power as we "split" in order to control the
  ## FWER when t=0 or t=1 (which has no errors anyway and is only shown here
  ## for completeness)

  expect_lt(max(res_t0[!(power_decay), false_error]), .05 + sim_err)
  ## And we don't have to use any local control, or any other tactic. Monotonicity, validity (i.e. p~U()), and gating are all we need
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps", false_error]), .05 + sim_err)
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps" & !(power_decay) & final_adj_method == "none", false_error]), .05 + sim_err)
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps" & !(power_decay) & final_adj_method == "none" & alpha_method == "fixed", false_error]), .05 + sim_err)

  ## Check on the bottom up approach
  expect_lt(max(res_t0$bottom_up_false_error), .05 + sim_err)

  ## Power does not have meaning when all hypotheses are true
  expect_equal(unique(res_t0[["power"]]), NaN)
  expect_equal(unique(res_t0[["bottom_up_power"]]), NaN)
})
