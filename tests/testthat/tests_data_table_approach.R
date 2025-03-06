context("Data.table-based Tree Simulation Functions")

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


test_that("simulating many p-values does what we expect", {
  set.seed(123456)

  ## We don't even need to adjust power as we "split" in order to control the
  ## FWER when t=0 or t=1 (which has no errors anyway and is only shown here
  ## for completeness)

  res1 <- simulate_many_runs_DT(
    n_sim = 1000, t = 0, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = FALSE, final_global_adj = "none", alpha_method = "fixed"
  )

  ## On average does only one single test and doesn't reject
  res1
  expect_true(is.numeric(res1))
  ## Should control res1 within simulation error
  sim_err <- 2 * sqrt(.05 * (1 - .05) / 1000)

  expect_lt(res1["false_error"], .05 + sim_err)
  expect_lt(res1["bottom_up_false_error"], .05 + sim_err)

  ## Power does not have meaning when all hypotheses are true
  expect_equal(res1[["power"]], NaN)
  expect_equal(res1[["bottom_up_power"]], NaN)

  ## Checking that the arguments work.
  alpha_methods <- c("fixed", "fixed_k_adj", "adjust_k_adj", "spending", "investing")
  final_adj_methods <- c("none", "fdr", "fwer")
  local_adj_methods <- c("local_hommel_all_ps", "local_unadj_all_ps")

  parms_t0 <- expand.grid(alpha_method = alpha_methods, final_adj_method = final_adj_methods, local_adj_method = local_adj_methods)

  ## TODO: Start here

  ## now for t=1
  set.seed(12345)
  res2 <- simulate_many_runs_DT(
    n_sim = 1000, t = 1, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = FALSE
  )

  ## Total nodes
  sum(3^seq(0, 3))
  (3^(3 + 1) - 1) / (3 - 1)
  ## Total leaves
  3^3

  ## On average we test about 9 nodes --- all of which are nonnull no surprise
  ## here. On average we make about 6 true discoveries. On average the
  ## proportion of true discoveries of non-null nodes tested is about .55. We
  ## do get to the leaves in this case --- tending to test about 2 of the 27
  ## leaves and rejecting correcting about 2.5 of them.

  res2
  ## No false positives possible here
  ## We should have more power than the bottom up approach
  ## Although we do fewer tests overall
  expect_lt(res2["false_error"], .05 + sim_err)
  expect_lt(res2["bottom_up_false_error"], .05 + sim_err)
  expect_lt(res2["bottom_up_power"], res2["power"])
  expect_lt(res2["bottom_up_power"], res2["leaf_power"])

  ### So, weak control works even when we don't split the data at each node.
  ### This is not realistic. But nice to know.

  ## Now "split the data" (i.e. reduce propostion of p<.05 from the nonnull
  ## nodes by increasing the first parameter of rbeta)

  set.seed(12345)
  res3 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = FALSE
  )
  res3
  ## Here we test fewer nodes
  expect_lt(res3["false_error"], .05 + sim_err)
  expect_lt(res3["bottom_up_false_error"], .05 + sim_err)
  expect_lt(res3["bottom_up_power"], res3["leaf_power"])

  set.seed(12345)
  res3a <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_simes, global_adj = "hommel", return_details = FALSE
  )
  res3a

  set.seed(12345)
  res3b <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE,
    local_adj_p_fn = local_simes, global_adj = "hommel", return_details = FALSE
  )
  res3b
  ## Notice that strong FWER control ot hold with adj_effN=FALSE. So, we need
  ## monotonicity and local gate but also need to reduce power when effects
  ## are mixed.
  res4 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE,
    local_adj_p_fn = local_hommel_all_ps, global_adj = "hommel", return_details = FALSE
  )
  res4
  expect_gt(res4["false_error"], .05 + sim_err)
  expect_lt(res4["bottom_up_false_error"], .05 + sim_err)

  ## Notice that the reduction in power is doing a lot of work here even
  ## without the local adjustment. For example, even with a large k and large l
  ## the algorithm to reduce N is nonlinear -- starts by equal splitting by k
  ## but never goes to 0 (it goes no further than 1/100 of the original).
  ## this next:

  # number of nodes is large here
  (5^(5 + 1) - 1) / (5 - 1)

  set.seed(123456)
  res5 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_unadj_all_ps,
    global_adj = "hommel", return_details = FALSE
  )
  res5

  ## So,the FWER is like .10 in a tree with 3000+ nodes even without local
  ## adjustment, just monotonicity and sample splitting

  expect_gt(res5["false_error"], .05 + sim_err)

  ## local hommel doesn't really help that much
  set.seed(123456)
  res6 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 5, max_level = 5,
    alpha = 0.05, N_total = 1000, beta_base = 0.1, adj_effN = FALSE,
    local_adj_p_fn = local_hommel_all_ps,
    global_adj = "hommel", return_details = FALSE
  )
  res6
  expect_gt(res6["false_error"], .05 + sim_err)



  ## So you need the splitting and monotonicity (which does the majority of the error control work),
  ## and then the local hommel p-value adjustment

  set.seed(123456)
  res7 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 5, max_level = 5,
    alpha = 0.05, N_total = 10000000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_hommel_all_ps,
    global_adj = "hommel", return_details = FALSE
  )
  res7
  ## Notice here strong control although the tests do not get down to the leaves
  expect_lt(res7["false_error"], .05 + sim_err)

  ## Now, look at the kind of alpha spending/investing idea -- moderate increase in power
  ## reduce tree size for speed
  set.seed(123456)
  res8 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 10000000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_hommel_all_ps,
    global_adj = "hommel", return_details = FALSE, alpha_method = "spending"
  )
  res8
  ## Notice here strong control although the tests do not get down to the leaves
  expect_lt(res8["false_error"], .05 + sim_err)
  # expect_gt(res8["power"], res7["power"])

  ### Now look at the alpha investing approach --- no real difference
  set.seed(123456)
  res9 <- simulate_many_runs_DT(
    n_sim = 1000, t = .5, k = 3, max_level = 3,
    alpha = 0.05, N_total = 10000000, beta_base = 0.1, adj_effN = TRUE,
    local_adj_p_fn = local_hommel_all_ps,
    global_adj = "hommel", return_details = FALSE, alpha_method = "investing"
  )
  res9
  ## Notice here strong control although the tests do not get down to the leaves
  expect_lt(res9["false_error"], .05 + sim_err)
})


##
