## This next takes too long for ordinary testing especially on CRAN

testthat::context("not for CRAN")

devtools::load_all()
library(dtplyr)
library(dplyr)

test_that("simulating many p-values does what we expect", {
  ## We are using k=4 and l=3 here which is 27 nodes
  # n_nodes_1 <- sum(nodes_per_level^seq(0, num_levels))
  # n_nodes_2 <- ((nodes_per_level^num_levels) - 1) / (nodes_per_level - 1)

  thek <- 4
  thel <- 3

  n_nodes_1 <- sum(thek^seq(0, thel))
  n_nodes_2 <- ((thek^(thel + 1)) - 1) / (thek - 1)
  n_nodes_1
  n_nodes_2
  ## Tests at each level
  4^c(0, 1, 2, 3)
  thek^seq(0, thel)
  num_leaves <- thek^thel
  num_leaves

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

  n_sims <- 1000
  ## Should control error rates etc. within simulation error and we are only doing 1000 sims here
  ## Using n_sims=1000 for speed
  sim_err <- 2 * sqrt(.5 * (1 - .5) / n_sims)
  .05 + sim_err

  set.seed(123456)
  res_t0_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = n_sims, t = 0, k = thek, max_level = thel,
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

  ## Notice some on the high side.
  summary(res_t0$false_error)
  table(res_t0$false_error < .05 + sim_err)
  res_t0[false_error > .07, ]

  expect_lt(max(res_t0$false_error), .05 + sim_err)

  ## Notice that we don't even need to adjust power as we "split" in order to
  ## control the FWER when t=0 or t=1 (which has no errors anyway and is only
  ## shown here for completeness)

  expect_lt(max(res_t0[!(adj_effN), false_error]), .05 + sim_err)
  ## And we don't have to use any local control, or any other tactic.
  ## Monotonicity, validity (i.e. p~U()), and gating are all we need
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps", false_error]), .05 + sim_err)
  expect_lt(max(res_t0[
    local_adj_method == "local_unadj_all_ps" & !(adj_effN) & final_adj_method == "none",
    false_error
  ]), .05 + sim_err)
  expect_lt(max(res_t0[local_adj_method == "local_unadj_all_ps" & !(adj_effN) &
    final_adj_method == "none" &
    alpha_method == "fixed", false_error]), .05 + sim_err)

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
      n_sim = n_sims, t = 1, k = thek, max_level = thel,
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

  summary(res_t1)

  ## The top down method rejects more nodes and more leaves Remember that here
  ## all leaves are false and should be rejected We test fewer leaves than the
  ## bottom up method: here, the maximum number of leaves tested is around 12
  ## (but there are 64 leaves.) That said, when we do test a leaf we tend to
  ## reject the hypothesis.

  expect_true(all(res_t1$bottom_up_power - res_t1$power < 0))
  expect_true(all(res_t1$bottom_up_power - res_t1$leaf_power < 0))

  res_t1[num_leaves_tested > 10, ]

  ## Now try it with t=.1 which is .1*.27 \approx 2 leaves

  set.seed(123456)
  res_t_some_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = n_sims, t = .1, k = thek, max_level = thel,
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
  res_t_some <- rbindlist(res_t_some_lst)

  summary(res_t_some)

  ## Which approaches control the FWER?
  parm_nms <- names(res_t_some)[1:4]

  res_t_some_ok_fwer <- droplevels(res_t_some[false_error <= .05 + sim_err, ])

  ## Multiple possibilities --- all of which control the FWER within simulation error
  lapply(res_t_some_ok_fwer[, .SD, .SDcols = parm_nms], table)

  ## Which have the highest power?
  res_t_some_ok_fwer %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)

  ## local_hommel_all_ps wins in general. Let's look for the set of approaches
  ## that we'd like to send to a bigger simulation

  ## res_t_some_ok_fwer[order(power, decreasing = TRUE), .SD, .SDcols = c("power", "leaf_power", "false_error", parm_nms)]

  ## What about with alpha fixed and no final adjustment? (the original idea)
  ## Below we see that we need data splitting and local hommel
  res_t_some_ok_fwer %>%
    filter(final_adj_method == "none" & alpha_method == "fixed") %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)
  ##       power leaf_power false_error alpha_method final_adj_method    local_adj_method adj_effN
  ##       <num>      <num>       <num>       <char>           <char>              <char>   <lgcl>
  ## 1: 0.7450000  1.0000000       0.022        fixed             none local_hommel_all_ps     TRUE
  ## 2: 0.6121262  0.7121212       0.065        fixed             none         local_simes     TRUE

  ## What among those without a final adjustment method but with a varying alpha and a version of sample splitting
  ## local_simes is clearly worst so exclude here
  res_t_some_ok_fwer %>%
    filter(final_adj_method == "none" & alpha_method != "fixed" & adj_effN & local_adj_method != "local_simes") %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)
  ##   power leaf_power false_error   alpha_method final_adj_method    local_adj_method adj_effN
  ##   <num>      <num>       <num>         <char>           <char>              <char>   <lgcl>
  ## 1: 0.771  0.7285714       0.056      investing             none local_hommel_all_ps     TRUE
  ## 2: 0.744  0.7941176       0.046       spending             none local_hommel_all_ps     TRUE
  ## 3: 0.733  1.0000000       0.022    fixed_k_adj             none local_hommel_all_ps     TRUE
  ## 4: 0.732  0.9795918       0.047 adaptive_k_adj             none  local_unadj_all_ps     TRUE
  ## 5: 0.729  1.0000000       0.004 adaptive_k_adj             none local_hommel_all_ps     TRUE

  ## So, worth exploring all alpha_methods, final_adj_methods, and local_hommel
  ## and local_unadj. Maybe also explore both with and without data splitting
  ## (idea is that one might imaging testing components of an index --- say an
  ## inde x with 10 variables in it. So we are not data splitting.)

  ## Compare to a situation with more k and more opportunities for ungating (t=.5)

  thek <- 10
  thel <- 4

  n_nodes_2 <- ((thek^(thel + 1)) - 1) / (thek - 1)
  n_nodes_2
  ## Tests at each level
  thek^seq(0, thel)
  num_leaves <- thek^thel
  num_leaves

  set.seed(123456)
  res_t_half_lst <- lapply(parms$idx, function(i) {
    x <- parms[.(i)]
    message(paste(c(i, x[1, ]), collapse = " "))
    tmp <- simulate_many_runs_DT(
      n_sim = n_sims, t = .5, k = 10, max_level = 4,
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
  res_t_half <- rbindlist(res_t_half_lst)

  save(res_t_half, res_t_some, file = "res_t_tests.rda")
  summary(res_t_half)

  ## Which approaches control the FWER?
  res_t_half_ok_fwer <- droplevels(res_t_half[false_error <= .05 + sim_err, ])

  ## Multiple possibilities:
  lapply(res_t_half_ok_fwer[, .SD, .SDcols = parm_nms], table)

  ## Which have the highest power?
  res_t_half_ok_fwer %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)

  ## local_hommel_all_ps wins in general with some local_unadj. Let's look for the set of approaches
  ## that we'd like to send to a bigger simulation

  ## What about with alpha fixed and no final adjustment? (the original idea)
  ## Below we see that we need data splitting and local hommel
  res_t_half_ok_fwer %>%
    filter(final_adj_method == "none" & alpha_method == "fixed") %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)
  ##
  ##       power leaf_power false_error alpha_method final_adj_method    local_adj_method adj_effN
  ##       <num>      <num>       <num>       <char>           <char>              <char>   <lgcl>
  ## 1: 0.7160000  1.0000000       0.000        fixed             none local_hommel_all_ps     TRUE
  ## 2: 0.4154847  0.3525777       0.005        fixed             none         local_simes     TRUE

  ## What among those without a final adjustment method but with a varying alpha and a version of sample splitting
  ## local_simes is clearly worst so exclude here
  res_t_half_ok_fwer %>%
    filter(final_adj_method == "none" & alpha_method != "fixed" & adj_effN & local_adj_method != "local_simes") %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)
  ##
  ##    power leaf_power false_error   alpha_method final_adj_method    local_adj_method adj_effN
  ##    <num>      <num>       <num>         <char>           <char>              <char>   <lgcl>
  ## 1: 0.756  1.0000000       0.000    fixed_k_adj             none local_hommel_all_ps     TRUE
  ## 2: 0.749  0.9989474       0.044 adaptive_k_adj             none  local_unadj_all_ps     TRUE
  ## 3: 0.734  0.5714286       0.001       spending             none local_hommel_all_ps     TRUE
  ## 4: 0.728        NaN       0.000 adaptive_k_adj             none local_hommel_all_ps     TRUE
  ## 5: 0.723  0.8666667       0.002      investing             none local_hommel_all_ps     TRUE

  ## What about using  a final global adjustment?

  res_t_half_ok_fwer %>%
    filter(final_adj_method != "none" & adj_effN & local_adj_method != "local_simes") %>%
    dplyr::select(one_of(c("power", "leaf_power", "false_error", parm_nms))) %>%
    arrange(desc(power), false_error)

  ## This suggests:

  ## (1) you can avoid local_adj_methods across all children if you use the adaptive_k_adj and sample splitting.
  ## (2) with no final global adjustment, you can have power with any of the above alpha methods (fixed and local_hommel), but also different alpha methods have slightly more power.
  ## (3) the final global adjustment also works with a wide variety of alpha_methods and including local adjustment methods.

  ## So let us inspect the following the best ideas for the broader sim to look at changes with k and l and t
  alpha_methods <- c("fixed", "fixed_k_adj", "adaptive_k_adj", "spending", "investing")
  final_adj_methods <- c("none", "fdr", "fwer")
  local_adj_methods <- c("local_hommel_all_ps", "local_unadj_all_ps")
  adj_effN <- c(TRUE, FALSE)
  ## and vary k, l, and t
})
