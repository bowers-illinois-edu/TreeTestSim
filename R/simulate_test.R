#' Simulate the Testing Procedure on a k-ary Tree with Optional Alpha Spending/Investing
#'
#' @description
#' This version of the simulation procedure incorporates an extra option for how the
#' significance level is allocated along the tree. The user may choose:
#' \itemize{
#'   \item \code{"fixed"}: use the same significance level α at every node (as before),
#'   \item \code{"spending"}: at each branch the available α is reduced by a fixed fraction,
#'   \item \code{"investing"}: a simple α–investing scheme where rejections earn a bonus.
#' }
#'
#' In this implementation an extra column \code{alpha_alloc} is added to each node.
#' The root is allocated the full α. Then, if a node is “active” (i.e. its p-value is below
#' its allocated threshold) its children are tested at a level that depends on the chosen method:
#' \itemize{
#'   \item For \code{"fixed"}, each child inherits the full α,
#'   \item For \code{"fixed_k_adj"}, each child is tested at \code{alpha/(1+alpha k)} which should be conservative,
#'   \item For \code{"adaptive_k_adj"}, each child is tested at \code{parent_p + (1-parent_p)(alpha/k)} which should be conservative,
#'   \item For \code{"spending"}, the children are tested at level \code{parent_alpha + max((alpha-parent_p),0)/k},
#'   \item For \code{"investing"}, the children are tested at level \code{parent_alpha + max((parent_alpha-parent_p),0)/k}.
#' }
#'
#' @param treeDT A data.table as produced by \code{generate_tree_DT()}.
#' @param alpha Numeric. The nominal significance level.
#' @param k Integer. The branching factor.
#' @param effN Numeric. The effective sample size at the root.
#' @param N_total Numeric. The total sample size at the root.
#' @param beta_base Numeric. The base parameter for the beta distribution.
#' @param adj_effN Logical. Whether to adjust the effective sample size at deeper levels.
#' @param local_adj_p_fn Function. A function that adjusts p-values at a node (e.g. \code{local_simes}).
#' @param global_adj Character. The method to adjust the leaves (e.g. "hommel").
#' @param alpha_method Character. One of \code{"fixed"}, \code{"spending"}, \code{"investing"}, \code{"fixed_k_adj"}, \code{"adaptive_k_adj"}.
#' @param return_details Logical. Whether to return the full simulated data.table.
#' @param final_global_adj Character. One of \code{"none"}, \code{"fdr"}, \code{"fwer"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{treeDT}{The data.table augmented with columns \code{p_val} and \code{alpha_alloc}.}
#'   \item{sim_res}{A data.table summarizing error rates and discoveries.}
#' }
#' @examples
#' dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.2)
#' res <- simulate_test_DT(dt,
#'   alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
#'   beta_base = 0.1, alpha_method = "spending", return_details = TRUE
#' )
#' @import data.table
#' @export
simulate_test_DT <- function(treeDT, alpha, k, effN, N_total, beta_base,
                             adj_effN = TRUE, local_adj_p_fn = local_simes, global_adj = "hommel",
                             alpha_method = "fixed", return_details = TRUE, final_global_adj = "none") {
  tree_sim <- copy(treeDT)
  tree_sim[, `:=`(p_val = NA_real_, alpha_alloc = NA_real_)]
  setkey(tree_sim, "node")

  tree_sim[node == 1, `:=`(
    alpha_alloc = alpha,
    p_val = fifelse(
      nonnull,
      rbeta(1, beta_base, 1),
      runif(1, min = 0, max = 1)
    )
  )]

  max_level <- max(tree_sim$level)

  # Bottom-up approach for comparison

  ### TODO: So we should set minimum power for each leaf (like the power you'd
  ### get with N_total/(k^max_level))
  ## And then maxpower at the overall test which is a function of
  ## pbeta(.05,a,1). Say, we set overall power at pbeta(alpha,a,1) choosing an
  ## a such that this equals .8, then we ask, say, what the N would have to be
  ## for this for a t-test maybe??? Basically, we can have minimum power (say,
  ## imagine relatively large leaves, say N=50, with 25 treated and 25
  ## controls, so pbeta(.05,a,1)==blah where blah is the power of a t.test like
  ## power.t.test()). For example, power.t.test(delta=.8,sd=1,sig.level=.05,n=25) -> power=.8
  ## imagine for the sake of these simulations that no block has fewer than that.
  ## So we get num_leaves*(25*2) as the total size of the N in the dataset. And this determines the top level power.
  ## At each split we divide the N_node/k -> convert to t.test power
  ## Effect at each node is the size weighted average of the block averages of that block (for delta, assume sd=1).



  effN_leaves <- N_total / (k^max_level)
  beta_eff_leaves_raw <- beta_base * sqrt(N_total / effN_leaves)
  cap_leaves <- 0.5
  eff_beta_leaves <- min(cap_leaves, beta_eff_leaves_raw)
  tree_sim[level == max_level, p_sim := fifelse(
    nonnull,
    rbeta(.N, eff_beta_leaves, 1),
    runif(.N, min = 0, max = 1)
  )]
  if (global_adj == "hommel") {
    global_adj_fn <- local_hommel_all_ps
  } else {
    global_adj_fn <- function(x) p.adjust(x, method = global_adj)
  }
  tree_sim[level == max_level, bottom_up_p_adj := global_adj_fn(p_sim)]
  bottom_up_false_error <- any(tree_sim[level == max_level & nonnull == FALSE, bottom_up_p_adj] <= alpha)
  bottom_up_true_discoveries <- sum(tree_sim[level == max_level & nonnull == TRUE, bottom_up_p_adj] <= alpha)
  bottom_up_power <- mean(tree_sim[level == max_level & nonnull == TRUE, bottom_up_p_adj] <= alpha)
  num_leaves <- nrow(tree_sim[level == max_level])

  # Top-down procedure
  for (l in 1:max_level) {
    active_parents <- tree_sim[level == (l - 1) & !is.na(p_val) & (p_val <= alpha_alloc), .(node, alpha_alloc, p_val)]
    if (nrow(active_parents) == 0) {
      break
    }

    for (i in seq_len(nrow(active_parents))) {
      parent_node <- active_parents$node[i]
      parent_alpha <- active_parents$alpha_alloc[i]
      parent_p <- active_parents$p_val[i]

      if (alpha_method == "fixed") {
        child_threshold <- parent_alpha
      } else if (alpha_method == "spending") {
        children_alpha_adjust <- max((alpha - parent_p), 0) / k
        child_threshold <- parent_alpha + children_alpha_adjust
      } else if (alpha_method == "investing") {
        children_alpha_adjust <- max((parent_alpha - parent_p), 0) / k
        child_threshold <- parent_alpha + children_alpha_adjust
      } else if (alpha_method == "fixed_k_adj") {
        child_threshold <- alpha / (1 + alpha * k)
      } else if (alpha_method == "adaptive_k_adj") {
        child_threshold <- parent_p + (1 - parent_p) * (alpha / k)
      } else {
        stop("alpha_method must be one of 'fixed', 'spending', 'investing', 'fixed_k_adj', or 'adaptive_k_adj'.")
      }

      child_rows <- tree_sim[parent == parent_node & level == l]
      if (nrow(child_rows) == 0) next

      if (adj_effN) {
        effN_current <- N_total / (k^l)
        beta_eff_raw <- beta_base * sqrt(N_total / effN_current)
        cap <- 0.5
        effective_beta <- min(beta_eff_raw, cap)
      } else {
        effective_beta <- beta_base
      }

      child_rows[is.na(p_sim), p_sim := fifelse(
        nonnull,
        parent_p + (1 - parent_p) * rbeta(.N, effective_beta, 1),
        runif(.N, min = parent_p, max = 1)
      )]

      child_rows[is.na(p_val), p_val := p_sim]

      child_rows[, local_adj_p := local_adj_p_fn(p_val), by = parent]

      child_rows[local_adj_p > child_threshold, p_val := NA_real_]
      child_rows[local_adj_p <= child_threshold, p_val := pmax(local_adj_p, p_val)]
      child_rows[, alpha_alloc := child_threshold]

      tree_sim[node %in% child_rows$node, `:=`(
        p_val = child_rows$p_val,
        alpha_alloc = child_rows$alpha_alloc
      )]
    }
  }

  if (final_global_adj == "fdr") {
    tree_sim[!is.na(p_val), p_val_final_adj := p.adjust(p_val, method = "BH")]
  } else if (final_global_adj == "fwer") {
    tree_sim[!is.na(p_val), p_val_final_adj := hommel(p_val)@adjusted]
  } else if (final_global_adj == "none") {
    tree_sim[, p_val_final_adj := p_val]
  } else {
    stop("final_global_adj must be one of 'fdr', 'fwer', or 'none'")
  }

  num_nodes_tested <- sum(!is.na(tree_sim$p_val))
  num_nonnull_nodes_tested <- sum(!is.na(tree_sim$p_val) & tree_sim$nonnull)
  false_error <- any(tree_sim[nonnull == FALSE & !is.na(p_val), p_val] <=
    tree_sim[nonnull == FALSE & !is.na(p_val), alpha_alloc])
  true_discoveries <- sum(tree_sim[nonnull == TRUE & !is.na(p_val), p_val] <=
    tree_sim[nonnull == TRUE & !is.na(p_val), alpha_alloc])
  power <- mean(tree_sim[nonnull == TRUE & !is.na(p_val), p_val] <=
    tree_sim[nonnull == TRUE & !is.na(p_val), alpha_alloc])
  ## This is a test
  num_leaves_tested <- sum(tree_sim[level == max_level, !is.na(p_val)])
  if (num_leaves_tested > 0) {
    leaf_power <- mean(tree_sim[level == max_level & nonnull == TRUE & !is.na(p_val), p_val] <= alpha, na.rm = TRUE)
    leaf_disc <- sum(tree_sim[level == max_level & nonnull == TRUE & !is.na(p_val), p_val] <= alpha, na.rm = TRUE)
  } else {
    leaf_power <- NA
    leaf_disc <- NA
  }

  sim_res <- data.table(
    num_nodes_tested = num_nodes_tested,
    num_nonnull_nodes_tested = num_nonnull_nodes_tested,
    false_error = false_error,
    true_discoveries = true_discoveries,
    power = power,
    num_leaves_tested = num_leaves_tested,
    leaf_power = leaf_power,
    leaf_disc = leaf_disc,
    bottom_up_false_error = bottom_up_false_error,
    bottom_up_true_discoveries = bottom_up_true_discoveries,
    bottom_up_power = bottom_up_power,
    num_leaves = num_leaves
  )

  if (return_details) {
    return(list(treeDT = tree_sim, sim_res = sim_res))
  } else {
    return(sim_res)
  }
}
