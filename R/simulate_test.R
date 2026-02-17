#' Convert Cohen's d to a Beta distribution shape parameter
#'
#' Given an effect size (Cohen's d), a sample size, and a significance level,
#' computes the Beta(a, 1) shape parameter \code{a} such that
#' \code{P(Beta(a,1) < alpha) = power}, where power comes from the normal
#' approximation for a two-sample test.
#'
#' @param effect_size Numeric. Cohen's d at each non-null leaf.
#' @param N_level Numeric (scalar or vector). Sample size at the relevant
#'   tree level.
#' @param alpha Numeric. The nominal significance level.
#' @return A numeric vector of Beta shape parameters (one per element of
#'   \code{N_level}).
#' @keywords internal
effect_size_to_beta <- function(effect_size, N_level, alpha = 0.05) {
  z_crit <- qnorm(1 - alpha / 2)
  power <- pnorm(effect_size * sqrt(N_level / 4) - z_crit)
  # Clamp: power can't drop below alpha (would give beta >= 1, no signal)
  # and can't hit exactly 1 (would give beta = 0)
  power <- pmax(power, alpha + 1e-6)
  power <- pmin(power, 1 - 1e-10)
  log(power) / log(alpha)
}

#' Simulate the Testing Procedure on a k-ary Tree with Optional Alpha Spending/Investing
#'
#' @description
#' This version of the simulation procedure incorporates an extra option for how the
#' significance level is allocated along the tree. The user may choose:
#' \itemize{
#'   \item \code{"fixed"}: use the same significance level \eqn{\alpha} at every node (as before),
#'   \item \code{"spending"}: at each branch the available \eqn{\alpha} is reduced by a fixed fraction,
#'   \item \code{"investing"}: a simple \eqn{\alpha}-investing scheme where rejections earn a bonus.
#' }
#'
#' In this implementation an extra column \code{alpha_alloc} is added to each node.
#' The root is allocated the full \eqn{\alpha}. Then, if a node is "active" (i.e. its p-value is below
#' its allocated threshold) its children are tested at a level that depends on the chosen method:
#' \itemize{
#'   \item For \code{"fixed"}, each child inherits the full \eqn{\alpha},
#'   \item For \code{"fixed_k_adj"}, each child is tested at \code{alpha/(1+alpha k)} which should be conservative,
#'   \item For \code{"adaptive_k_adj"}, each child is tested at \code{parent_p + (1-parent_p)(alpha/k)} which should be conservative,
#'   \item For \code{"spending"}, the children are tested at level \code{parent_alpha + max((alpha-parent_p),0)/k},
#'   \item For \code{"investing"}, the children are tested at level \code{parent_alpha + max((parent_alpha-parent_p),0)/k}.
#' }
#'
#' @param treeDT A data.table as produced by \code{generate_tree_DT()}.
#' @param alpha Numeric. The nominal significance level.
#' @param k Integer. The branching factor.
#' @param N_total Numeric. The total sample size at the root.
#' @param effect_size Numeric. Cohen's d — the standardized effect size at
#'   each non-null leaf. Together with \code{N_total} and the tree structure,
#'   this determines power at every level of the tree.
#' @param power_decay Logical. If \code{TRUE} (default), the strength of
#'   non-null p-value generation decreases at deeper levels (reflecting the
#'   smaller per-node sample size \code{N_total / k^level}). If \code{FALSE},
#'   the root-level power is used at all depths — useful for modeling
#'   situations where the effect size itself grows to compensate for smaller
#'   samples.
#' @param local_adj_p_fn Function. A function that adjusts p-values at a node (e.g. \code{local_simes}).
#' @param global_adj Character. The method to adjust the leaves (e.g. "hommel").
#' @param alpha_method Character. One of \code{"fixed"}, \code{"spending"},
#'   \code{"investing"}, \code{"fixed_k_adj"}, \code{"adaptive_k_adj"}, or
#'   \code{"adaptive_power"}. The \code{"adaptive_power"} method uses
#'   \code{\link[manytestsr]{compute_adaptive_alphas}} to precompute a
#'   depth-indexed alpha schedule based on estimated power decay through the
#'   tree.
#' @param return_details Logical. Whether to return the full simulated data.table.
#' @param final_global_adj Character. One of \code{"none"}, \code{"fdr"}, \code{"fwer"}.
#' @param monotonicity Logical. TRUE if we require child nodes p-values to be
#'   no larger than those of parent nodes. FALSE child nodes could have smaller
#'   p-values.
#'
#' @return A list with components:
#' \describe{
#'   \item{treeDT}{The data.table augmented with columns \code{p_val} and \code{alpha_alloc}.}
#'   \item{sim_res}{A data.table summarizing error rates, discoveries, and
#'     \code{root_power} (the theoretical power at the root level).}
#' }
#' @examples
#' dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.2)
#' res <- simulate_test_DT(dt,
#'   alpha = 0.05, k = 3, N_total = 1000,
#'   effect_size = 0.5, alpha_method = "spending", return_details = TRUE
#' )
#' @import data.table
#' @importFrom stats pnorm qnorm
#' @export
simulate_test_DT <- function(treeDT, alpha, k, N_total, effect_size,
                             power_decay = TRUE, local_adj_p_fn = local_simes,
                             global_adj = "hommel",
                             alpha_method = "fixed", return_details = TRUE,
                             final_global_adj = "none",
                             monotonicity = TRUE) {
  tree_sim <- copy(treeDT)
  tree_sim[, `:=`(p_val = NA_real_, alpha_alloc = NA_real_)]
  setkey(tree_sim, "node")

  ## Root p-value: beta parameter from effect_size at full sample size
  root_beta <- effect_size_to_beta(effect_size, N_total, alpha)

  tree_sim[node == 1, `:=`(
    alpha_alloc = alpha,
    p_val = fifelse(
      nonnull,
      rbeta(1, root_beta, 1),
      runif(1, min = 0, max = 1)
    )
  )]

  max_level <- max(tree_sim$level)

  # Bottom-up approach for comparison: use leaf-level sample size
  eff_beta_leaves <- effect_size_to_beta(effect_size, N_total / (k^max_level), alpha)
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

  # Dual bottom-up: always compute both hommel and BH on leaf p_sim
  tree_sim[level == max_level, bu_hommel_p := local_hommel_all_ps(p_sim)]
  bu_hommel_false_error <- any(tree_sim[level == max_level & nonnull == FALSE, bu_hommel_p] <= alpha)
  bu_hommel_true_disc <- sum(tree_sim[level == max_level & nonnull == TRUE, bu_hommel_p] <= alpha)
  bu_hommel_power <- mean(tree_sim[level == max_level & nonnull == TRUE, bu_hommel_p] <= alpha)

  tree_sim[level == max_level, bu_bh_p := p.adjust(p_sim, method = "BH")]
  bu_bh_false_error <- any(tree_sim[level == max_level & nonnull == FALSE, bu_bh_p] <= alpha)
  bu_bh_true_disc <- sum(tree_sim[level == max_level & nonnull == TRUE, bu_bh_p] <= alpha)
  bu_bh_power <- mean(tree_sim[level == max_level & nonnull == TRUE, bu_bh_p] <= alpha)

  # Precompute adaptive alpha schedule if needed
  if (alpha_method %in% c("adaptive_power", "adaptive_power_pruned")) {
    alpha_schedule <- manytestsr::compute_adaptive_alphas(
      k = k, delta_hat = effect_size / 2, N_total = N_total,
      max_depth = max_level + 1L, thealpha = alpha
    )
  }

  # For pruned method, pre-compute per-depth path power (product of ancestor
  # thetas). These depend only on effect size and sample sizes, not on which
  # branches survive, so we compute them once.
  if (alpha_method == "adaptive_power_pruned") {
    z_crit <- qnorm(1 - alpha / 2)
    delta_hat <- effect_size / 2
    path_power_by_depth <- numeric(max_level + 1L)
    path_power_by_depth[1] <- 1.0  # root has no ancestors
    for (d in 2:(max_level + 1L)) {
      theta_prev <- pnorm(delta_hat * sqrt(N_total / k^(d - 2)) - z_crit)
      theta_prev <- max(min(theta_prev, 1), 0)
      path_power_by_depth[d] <- path_power_by_depth[d - 1] * theta_prev
    }
  }

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
      } else if (alpha_method == "adaptive_power") {
        # Tree levels are 0-indexed (root = 0); compute_adaptive_alphas is
        # 1-indexed (root = 1). Children at tree level l -> depth l + 1.
        child_threshold <- alpha_schedule[l + 1]
      } else if (alpha_method == "adaptive_power_pruned") {
        # Same lookup as adaptive_power; the schedule is updated after
        # each level to reflect branch pruning (see below).
        child_threshold <- alpha_schedule[l + 1]
      } else {
        stop("alpha_method must be one of 'fixed', 'spending', 'investing', 'fixed_k_adj', 'adaptive_k_adj', 'adaptive_power', or 'adaptive_power_pruned'.")
      }

      child_rows <- tree_sim[parent == parent_node & level == l]
      if (nrow(child_rows) == 0) next

      ## Power-calibrated beta at this level
      if (power_decay) {
        effective_beta <- effect_size_to_beta(effect_size, N_total / (k^l), alpha)
      } else {
        ## No decay: use root-level power at all depths
        effective_beta <- root_beta
      }

      if (monotonicity) {
        ## Child p-values are at least as large as the parent's
        child_rows[is.na(p_sim), p_sim := fifelse(
          nonnull,
          parent_p + (1 - parent_p) * rbeta(.N, effective_beta, 1),
          runif(.N, min = parent_p, max = 1)
        )]
      } else {
        parent_p_tmp <- 0
        child_rows[is.na(p_sim), p_sim := fifelse(
          nonnull,
          parent_p_tmp + (1 - parent_p_tmp) * rbeta(.N, effective_beta, 1),
          runif(.N, min = parent_p_tmp, max = 1)
        )]
      }

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

    # Branch pruning: after all parents at level l are processed, count
    # how many children survived and recompute the alpha schedule for
    # deeper levels. Fewer surviving nodes = lower error load = more
    # alpha budget per test at subsequent depths.
    if (alpha_method == "adaptive_power_pruned" && l < max_level) {
      n_survived <- tree_sim[level == l & !is.na(p_val) & (p_val <= alpha_alloc), .N]
      if (n_survived > 0L) {
        for (d in (l + 2):(max_level + 1L)) {
          # In pruned tree: n_survived * k^{d-1-l} nodes at depth d
          # (vs k^{d-1} in full tree). path_power is unchanged.
          n_nodes_pruned <- n_survived * k^(d - 1L - l)
          sum_path_power <- n_nodes_pruned * path_power_by_depth[d]
          alpha_schedule[d] <- min(alpha, alpha / sum_path_power)
        }
      }
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
  num_leaves_tested <- sum(tree_sim[level == max_level, !is.na(p_val)])
  if (num_leaves_tested > 0) {
    leaf_power <- mean(tree_sim[level == max_level & nonnull == TRUE & !is.na(p_val), p_val] <= alpha, na.rm = TRUE)
    leaf_disc <- sum(tree_sim[level == max_level & nonnull == TRUE & !is.na(p_val), p_val] <= alpha, na.rm = TRUE)
  } else {
    leaf_power <- NA
    leaf_disc <- NA
  }

  # Theoretical power at the root level
  root_power <- pnorm(effect_size * sqrt(N_total / 4) - qnorm(1 - alpha / 2))

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
    num_leaves = num_leaves,
    bu_hommel_false_error = bu_hommel_false_error,
    bu_hommel_true_disc = bu_hommel_true_disc,
    bu_hommel_power = bu_hommel_power,
    bu_bh_false_error = bu_bh_false_error,
    bu_bh_true_disc = bu_bh_true_disc,
    bu_bh_power = bu_bh_power,
    root_power = root_power
  )

  if (return_details) {
    return(list(treeDT = tree_sim, sim_res = sim_res))
  } else {
    return(sim_res)
  }
}
