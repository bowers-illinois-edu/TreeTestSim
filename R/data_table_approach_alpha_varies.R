#' Simulate the Testing Procedure on a k-ary Tree with Optional Alpha Spending/Investing
#'
#' @description
#' This version of the simulation procedure incorporates an extra option for how the
#' significance level is allocated along the tree. The user may choose:
#'
#' \itemize{
#'   \item \code{"fixed"}: use the same significance level α at every node (as before),
#'   \item \code{"spending"}: at each branch the available α is reduced by a fixed fraction,
#'   \item \code{"investing"}: a simple α–investing scheme where rejections earn a bonus.
#' }
#'
#' In this implementation an extra column \code{alpha_alloc} is added to each node.
#' The root is allocated the full α. Then, if a node is “active” (i.e. its p-value is below
#' its allocated threshold) its children are tested at a level that depends on the chosen method:
#'
#' \itemize{
#'   \item For \code{"fixed"}, each child inherits the full α,
#'   \item For \code{"fixed_k_adj"}, each child is tested at \code{alpha/(1+alpha k)} which should be conservative
#'   \item For \code{"adaptive_k_adj"}, each child is tested at \code{parent_p + (1-parent_p)(alpha/k)} which should be conservative
#'   \item For \code{"spending"}, the children are tested at level \code{parent_alpha + max((alpha-parent_p),0)/k}
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
#' @param local_adj_p_fn Function. A function (e.g. \code{local_simes} or \code{local_hommel_all_ps} or \code{local_unadj_all_ps})
#'   that adjusts p-values at a node. The \code{local_unadj_all_ps} function does no local adjustment.
#' @param global_adj Character. The method to adjust the leaves (e.g. "hommel").
#' @param alpha_method Character. One of \code{"fixed"}, \code{"spending"}, \code{"investing"}, \code{"fixed_k_adj"}, \code{"adaptive_k_adj"}.
#' @param final_global_adj Character. One of \code{"none"}, \code{"fdr"}, \code{"fwer"}. The \code{"fdr"} and \code{"fwer"} options should only be used when \code{local_adj_p_fn} is \code{local_unadj_all_ps}.
#' @param return_details Logical. Whether to return the full simulated data.table.
#'
#' @return A list with components:
#' \describe{
#'   \item{treeDT}{The data.table augmented with columns \code{p_val} and \code{alpha_alloc}.}
#'   \item{sim_res}{A data.table summarizing error rates and discoveries.}
#' }
#'
#' @examples
#' dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.2)
#' res <- simulate_test_DT(dt,
#'   alpha = 0.05, k = 3, effN = 1000, N_total = 1000,
#'   beta_base = 0.1, alpha_method = "spending", spend_frac = 0.5
#' )
#' print(res$treeDT)
#'
#' @export
simulate_test_DT <- function(treeDT, alpha, k, effN, N_total, beta_base,
                             adj_effN = TRUE, local_adj_p_fn = local_simes, global_adj = "hommel",
                             alpha_method = "fixed", return_details = TRUE, final_global_adj = "none") {
  # Work on a copy so that the original tree is preserved.
  tree_sim <- copy(treeDT)
  tree_sim[, `:=`(p_val = NA_real_, alpha_alloc = NA_real_)]
  setkey(tree_sim, "node")

  # ---- Level 0 (the root) ----
  tree_sim[node == 1, alpha_alloc := alpha]
  tree_sim[node == 1, p_val := fifelse(
    nonnull,
    rbeta(1, beta_base, 1),
    runif(1, min = 0, max = 1)
  )]

  max_level <- max(tree_sim$level)

  # ---- Bottom-up approach for leaves (for comparison) ----
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
    global_adj_fn <- function(x) {
      p.adjust(x, method = global_adj)
    }
  }
  tree_sim[level == max_level, bottom_up_p_adj := global_adj_fn(p_sim)]
  bottom_up_false_error <- any(tree_sim[level == max_level & nonnull == FALSE, bottom_up_p_adj] <= alpha)
  bottom_up_true_discoveries <- sum(tree_sim[level == max_level & nonnull == TRUE, bottom_up_p_adj] <= alpha)
  bottom_up_power <- mean(tree_sim[level == max_level & nonnull == TRUE, bottom_up_p_adj] <= alpha)
  num_leaves <- nrow(tree_sim[level == max_level])

  # ---- Process level-by-level for the top-down procedure ----
  # We now update each branch separately so that each active parent “passes along” its own α–budget.
  for (l in 1:max_level) {
    # Identify active parents at level l-1: those that were tested and rejected at their allocated level.
    active_parents <- tree_sim[level == (l - 1) & !is.na(p_val) & (p_val <= alpha_alloc), .(node, alpha_alloc, p_val)]
    if (nrow(active_parents) == 0) break # no branch continues.

    # Process each active parent separately.
    for (i in seq_len(nrow(active_parents))) {
      parent_node <- active_parents$node[i]
      parent_alpha <- active_parents$alpha_alloc[i]
      parent_p <- active_parents$p_val[i]

      ## TODO: move this to a function
      # Determine children threshold based on chosen alpha method.
      if (alpha_method == "fixed") {
        child_threshold <- parent_alpha # no change: always use the fixed α.
      } else if (alpha_method == "spending") {
        # child_threshold <- spend_frac * parent_alpha

        ## Here if parent_p = .01 and parent_alpha=.05, we have .05-.01 = .04 to
        ## divide among the k children. Actually we just want alpha not parent_alpha
        ## otherwise the adjustment will grow throughout the tree.
        children_alpha_adjust <- max((alpha - parent_p), 0) / k

        child_threshold <- parent_alpha + children_alpha_adjust
      } else if (alpha_method == "investing") {
        ## child_threshold <- invest_frac * parent_alpha + invest_bonus
        children_alpha_adjust <- max((parent_alpha - parent_p), 0) / k
        child_threshold <- parent_alpha + children_alpha_adjust
      } else if (alpha_method == "fixed_k_adj") {
        child_threshold <- alpha / (1 + alpha * k)
      } else if (alpha_method == "adaptive_k_adj") {
        child_threshold <- parent_p + (1 - parent_p) * (alpha / k)
      } else {
        stop("alpha_method must be one of 'fixed', 'spending', or 'investing' or 'fixed_k_adj' or 'adaptive_k_adj'.")
      }

      # Get children of the current parent.
      child_rows <- tree_sim[parent == parent_node & level == l]
      if (nrow(child_rows) == 0) next

      # Adjust effective beta for current level.
      if (adj_effN) {
        effN_current <- N_total / (k^l)
        beta_eff_raw <- beta_base * sqrt(N_total / effN_current)
        cap <- 0.5
        effective_beta <- min(beta_eff_raw, cap)
      } else {
        effective_beta <- beta_base
      }

      # For children: if p_sim is not yet simulated, simulate it using parent's p-value.
      ## This also implements monotonicity
      child_rows[is.na(p_sim), p_sim := fifelse(
        nonnull,
        parent_p + (1 - parent_p) * rbeta(.N, effective_beta, 1),
        runif(.N, min = parent_p, max = 1)
      )]

      ## p_sim is like the temporary p_value. It exists here because we did the
      ## bottom up tests first on all of the leaves

      child_rows[is.na(p_val), p_val := p_sim]

      # Locally adjust the group of children using the specified local adjustment function.
      child_rows[, local_adj_p := local_adj_p_fn(p_val), by = parent]

      # Gate the branch: if the locally adjusted p-value exceeds the branch’s threshold,
      # set p_val to NA (i.e. stop testing that branch).
      child_rows[local_adj_p > child_threshold, p_val := NA_real_]
      # Otherwise, take the maximum of local_adj_p and the original p_val.
      child_rows[local_adj_p <= child_threshold, p_val := pmax(local_adj_p, p_val)]

      # Update the α–budget for these children.
      child_rows[, alpha_alloc := child_threshold]

      # Write the updated values back into tree_sim.
      tree_sim[node %in% child_rows$node, `:=`(
        p_val = child_rows$p_val,
        alpha_alloc = child_rows$alpha_alloc,
        local_adj_p = child_rows$local_adj_p
      )]
    } # end loop over active parents at level l
  } # end levels

  ## This next does a final global adjustment given the p-values that survive the gating and monotonicity
  if (final_global_adj == "fdr") {
    tree_sim[!is.na(p_val), p_val_final_adj := p.adjust(p_val, method = "BH")]
  } else if (final_global_adj == "fwer") {
    ## This next might be slower than the hommel() function
    ## tree_sim[, p_val_final_adj := p.adjust(p_val,method="hommel")
    tree_sim[!is.na(p_val), p_val_final_adj := hommel(p_val)@adjusted]
  } else if (final_global_adj == "none") {
    ## No adjustment
    tree_sim[, p_val_final_adj := p_val]
  } else {
    stop("final_global_adj must be one of 'fdr', 'fwer', or 'none'")
  }

  # ---- Summarize errors and discoveries ----
  num_nodes_tested <- sum(!is.na(tree_sim$p_val))
  num_nonnull_nodes_tested <- sum(!is.na(tree_sim$p_val) & tree_sim$nonnull == TRUE)
  false_error <- any(tree_sim[nonnull == FALSE & !is.na(p_val), p_val] <=
    tree_sim[nonnull == FALSE & !is.na(p_val), alpha_alloc])
  true_discoveries <- sum(tree_sim[nonnull == TRUE & !is.na(p_val), p_val] <=
    tree_sim[nonnull == TRUE & !is.na(p_val), alpha_alloc])
  power <- mean(tree_sim[nonnull == TRUE & !is.na(p_val), p_val] <=
    tree_sim[nonnull == TRUE & !is.na(p_val), alpha_alloc])
  ## This is a test
  num_leaves_tested <- sum(tree_sim[level == max_level, !is.na(p_val)])
  stopifnot(num_leaves_tested <= num_leaves)
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

#' Simulate the Family-Wise Error Rate (FWER) using a Data.table-based Tree with Alpha Options
#'
#' @description
#' Runs the hierarchical testing procedure on complete k-ary trees (generated via \code{generate_tree_DT()})
#' \code{n_sim} times and returns the estimated FWER. The user may choose whether to use a fixed α,
#' α–spending, or a rudimentary α–investing scheme.
#'
#' @param n_sim Integer. The number of simulation replicates.
#' @param t Numeric in [0,1]. The probability that a leaf is non-null.
#' @param k Integer. The branching factor.
#' @param max_level Integer. The maximum level (depth) of the tree.
#' @param alpha Numeric. The nominal significance level.
#' @param N_total Numeric. The total sample size at the root.
#' @param beta_base Numeric. The base parameter for the Beta distribution.
#' @param return_details Logical.
#' @param alpha_method Character. One of \code{"fixed"}, \code{"spending"}, or \code{"investing"}.
#' @param spend_frac Numeric. (For spending) Fraction of parent's α allocated to children.
#' @param invest_frac Numeric. (For investing) Fraction of parent's α allocated to children.
#' @param invest_bonus Numeric. (For investing) Bonus added to the allocation.
#' @param ... Additional arguments passed to \code{simulate_test_DT()}.
#'
#' @return A numeric value in [0,1] giving the estimated FWER (plus additional simulation summary metrics).
#'
#' @examples
#' fwer <- simulate_many_runs_DT(
#'   n_sim = 1000, t = 0.2, k = 3, max_level = 3,
#'   alpha = 0.05, N_total = 1000, beta_base = 0.1,
#'   alpha_method = "investing", invest_frac = 0.5, invest_bonus = 0.01
#' )
#' print(fwer)
#'
#' @importFrom parallel mclapply
#' @export
simulate_many_runs_DT <- function(n_sim, t, k, max_level, alpha, N_total, beta_base = 0.1,
                                  adj_effN = TRUE, local_adj_p_fn = local_simes, return_details = FALSE,
                                  global_adj = "hommel", alpha_method = "fixed", final_global_adj = "none", multicore = FALSE) {
  treeDT <- generate_tree_DT(max_level, k, t)

  if (multicore) {
    ncores <- future::availableCores()
  } else {
    ncores <- 1
  }

  res <- mclapply(
    seq_len(n_sim),
    function(i) {
      simulate_test_DT(treeDT, alpha, k,
        effN = N_total, N_total = N_total, beta_base = beta_base,
        adj_effN = adj_effN, local_adj_p_fn = local_adj_p_fn, global_adj = global_adj,
        alpha_method = alpha_method,
        return_details = return_details, final_global_adj = final_global_adj
      )
    },
    mc.cores = ncores, mc.set.seed = TRUE
  )

  ## res_dt <- rbindlist(lapply(res, function(x) x$sim_res))
  res_dt <- rbindlist(res)
  mean_res <- unlist(res_dt[, lapply(.SD, mean, na.rm = TRUE)])
  return(mean_res)
}


#' Generate a Complete k-ary Tree as a data.table
#'
#' @description
#' Constructs a complete k-ary tree with levels 0 to \code{max_level}.
#' The tree is stored in a data.table with one row per node and columns for:
#' \code{node} (the node's index),
#' \code{level} (its level in the tree, with 0 = root),
#' \code{parent} (the index of its parent; NA for the root),
#' and \code{nonnull} (a logical indicating whether the node is non-null).
#'
#' At the leaves (nodes at \code{level == max_level}), the \code{nonnull} flag
#' is sampled as \code{TRUE} with probability \code{t}. For internal nodes,
#' \code{nonnull} is set to \code{TRUE} if any of its children are non-null.
#'
#' @param max_level Integer. The maximum level of the tree (with the root at level 0).
#' @param k Integer. The branching factor (each internal node has exactly \code{k} children).
#' @param t Numeric in [0,1]. The probability that a leaf is non-null (has a treatment effect).
#'
#' @return A data.table with columns \code{node}, \code{level}, \code{parent}, and \code{nonnull}.
#'
#' @examples
#' dt <- generate_tree_DT(max_level = 3, k = 3, t = 0.2)
#' print(dt)
#'
#' @export
generate_tree_DT <- function(max_level, k, t) {
  # Total number of nodes in a complete k-ary tree:
  total_nodes <- sum(k^(0:max_level))

  # Initialize the data.table with node IDs.
  treeDT <- data.table(node = 1:total_nodes)

  # Compute levels by assigning nodes in breadth-first order.
  levels <- integer(total_nodes)
  cum <- cumsum(k^(0:max_level))
  start <- 1
  for (l in 0:max_level) {
    end <- cum[l + 1]
    levels[start:end] <- l
    start <- end + 1
  }
  treeDT[, level := levels]

  # Compute parent indices.
  # In a complete k-ary tree stored in breadth-first order:
  # For node i (i > 1), parent = floor((i - 2) / k) + 1.
  treeDT[, parent := ifelse(node == 1, NA_integer_, floor((node - 2) / k) + 1)]

  # Allocate the "nonnull" flag.
  # For leaves (level == max_level), sample from Bernoulli(t).
  ## treeDT[level == max_level, nonnull := runif(.N) < t]
  # Complete Random assignment at leaves
  treeDT[, nonnull := NA]
  treeDT[level == max_level, nonnull := {
    n_leaves <- .N
    m_leaves <- floor(t * n_leaves)
    sample(rep(c(FALSE, TRUE), c(n_leaves - m_leaves, m_leaves)))
  }]

  # Propagate nonnull flags upward (for internal nodes).
  if (max_level > 0) {
    for (l in (max_level - 1):0) {
      # For nodes at level l, set nonnull = TRUE if any of its children are TRUE.
      treeDT[level == l, nonnull := {
        sapply(node, function(i) {
          # Compute children indices for node i:
          children <- ((i - 1) * k + 2):((i - 1) * k + k + 1)
          children <- children[children <= total_nodes]
          any(treeDT[children, nonnull], na.rm = TRUE)
        })
      }]
    }
  } else {
    # Only one node.
    treeDT[node == 1, nonnull := runif(1) < t]
  }

  return(treeDT)
}

#' @title Compute local Simes p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the Simes p-value
#'   \eqn{\min_{i=1\ldots k} \{ (k/i) * p_{(i)} \}}, where \eqn{p_{(1)} \le \ldots \le p_{(k)}}.
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used in this function)
#'
#' @details The Simes p-value is a valid test for the intersection hypothesis
#'   \eqn{H_{\cap} : \text{all child hypotheses are null}}, under certain independence
#'   or positive dependence conditions among the p-values.
#'
#' @return A single numeric value: the Simes combination p-value.
#'
#' @examples
#' local_simes(c(0.01, 0.04, 0.10, 0.20))
#' @export
local_simes <- function(pvals_children, alpha = .05) {
  k <- length(pvals_children)
  sort_p <- sort(pvals_children)
  i_seq <- seq_len(k)
  simes_vals <- (k / i_seq) * sort_p
  return(min(simes_vals))
}

#' @title Compute local Hommel p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the Hommell adjusted p-values
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used in this function)
#'
#' @details TODO The Hommel adjustment and the hommel package is faster than p.adjust().
#'
#' @return A vector of adjusted p-values
#'
#' @examples
#' local_simes(c(0.01, 0.04, 0.10, 0.20))
#'
#' @importFrom hommel hommel
#' @export
local_hommel_all_ps <- function(pvals_children, alpha = alpha) {
  adj_p_vals <- hommel(pvals_children)@adjusted
  return(adj_p_vals)
}

#' @title Unadjusted local step
#'
#' @description Given \eqn{k} child p-values, return the highest p-value below alpha. This is the same as respecting the global gating criteria of not continuing to test along any path where p > alpha, but doesn't do a local, parent-level adjustment. This is mostly useful for simulation and study purposes although with k=2 and low l  (for example, less than 5 or so) this will also tend to control the FWER.
#'
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha
#'
#' @details The idea is to not adjust so that we can compare the consequences of not adjusting to adjusting
#'
#' @return A single numeric value: the highest p<=alpha among the children of the parent node or the minimum of the p-values of the children if they are all greater than alpha
#'
#' @examples
#' local_min_p(c(0.01, 0.04, 0.10, 0.20)) ## should return .04
#' local_min_p(c(0.10, 0.20)) ## should return .04
#'
#' @export
local_min_p <- function(pvals_children, alpha = .05) {
  p_le_alpha <- pvals_children[pvals_children <= alpha]
  if (length(p_le_alpha) == 0) {
    return(min(pvals_children))
  } else {
    return(max(p_le_alpha))
  }
}

#' @title Unadjusted local step
#'
#' @description Given \eqn{k} child p-values just pass them through for simulation purposes
#'
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha
#'
#' @details The idea is to not adjust so that we can compare the consequences of not adjusting to adjusting
#'
#' @return A single numeric value: just passes through the unadjusted p values
#'
#' @export
local_unadj_all_ps <- function(pvals_children, alpha = .05) {
  pvals_children
}


#' @title Compute local BH FDR adjusted p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the BH/FDR adjusted p-values
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used in this function)
#'
#' @details TODO
#'
#' @return A vector of adjusted p-values
#'
#' @examples
#' local_bh_ps(c(0.01, 0.04, 0.10, 0.20))
#'
#' @importFrom hommel hommel
#' @export
local_bh_all_ps <- function(pvals_children, alpha = alpha) {
  adj_p_vals <- p.adjust(pvals_children, method = "BH")
  return(adj_p_vals)
}


#' Sample from a Scaled Beta Distribution
#'
#' This function generates random samples from a Beta distribution and then linearly transforms
#' the samples so that they lie in the interval \eqn{[a, 1]}, where \eqn{a} is a user-defined minimum.
#'
#' @param a Numeric. The lower bound of the output interval. Must be between 0 (inclusive) and 1 (exclusive).
#' @param shape1 Numeric. The first shape parameter of the Beta distribution.
#' @param shape2 Numeric. The second shape parameter of the Beta distribution.
#' @param n Integer. The number of random samples to generate (default is 1).
#'
#' @return A numeric vector containing \eqn{n} samples that have been scaled to lie in the interval \eqn{[a, 1]}.
#'
#' @details The function first generates \eqn{n} samples from a standard Beta distribution on the interval \eqn{[0, 1]}
#' using \code{rbeta(n, shape1, shape2)}. It then performs a linear transformation of the form:
#' \deqn{x = a + (1 - a) \times \text{beta_sample}}
#' so that when \code{beta_sample} equals 0, the transformed value is \eqn{a}, and when \code{beta_sample} equals 1,
#' the output is 1.
#'
#' @examples
#' # Generate 10 samples with a lower bound of 0.2, and Beta parameters shape1 = 2 and shape2 = 5:
#' samples <- r_scaled_beta(0.2, 2, 5, n = 10)
#' print(samples)
#'
#' @export
r_scaled_beta <- function(a, shape1, shape2, n = 1) {
  # Ensure 'a' is a numeric value within the interval [0, 1).
  if (!is.numeric(a) || a < 0 || a >= 1) {
    stop("Parameter 'a' must be numeric and in the interval [0, 1).")
  }

  # Generate 'n' samples from the standard Beta distribution on [0, 1].
  beta_samples <- rbeta(n, shape1, shape2)

  # Linearly transform the Beta samples to the interval [a, 1].
  # The transformation is: scaled_sample = a + (1 - a) * beta_sample.
  # When beta_sample = 0, scaled_sample = a, and when beta_sample = 1, scaled_sample = 1.
  scaled_samples <- a + (1 - a) * beta_samples

  # Return the transformed samples.
  return(scaled_samples)
}
