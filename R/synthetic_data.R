#' Generate a Synthetic Block-Randomized Experiment with k-ary Tree Structure
#'
#' Creates individual-level and block-level data for a complete k-ary tree
#' experiment. Leaf nodes are blocks with \code{N_per_block} observations each,
#' block-randomized to treatment and control. The block hierarchy is encoded in
#' a factor variable (\code{lvls_fac}) suitable for use with
#' \code{\link[manytestsr]{splitSpecifiedFactorMulti}}.
#'
#' @param k Integer. Branching factor (number of children per internal node).
#' @param max_level Integer. Tree depth (number of levels below the root).
#'   Total number of leaf blocks is \code{k^max_level}.
#' @param N_per_block Integer. Number of observations per leaf block
#'   (default 50). Must be even for balanced randomization.
#' @param t Numeric between 0 and 1. Probability that a leaf is non-null (has a
#'   treatment effect). Propagated bottom-up: internal nodes are non-null if
#'   any descendant leaf is non-null.
#' @param block_sd Numeric. Standard deviation of block-level baseline
#'   means (default 0.5). Individual outcomes are
#'   \code{y0 ~ N(block_mean, 1)}.
#'
#' @return A list with components:
#' \describe{
#'   \item{idat}{Individual-level data.table with columns: \code{bF} (block
#'     factor), \code{id}, \code{y0} (potential outcome under control),
#'     \code{trt} (0/1 treatment), \code{trtF} (factor treatment),
#'     \code{blockF} (same as bF), \code{lvls_fac}, \code{nonnull}.}
#'   \item{bdat}{Block-level data.table with columns: \code{bF}, \code{blockF},
#'     \code{nonnull}, \code{lvls_fac}, \code{nb}, \code{hwt} (harmonic mean
#'     weight).}
#' }
#'
#' @examples
#' synth <- generate_synthetic_experiment(k = 4, max_level = 2, N_per_block = 50, t = 0.5)
#' dim(synth$idat)  # 800 x ...
#' dim(synth$bdat)  # 16 x ...
#'
#' @import data.table
#' @export
generate_synthetic_experiment <- function(k, max_level, N_per_block = 50L,
                                          t = 0.5, block_sd = 0.5) {
  stopifnot(k >= 2, max_level >= 1, N_per_block >= 4, N_per_block %% 2 == 0,
            t >= 0, t <= 1, block_sd >= 0)

  # 1. Generate tree structure with known non-null nodes
  treeDT <- generate_tree_DT(max_level = max_level, k = k, t = t)

  # 2. Build lvls_fac from node ancestry paths
  parent_lookup <- setNames(treeDT$parent, treeDT$node)

  get_path <- function(node) {
    path <- integer(0)
    current <- node
    while (!is.na(parent_lookup[as.character(current)])) {
      path <- c(path, current)
      current <- parent_lookup[as.character(current)]
    }
    rev(path)
  }

  treeDT[, lvls_fac := factor(vapply(node, function(n) {
    paste(get_path(n), collapse = ".")
  }, character(1)))]

  # 3. Extract leaf nodes as block-level data
  bdt <- treeDT[level == max_level, .(node, nonnull, lvls_fac)]
  bdt[, bF := factor(node)]
  bdt[, blockF := bF]
  bdt[, nb := N_per_block]

  N_total <- nrow(bdt) * N_per_block
  pb <- 0.5  # balanced randomization
  bdt[, hwt := (nb / N_total) * pb * (1 - pb)]

  setkey(bdt, bF)

  # 4. Create individual-level data
  bdt[, bary0 := rnorm(.N, mean = 0, sd = block_sd)]

  idat <- data.table(bF = factor(rep(bdt$node, each = N_per_block)))
  setkey(idat, bF)
  idat[, id := seq_len(.N)]

  # Merge block-level info
  idat <- merge(idat, bdt[, .(bF, bary0, nonnull, lvls_fac)], by = "bF")

  # Individual outcomes: y0 ~ N(block_mean, 1)
  idat[, y0 := rnorm(.N, mean = bary0, sd = 1)]

  # Block-randomized treatment: half treated, half control per block
  idat[, trt := sample(rep(c(0L, 1L), N_per_block / 2)), by = bF]
  idat[, trtF := factor(trt)]
  idat[, blockF := bF]

  # Clean up
  idat[, bary0 := NULL]
  bdt[, bary0 := NULL]
  bdt[, node := NULL]

  return(list(idat = idat, bdat = bdt))
}


#' Run All Testing Methods on a Synthetic Experiment
#'
#' Generates treatment effects on the synthetic data and runs all six testing
#' methods (four top-down, two bottom-up) via \code{\link{padj_test_fn}}.
#'
#' @param idat Individual-level data.table from
#'   \code{\link{generate_synthetic_experiment}}.
#' @param bdat Block-level data.table from
#'   \code{\link{generate_synthetic_experiment}}.
#' @param k Integer. Branching factor (must match the data).
#' @param max_level Integer. Tree depth (must match the data).
#' @param tau_size Numeric. Effect size parameter passed to \code{tau_fn}.
#' @param prop_blocks_0 Numeric between 0 and 1. Proportion of blocks with zero
#'   treatment effect. When the synthetic data has a known non-null pattern
#'   (from \code{generate_synthetic_experiment}), set to 0 and let the
#'   \code{nonnull} column control which blocks have effects.
#' @param tau_fn Function. Treatment effect generator (default
#'   \code{\link{tau_norm}}). See \code{\link{create_effects}}.
#' @param nsims Integer. Number of randomization reruns per method (default 200).
#' @param ncores_sim Integer. Number of cores for parallel reruns (default 1).
#' @param thealpha Numeric. Nominal significance level (default 0.05).
#' @param adaptive_tau Deprecated, ignored. The error load check in
#'   \code{compute_adaptive_alphas} replaces the old tau threshold.
#' @param methods Character vector. Which methods to run. Default: all six.
#'
#' @return A data.table with one row per method, containing columns
#'   \code{method} plus all error/power metrics from \code{\link{calc_errs_new}}.
#'
#' @import data.table
#' @importFrom manytestsr find_blocks splitSpecifiedFactorMulti alpha_adaptive
#' @export
run_synthetic_comparison <- function(idat, bdat, k, max_level,
                                     tau_size, prop_blocks_0 = 0,
                                     tau_fn = tau_norm, nsims = 200L,
                                     ncores_sim = 1L, thealpha = 0.05,
                                     adaptive_tau = 0.1,
                                     methods = c("td_unadj", "td_hommel",
                                                 "td_adaptive", "td_adaptive_hommel",
                                                 "bu_hommel", "bu_bh")) {

  N_total <- nrow(idat)

  # Build the adaptive alpha closure if needed
  adaptive_afn <- NULL
  if (any(grepl("td_adaptive", methods))) {
    adaptive_afn <- manytestsr::alpha_adaptive(
      k = k, delta_hat = max(tau_size, 0.01), N_total = N_total
    )
  }

  # Method -> parameter mapping
  method_specs <- list(
    td_unadj = list(
      p_adj_method = "split", splitfn = manytestsr::splitSpecifiedFactorMulti,
      afn = NULL, local_adj_p_fn = local_unadj_all_ps
    ),
    td_hommel = list(
      p_adj_method = "split", splitfn = manytestsr::splitSpecifiedFactorMulti,
      afn = NULL, local_adj_p_fn = local_hommel_all_ps
    ),
    td_adaptive = list(
      p_adj_method = "split", splitfn = manytestsr::splitSpecifiedFactorMulti,
      afn = adaptive_afn, local_adj_p_fn = local_unadj_all_ps
    ),
    td_adaptive_hommel = list(
      p_adj_method = "split", splitfn = manytestsr::splitSpecifiedFactorMulti,
      afn = adaptive_afn, local_adj_p_fn = local_hommel_all_ps
    ),
    bu_hommel = list(
      p_adj_method = "hommel", splitfn = NULL,
      afn = NULL, local_adj_p_fn = local_unadj_all_ps
    ),
    bu_bh = list(
      p_adj_method = "fdr", splitfn = NULL,
      afn = NULL, local_adj_p_fn = local_unadj_all_ps
    )
  )

  results_list <- list()

  for (m in methods) {
    spec <- method_specs[[m]]
    if (is.null(spec)) {
      warning("Unknown method: ", m, ". Skipping.")
      next
    }

    res <- padj_test_fn(
      idat = idat, bdat = bdat,
      blockid = "bF", trtid = "trt",
      fmla = Y ~ trtF | bF,
      ybase = "y0",
      prop_blocks_0 = prop_blocks_0,
      tau_fn = tau_fn, tau_size = tau_size,
      pfn = manytestsr::pOneway,
      afn = spec$afn,
      local_adj_p_fn = spec$local_adj_p_fn,
      p_adj_method = spec$p_adj_method,
      splitfn = spec$splitfn,
      splitby = if (!is.null(spec$splitfn)) "lvls_fac" else NULL,
      nsims = nsims,
      ncores = 1,
      ncores_sim = ncores_sim,
      thealpha = thealpha,
      stop_splitby_constant = TRUE,
      return_bottom_up = FALSE,
      blocksize = "hwt"
    )

    # padj_test_fn returns averaged error rates across nsims
    if (is.data.table(res)) {
      summary_row <- res[, lapply(.SD, mean, na.rm = TRUE)]
    } else {
      summary_row <- as.data.table(as.list(res))
    }
    summary_row[, method := m]
    results_list[[m]] <- summary_row
  }

  return(rbindlist(results_list, fill = TRUE))
}
