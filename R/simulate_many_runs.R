 #' Simulate the Family-Wise Error Rate (FWER) using a Data.table-based Tree with Alpha Options
 #'
 #' @description
 #' Runs the hierarchical testing procedure on complete k-ary trees (generated via \code{generate_tree_DT()})
 #' \code{n_sim} times and returns the estimated FWER. The user may choose whether to use a fixed \eqn{α},
 #' \eqn{α}–spending, or a rudimentary \eqn{α}–investing scheme.
 #'
 #' @param n_sim Integer. The number of simulation replicates.
 #' @param t Numeric in [0,1]. The probability that a leaf is non-null.
 #' @param k Integer. The branching factor.
 #' @param max_level Integer. The maximum level (depth) of the tree.
 #' @param alpha Numeric. The nominal significance level.
 #' @param N_total Numeric. The total sample size at the root.
 #' @param beta_base Numeric. The base parameter for the Beta distribution.
 #' @param adj_effN Logical. Whether to adjust the effective sample size at deeper levels.
 #' @param local_adj_p_fn Function. A function (e.g. \code{local_simes}) to adjust p-values at a node.
 #' @param return_details Logical. Whether to return the full simulated data.table.
 #' @param global_adj Character. The method to adjust the leaves (e.g. "hommel").
 #' @param alpha_method Character. One of \code{"fixed"}, \code{"spending"}, \code{"investing"}.
 #' @param final_global_adj Character. One of \code{"none"}, \code{"fdr"}, \code{"fwer"}.
 #' @param multicore Logical. Whether to use multiple cores (via \code{parallel::mclapply}).
 #' @param ... Additional arguments passed to \code{simulate_test_DT()}.
 #'
 #' @return A named numeric vector giving averaged error rates and discovery metrics across simulations.
 #'
 #' @importFrom parallel mclapply
 #' @export
 simulate_many_runs_DT <- function(n_sim, t, k, max_level, alpha, N_total, beta_base = 0.1,
                                   adj_effN = TRUE, local_adj_p_fn = local_simes, return_details = FALSE,
                                   global_adj = "hommel", alpha_method = "fixed", final_global_adj = "none", multicore = FALSE, ...) {
   treeDT <- generate_tree_DT(max_level, k, t)

   if (multicore) {
     ncores <- parallel::detectCores()
   } else {
     ncores <- 1
   }

   res <- mclapply(
     seq_len(n_sim),
     function(i) {
       simulate_test_DT(treeDT, alpha, k,
                        effN = N_total, N_total = N_total, beta_base = beta_base,
                        adj_effN = adj_effN, local_adj_p_fn = local_adj_p_fn, global_adj = global_adj,
                        alpha_method = alpha_method, return_details = return_details,
                        final_global_adj = final_global_adj, ...)
     },
     mc.cores = ncores, mc.set.seed = TRUE
   )

   res_dt <- data.table::rbindlist(res)
   mean_res <- unlist(res_dt[, lapply(.SD, mean, na.rm = TRUE)])
   return(mean_res)
 }