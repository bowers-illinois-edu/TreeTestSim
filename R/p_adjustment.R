#' @title Compute local Simes p-value for a Vector of Child p-values
#'
#' @description Given \eqn{k} child p-values, computes the Simes p-value
#'   \eqn{\min_{i=1\ldots k} \{ (k/i) * p_{(i)} \}}, where \eqn{p_{(1)} \le \ldots \le p_{(k)}}.
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used in this function)
#'
#' @return A single numeric value: the Simes combination p-value.
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
#' @description Given \eqn{k} child p-values, computes the Hommel-adjusted p-values
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used in this function)
#'
#' @return A vector of adjusted p-values
#' @importFrom hommel hommel
#' @examples
#' local_hommel_all_ps(c(0.01, 0.04, 0.10, 0.20))
#' @export
local_hommel_all_ps <- function(pvals_children, alpha = .05) {
  adj_p_vals <- hommel(pvals_children)@adjusted
  return(adj_p_vals)
}

#' @title Unadjusted local minimal p-value
#'
#' @description Given \eqn{k} child p-values, returns the highest p-value below alpha;
#'   if none are below alpha returns the smallest p-value.
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha
#'
#' @return A single numeric value.
#' @examples
#' local_min_p(c(0.01, 0.04, 0.10, 0.20))
#' local_min_p(c(0.10, 0.20))
#' @export
local_min_p <- function(pvals_children, alpha = .05) {
  p_le_alpha <- pvals_children[pvals_children <= alpha]
  if (length(p_le_alpha) == 0) {
    return(min(pvals_children))
  } else {
    return(max(p_le_alpha))
  }
}

#' @title Unadjusted local step (pass-through)
#'
#' @description Returns the input p-values unmodified.
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used)
#'
#' @return A numeric vector: the same as input.
#' @examples
#' local_unadj_all_ps(c(0.01, 0.04))
#' @export
local_unadj_all_ps <- function(pvals_children, alpha = .05) {
  pvals_children
}

#' @title BH local adjustment
#'
#' @description Performs Benjamini-Hochberg adjustment on a vector of p-values.
#'
#' @param pvals_children Numeric vector of child p-values.
#' @param alpha Numeric scalar of alpha (not used)
#'
#' @return A numeric vector of BH-adjusted p-values.
#' @examples
#' local_bh_all_ps(c(0.01, 0.04, 0.10, 0.20))
#' @export
local_bh_all_ps <- function(pvals_children, alpha = .05) {
  p.adjust(pvals_children, method = "BH")
}
