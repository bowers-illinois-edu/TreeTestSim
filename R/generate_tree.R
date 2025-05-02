#' Generate a Complete k-ary Tree as a data.table
#'
#' @description
#' Constructs a complete k-ary tree with levels 0 to \code{max_level}.
#' The tree is stored in a data.table with one row per node and columns for:
#' \code{node} (the node's index),
#' \code{level} (its level in the tree, with 0 = root),
#' \code{parent} (the index of its parent; NA for the root), and
#' \code{nonnull} (a logical indicating whether the node is non-null).
#'
#' At the leaves (nodes at \code{level == max_level}), the \code{nonnull} flag
#' is sampled as \code{TRUE} with probability \code{t}. For internal nodes,
#' \code{nonnull} is set to \code{TRUE} if any of its children are non-null.
#'
#' @param max_level Integer. The maximum level of the tree (with the root at level 0).
#' @param k Integer. The branching factor (each internal node has exactly \code{k} children).
#' @param t Numeric in \eqn{[0,1]}. The probability that a leaf is non-null (has a treatment effect).
#'
#' @return A data.table with columns \code{node}, \code{level}, \code{parent}, and \code{nonnull}.
#' @export
generate_tree_DT <- function(max_level, k, t) {
  total_nodes <- sum(k^(0:max_level))
  treeDT <- data.table(node = 1:total_nodes)
  levels <- integer(total_nodes)
  cum_nodes <- cumsum(k^(0:max_level))
  start <- 1
  for (l in 0:max_level) {
    end <- cum_nodes[l + 1]
    levels[start:end] <- l
    start <- end + 1
  }
  treeDT[, level := levels]
  treeDT[, parent := ifelse(node == 1, NA_integer_, floor((node - 2) / k) + 1)]
  treeDT[, nonnull := NA]
  treeDT[level == max_level, nonnull := {
    n_leaves <- .N
    m_leaves <- floor(t * n_leaves)
    sample(rep(c(FALSE, TRUE), c(n_leaves - m_leaves, m_leaves)))
  }]
  if (max_level > 0) {
    for (l in (max_level - 1):0) {
      treeDT[level == l, nonnull := {
        sapply(node, function(i) {
          children <- ((i - 1) * k + 2):((i - 1) * k + k + 1)
          children <- children[children <= total_nodes]
          any(treeDT[children, nonnull], na.rm = TRUE)
        })
      }]
    }
  } else {
    treeDT[node == 1, nonnull := runif(1) < t]
  }
  return(treeDT)
}
