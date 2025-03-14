% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeSimFWER.R
\name{simulate_test}
\alias{simulate_test}
\title{Simulate the Testing Procedure on a Tree Node}
\usage{
simulate_test(node, parent_p, alpha, k, effN, N_total, beta_base)
}
\arguments{
\item{node}{A list representing a node (or subtree) created by \code{generate_tree()}.}

\item{parent_p}{A numeric value giving the lower bound for the p-value inherited from the parent.}

\item{alpha}{A numeric value for the significance level.}

\item{k}{An integer for the number of children per node.}

\item{effN}{A numeric value for the effective sample size available at the current node.}

\item{N_total}{A numeric value for the total sample size available at the root.}

\item{beta_base}{A numeric value for the base parameter of the beta distribution (used at the root).}
}
\value{
A list with components:
\describe{
\item{p_val}{The simulated p-value at the current node.}
\item{false_error}{A logical value: \code{TRUE} if a type I error occurred somewhere in the branch, otherwise \code{FALSE}.}
}
}
\description{
Recursively simulate the testing procedure along a tree branch. At each node, a p-value is
generated using a uniform distribution if the node is null or a rescaled beta distribution if non-null.
The function applies a global gate by computing a Simes test over the children.
}
\examples{
tree <- generate_tree(0, max_level = 3, k = 3, t = 0.2)
res <- simulate_test(tree, parent_p = 0, alpha = 0.05, k = 3, effN = 1000, N_total = 1000, beta_base = 0.1)
print(res)

}
