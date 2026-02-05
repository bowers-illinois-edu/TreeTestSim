# Package-level declarations to satisfy R CMD check

#' @importFrom stats as.formula median p.adjust pbeta qnorm rbeta rnorm runif sd setNames
#' @keywords internal
"_PACKAGE"

# Declare data.table column names used in non-standard evaluation
# This prevents "no visible binding for global variable" notes
utils::globalVariables(c(

  # Tree structure columns
  "node", "level", "parent", "nonnull",
  # Simulation columns
  "p_val", "p_sim", "alpha_alloc", "local_adj_p", "p_val_final_adj",
  "bottom_up_p_adj", "bu_hommel_p", "bu_bh_p",
  # Block/individual data columns
  "bF", "blockF", "lvls_fac", "nb", "hwt", "bary0",
  "id", "y0", "trt", "trtF", "Y",
  # Effect creation columns
  "y1sim", "y1new", "trueblocks", "trueate", "truetaui",
  # Testing columns
  "p", "max_p", "max_alpha", "hit", "hitb", "hitb2", "single_hit",
  "blocksbygroup", "true0", "truenot0", "a",
  # Simulation helper columns
  "newZ", "newZF", "newcov", "nodenum_current",
  # Summary columns
  "method", "deterrs",
  # data.table's . function
  "."
))
