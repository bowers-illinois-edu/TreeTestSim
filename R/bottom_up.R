#' Test every block and adjust the p-values
#'
#' This function tests a hypothesis (of no effects currently) within each and
#' every block and then adjusts the results to control the False Discovery Rate
#' or Family-wise Error Rate  using the base R \code{\link{p.adjust}} function.
#'
#' @param idat Data at the unit level as a data.table.
#' @param bdat Data at the block level as a data.table.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param p_adj_method  A string indicating which method `p.adjust` should use. Default is "fdr".
#' @param simthresh Below which number of total observations should the p-value functions use permutations rather than asymptotic approximations
#' @param sims Number of permutations for permutation-based testing when the total number of cases is below the threshold.
#' @param fmla A formula with `outcome~treatment assignment | block` where treatment assignment and block must be factors. The tests will run within each block so that `|block` part is not strictly important.
#' @param parallel Should the pfn use multicore processing for permutation based testing? Default is no. But could be "snow" or "multicore" following `approximate` in the coin package.
#' @param ncores The number of cores to use in the p-value creation function
#' @param copydts TRUE or FALSE. TRUE if using standalone. Could be FALSE if copied objects are being sent to this function from other functions. Basically the question is whether we want this function to change the block or individual level datasets in the main environment or return a new object.
#' @return A data.table containing information about the sequence tests. The column `max_p` is the adjusted p-value and the column `p` is the unadjusted p-value.
#' @importFrom stringi stri_count_fixed stri_split_fixed stri_split stri_sub
#' @import data.table
#' @export
adjust_block_tests <- function(idat, bdat, blockid = "block", pfn, p_adj_method = "fdr", simthresh = 20,
                               sims = 1000, fmla = YContNorm ~ trtF,
                               parallel = "multicore", ncores = 4, copydts = FALSE) {
  if (copydts) {
    bdat <- copy(bdat)
    idat <- copy(idat)
  }
  theps <- idat[, .(p = pfn(fmla = fmla, dat = .SD, parallel = parallel, ncpu = ncores, simthresh = simthresh, sims = sims)), by = blockid]
  setkeyv(theps, blockid)
  setkeyv(bdat, blockid)
  res <- merge(theps, bdat)
  res[, max_p := p.adjust(p, method = p_adj_method)]
  res[, nodenum_current := get(blockid)]
  return(res)
}
