# Functions to assess the procedures via simulation. Mostly for use in the paper itself.

#' Create treatment effect, then add tau and test using the SIUP method
#'
#' This function returns a function that carries with it the environment containing the arguments.
#' The idea is to make parallelization easier.
#'
#' @param idat Individual level data
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with \code{outcome~treatment assignment  | block}
#' where treatment assignment and block must be factors.
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all

#' @param tau_fn Is a function that turns ybase into the potential outcome
#' under treatment --- it is a treatment effect creating function.

#' @param tau_size Is the parameter for the tau_fn --- like the true average
#' effect size within a block.

#' @param by_block Is an argument to [create_effects] to create true effects by
#' block or across the whole dataset

#' @param pfn A function to produce pvalues --- using idat.

#' @param afn A function to adjust alpha at each step. Takes one or more
#' p-values plus a stratum or batch indicator.

#' @param local_adj_p_fn A function used to adjust p-values at each step (the
#' p-values of the children of a given parent.)

#' @param p_adj_method Is "split" to use [manytestsr::find_blocks] for top-down
#' testing and "fdr" or "holm" etc to use [stats::p.adjust] and do the test in
#' every block.

#' @param nsims Is the number of simulations to run --- each simulation uses
#' the same treatment effects be re-assigns treatment (re-shuffles treatment
#' and re-reveals the observed outcomes as a function of the potential
#' outcomes)

#' @param ncores Tells p-value functions how many cores to use. Mostly ignored
#' in use of this function because we are tending to parallelize at higher
#' loops.

#' @param ncores_sim Number of cores used to repeat the simulation.
#' \code{ncores} should be 1 if this is more than one.

#' @param splitfn A function to split the data into two pieces --- using bdat

#' @param covariate is the name of a covariate to be used in created covariate
#' dependent treatment effects. If NULL then the tau_fn should not use a
#' covariate. If "newcov", then create a new covariate with a known (moderate)
#' relationship with the potential outcome under control. This relationship is
#' currently fixed with an R^2 of about .1.

#' @param splitby A string indicating which column in bdat contains a variable
#' to guide splitting (for example, a column with block sizes or block harmonic
#' mean weights or a column with a covariate (or a function of covariates))

#' @param thealpha Is the error rate for a given test (for cases where alphafn
#' is NULL, or the starting alpha for alphafn not null)

#' @param blocksize The character name of the column that measures the block
#' size.

#' @param stop_splitby_constant TRUE is the algorithm should stop splitting
#' when the splitting criteria is constant within set/parent or whether it
#' should continue but split randomly.

#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"), and a
#' node level dataset  ("detnodes"). Default here is FALSE. Only use TRUE when
#' not using simulations.

#' @param bottom_up_adj Is a string for the `p.adjust` function such as
#' "hommel" or "fdr". When doing simulations where `p_adj_method` include such
#' approaches, care is needed.(TODO fix this)

#' @return A pvalue for each block
#' @import manytestsr
#' @export
padj_test_fn <- function(idat, bdat, blockid, trtid = "trt", fmla = Y ~ trtF | blockF, ybase,
                         prop_blocks_0, tau_fn, tau_size, by_block = TRUE,
                         pfn, afn, local_adj_p_fn = local_unadj_all_ps, p_adj_method, nsims, ncores = 1,
                         ncores_sim = 1, bottom_up_adj = "hommel",
                         splitfn = NULL, covariate = NULL, splitby = NULL, thealpha = .05, blocksize = "hwt",
                         stop_splitby_constant = TRUE, return_details = FALSE, return_bottom_up = TRUE) {
  if (!is.null(afn) & is.character(afn)) {
    if (afn == "NULL") {
      afn <- NULL
    } else {
      afn <- get(afn)
    }
  }

  datnew <- copy(idat)
  bdatnew <- copy(bdat)

  setkeyv(datnew, blockid)
  setkeyv(bdatnew, blockid)

  ## make data the same way for each design. This is a hack to have the seed within the function.
  set.seed(12345)

  ## If covariate="newcov" then make a covariate with a known relationship with
  ## the potential outcome to control (this to avoid problems with some created
  ## covariates perfectly predicting the outcome under control)
  if ((!is.null(covariate) && covariate == "newcov") || (!is.null(splitby) && splitby == "newcov")) {
    datnew[, newcov := {
      tmp <- .01 * sd(get(ybase)) * get(ybase) + rnorm(.N, mean = 0, sd = sd(get(ybase)))
      ifelse(tmp < 0, 0, tmp)
    }, by = blockid]
    newcovb <- datnew[, .(newcov = mean(newcov)), by = blockid]
    bdatnew <- bdatnew[newcovb]
  }

  ## Create potential outcome to treatment using a tau_fn and tau_size etc.
  datnew$y1new <- create_effects(
    idat = datnew, ybase = ybase, blockid = blockid,
    tau_fn = tau_fn, tau_size = tau_size,
    prop_blocks_0 = prop_blocks_0, covariate = covariate, by_block = by_block
  )
  datnew[, trueblocks := ifelse(abs(y1new - get(ybase)) <= .Machine$double.eps, 0, 1)]
  datnew[, trueate := mean(y1new - get(ybase)), by = blockid]
  bdat_effects <- datnew[, .(
    trueblocks = unique(trueblocks),
    trueate = mean(trueate)
  ), by = blockid]
  stopifnot(all(bdat_effects$trueblocks %in% c(0, 1)))
  setkeyv(bdat_effects, blockid)
  bdatnew <- bdatnew[bdat_effects, ]
  ## TODO: probably dont need both nonnull and trueblocks eventually
  bdatnew[, nonnull := (trueblocks == 1)]

  if (p_adj_method == "split") {
    reveal_and_test_fn <- reveal_po_and_test_siup
  } else {
    reveal_and_test_fn <- reveal_po_and_test
  }

  if (ncores_sim > 1) {
    p_sims_lst <- mclapply(1:nsims, function(i) {
      reveal_and_test_fn(
        idat = datnew, bdat = bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new", fmla = fmla,
        ybase = ybase, prop_blocks_0 = prop_blocks_0, tau_fn = tau_fn, tau_size = tau_size, pfn = pfn, afn = afn,
        local_adj_p_fn = local_adj_p_fn, p_adj_method = p_adj_method,
        splitfn = splitfn, splitby = splitby, thealpha = thealpha,
        stop_splitby_constant = stop_splitby_constant, ncores = ncores, return_details = return_details,
        return_bottom_up = return_bottom_up, bottom_up_adj = bottom_up_adj,
        blocksize = blocksize
      )
    }, mc.cores = ncores_sim)
  } else {
    p_sims_lst <- replicate(nsims, reveal_and_test_fn(
      idat = datnew, bdat =
        bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new", fmla = fmla,
      ybase = ybase, prop_blocks_0 = prop_blocks_0, tau_fn = tau_fn, tau_size =
        tau_size, pfn = pfn, afn = afn, local_adj_p_fn = local_adj_p_fn, p_adj_method = p_adj_method, splitfn =
        splitfn, splitby = splitby, thealpha = thealpha, stop_splitby_constant =
        stop_splitby_constant, ncores = ncores, return_details = return_details, return_bottom_up = return_bottom_up,
      blocksize = blocksize, bottom_up_adj = bottom_up_adj
    ), simplify = FALSE)
  }

  if (length(p_sims_lst) == 1) {
    ## If we are using this function to apply to a given dataset, and so want
    ## details about the blocks, the results are a list of objects which summarize errors and discoveries
    p_sims <- p_sims_lst[[1]]
  } else {
    p_sims <- data.table::rbindlist(p_sims_lst)
  }

  return(p_sims)
}

#' Repeat experiment, reveal treatment effects from the potential outcomes, test within each block, summarize

#'
#' The function does hypothesis tests within blocks and then summarizes the
#' results  of this testing across the blocks. It is designed for use by
#' padj_test_fn.

#' @param idat Individual level data
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with outcome~treatment assignment  where treatment assignment is a factor
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param y1var Is the name of the potential outcome to treatment
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all

#' @param tau_fn Is a function that turns ybase into the potential outcome
#' under treatment --- it is a treatment effect creating function.

#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param p_adj_method Is an input to p.adjust like "fdr" or "holm" to use adjust_block_tests
#' @param afn NULL. Not used. Included here to enable same functions to call this function and reveal_po_and_test_siup

#' @param splitfn Must be null. Only exists so that we can use the same efffect
#' creation and testing function for both approaches.

#' @param splitby Must be null. Only exists so that we can use the same efffect
#' creation and testing function for both approaches.

#' @param thealpha Is the error rate for a given test (for cases where alphafn
#' is NULL, or the starting alpha for alphafn not null)

#' @param copydts TRUE or FALSE. TRUE if using [manytestsr::find_blocks]
#' standalone. FALSE if copied objects are being sent to find_blocks from other
#' functions.

#' @param stop_splitby_constant FALSE (not used here because this algorithm
#' does not split) is the algorithm should stop splitting when the splitting
#' criteria is constant within set/parent or whether it should continue but
#' split randomly.

#' @param ncores The number of cores or threads to use for the test statistic creation and possible permutation testing
#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"), and a
#' node level dataset  ("detnodes"). Default here is FALSE. Only use TRUE when
#' not using simulations.

#' @return False positive proportion out of the tests across the blocks, The
#' false discovery rate (proportion rejected of false nulls out of all
#' rejections), the power of the adjusted tests across blocks (the proportion
#' of correctly rejected hypotheses out of all correct hypotheses --- in this
#' case correct means non-null), and power of the unadjusted test (proportion
#' correctly rejected out of  all correct hypothesis, but using unadjusted
#' p-values).

#' @export
reveal_po_and_test <- function(idat, bdat, blockid, trtid, fmla = NULL, ybase, y1var,
                               prop_blocks_0, tau_fn, tau_size, pfn, p_adj_method = "fdr",
                               afn = NULL, splitfn = NULL, splitby = NULL, thealpha = .05, copydts = FALSE,
                               stop_splitby_constant = FALSE, ncores = 1, return_details = FALSE, local_adj_p_fn, return_bottom_up, blocksize, bottom_up_adj) {
  stopifnot(is.null(splitfn) | splitfn == "NULL")
  # make no effects within block by shuffling treatment, this is the engine of variability in the sim
  idat[, newZ := sample(get(trtid)), by = blockid]
  # reveal relevant potential outcomes with possible known effect made manifest in y1var
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)]
  # the pvalue functions want a factor
  idat[, newZF := factor(newZ)]
  fmla <- Y ~ newZF

  if (ncores > 1) {
    parallel <- "multicore"
  }
  if (ncores == 1) {
    parallel <- "no"
  }

  res <- adjust_block_tests(
    idat = idat, bdat = bdat, blockid = blockid, p_adj_method = p_adj_method,
    pfn = pfn, fmla = fmla,
    copydts = copydts, ncores = ncores, parallel = parallel
  )

  res[, hit := max_p < thealpha]

  errs <- calc_errs_new(
    testobj = res,
    truevar_name = "trueate",
    trueeffect_tol = .Machine$double.eps,
    blockid = blockid,
    return_details = return_details
  )

  return(errs)
}

#' Repeat experiment, reveal treatment effects from the potential outcomes, test within partitions, summarize
#'

#' The function does hypothesis tests within partitions of blocks (including
#' individual blocks depending on the splitting algorithm) and then summarizes
#' the results  of this testing across the blocks. It very much depends  on
#' padj_test_fn.
#'

#' @param idat Data at the unit level.
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable

#' @param fmla A formula with outcome~treatment assignment  | block where
#' treatment assignment and block must be factors. (NOT USED HERE)

#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param y1var Is the name of the potential outcome to treatment
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all

#' @param tau_fn Is a function that turns ybase into the potential outcome
#' under treatment --- it is a treatment effect creating function.

#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param afn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
#' @param p_adj_method Must be "split" here.

#' @param copydts TRUE or FALSE. TRUE if using find_blocks standalone. FALSE if
#' copied objects are being sent to find_blocks from other functions.

#' @param splitfn A function to split the data into two pieces --- using bdat

#' @param splitby A string indicating which column in bdat contains a variable
#' to guide splitting (for example, a column with block sizes or block harmonic
#' mean weights or a column with a covariate (or a function of covariates))

#' @param thealpha Is the error rate for a given test (for cases where alphafn
#' is NULL, or the starting alpha for alphafn not null)

#' @param stop_splitby_constant TRUE is the algorithm should stop splitting
#' when the splitting criteria is constant within set/parent or whether it
#' should continue but split randomly.

#' @param ncores The number of cores or threads to use for the test statistic creation and possible permutation testing
#' @param blocksize The character name of the column that measures the block size.
#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"), and a
#' node level dataset  ("detnodes"). Default here is FALSE. Only use TRUE when
#' not using simulations.

#' @return False positive proportion out of the tests across the blocks, The
#' false discovery rate (proportion rejected of false nulls out of all
#' rejections), the power of the adjusted tests across blocks (the proportion
#' of correctly rejected hypotheses out of all correct hypotheses --- in this
#' case correct means non-null), and power of the unadjusted test (proportion
#' correctly rejected out of  all correct hypothesis, but using unadjusted
#' p-values).

#' @export
reveal_po_and_test_siup <- function(idat, bdat, blockid, trtid, fmla = Y ~ newZF | blockF, ybase, y1var,
                                    prop_blocks_0, tau_fn, tau_size, pfn, afn, local_adj_p_fn, p_adj_method = "split",
                                    copydts = FALSE, splitfn, splitby, thealpha = .05,
                                    stop_splitby_constant = TRUE, ncores = 1, truevar_name = "nonnull",
                                    blocksize = "hwt", return_details = FALSE, return_bottom_up = TRUE, bottom_up_adj) {
  if (!is.null(afn) & is.character(afn)) {
    if (afn == "NULL") {
      afn <- NULL
    } else {
      afn <- get(afn)
    }
  }

  if (ncores > 1) {
    parallel <- "multicore"
  } else {
    parallel <- "no"
  }
  # 	message(tau_size)
  # 	make relationship between treatment and outcome 0 on average in preparation for power analysis and error rate sims
  idat[, newZ := sample(get(trtid)), by = blockid]

  # Then reveal an observed outcome that contains a treatment effect from y1var
  # as a function of newZ. So the treatment effect is known.

  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)]

  # bdat in find_blocks doesn't really  carry important information other than block id
  # setkeyv(bdat, blockid)
  idat[, newZF := factor(newZ)]
  fmla <- as.formula(paste("Y~newZF|", blockid, sep = ""))

  res <- find_blocks(
    idat = idat, bdat = bdat, blockid = blockid,
    splitfn = splitfn, pfn = pfn, alphafn = afn,
    local_adj_p_fn = local_adj_p_fn,
    thealpha = thealpha, fmla = fmla,
    parallel = parallel, ncores = ncores, copydts = copydts, splitby = splitby,
    blocksize = blocksize,
    stop_splitby_constant = stop_splitby_constant, return_what = c("blocks", "nodes")
  )

  if (return_bottom_up) {
    res_bottom_up <- adjust_block_tests(
      idat = idat, bdat = bdat, blockid = blockid,
      p_adj_method = bottom_up_adj, pfn = pfn, fmla = fmla,
      copydts = copydts, ncores = ncores, parallel = parallel
    )

    res_bottom_up[, hit := max_p < thealpha]

    errs_bottom_up <- calc_errs_new(
      testobj = res_bottom_up,
      truevar_name = "trueate",
      trueeffect_tol = .Machine$double.eps,
      blockid = blockid,
      return_details = return_details
    )
  }

  res_tree <- make_results_tree(orig_res = res$bdat, block_id = blockid, truevar_name = truevar_name)

  ## errs_tree <- calc_errs_new(
  ##  testobj = res$bdat,
  ##  truevar_name = "trueate",
  ##  trueeffect_tol = .Machine$double.eps,
  ##  blockid = blockid,
  ##  return_details = return_details
  ## )
  errs_tree <- res_tree$test_summary

  if (return_details && !return_bottom_up) {
    results <- list(tree = res_tree, res = res)
    return(results)
  }
  if (!return_details && !return_bottom_up) {
    return(errs_tree)
  }
  if (!return_details && return_bottom_up) {
    names(errs_bottom_up) <- paste0("bot_", names(errs_bottom_up))
    results <- cbind(errs_tree = errs_tree, errs_bottom_up = errs_bottom_up)
    return(results)
  }
  if (return_details && return_bottom_up) {
    results <- list(
      tree = res_tree, errs_tree = errs_tree,
      bottom_up = res_bottom_up, errs_bottom_up = errs_bottom_up
    )
    return(results)
  }
}
