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
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param by_block Is an argument to [create_effects] to create true effects by block or across the whole dataset
#' @param pfn A function to produce pvalues --- using idat.
#' @param afn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
#' @param p_adj_method Is "split" to use [manytestsr::find_blocks] for top-down testing and "fdr" or "holm" etc to use [stats::p.adjust] and do the test in every block.
#' @param nsims Is the number of simulations to run --- each simulation uses the same treatment effects be re-assigns treatment (re-shuffles treatment and re-reveals the observed outcomes as a function of the potential outcomes)
#' @param ncores Tells p-value functions how many cores to use. Mostly ignored in use of this function because we are tending to parallelize at higher loops.
#' @param splitfn A function to split the data into two pieces --- using bdat
#' @param covariate is the name of a covariate to be used in created covariate dependent treatment effects. If NULL then the tau_fn should not use a covariate. If "newcov", then create a new covariate with a known (moderate) relationship with the potential outcome under control. This relationship is currently fixed with an R^2 of about .1.
#' @param splitby A string indicating which column in bdat contains a variable to guide splitting (for example, a column with block sizes or block harmonic mean weights or a column with a covariate (or a function of covariates))
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @param blocksize The character name of the column that measures the block size.
#' @param stop_splitby_constant TRUE is the algorithm should stop splitting when the splitting criteria is constant within set/parent or whether it should continue but split randomly.
#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"), and a
#' node level dataset  ("detnodes"). Default here is FALSE. Only use TRUE when
#' not using simulations.
#' @return A pvalue for each block
#' @export
padj_test_fn <- function(idat, bdat, blockid, trtid = "trt", fmla = Y ~ trtF | blockF, ybase,
                         prop_blocks_0, tau_fn, tau_size, by_block = TRUE, pfn, afn, p_adj_method, nsims, ncores = 1,
                         splitfn = NULL, covariate = NULL, splitby = NULL, thealpha = .05, blocksize = "hwt",
                         stop_splitby_constant = TRUE, return_details = FALSE) {
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
  # bdatnew$trueblocks <- bdat_effects$trueblocks

  if (p_adj_method == "split") {
    reveal_and_test_fn <- reveal_po_and_test_siup
  } else {
    reveal_and_test_fn <- reveal_po_and_test
  }

  p_sims_lst <- replicate(nsims, reveal_and_test_fn(
    idat = datnew, bdat =
      bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new", fmla = fmla,
    ybase = ybase, prop_blocks_0 = prop_blocks_0, tau_fn = tau_fn, tau_size =
      tau_size, pfn = pfn, afn = afn, p_adj_method = p_adj_method, splitfn =
      splitfn, splitby = splitby, thealpha = thealpha, stop_splitby_constant =
      stop_splitby_constant, ncores = ncores, return_details = return_details,
    blocksize = blocksize
  ), simplify = FALSE)

  if (length(p_sims_lst) == 1) {
    ## If we are using this function to apply to a given dataset, and so want
    ## details about the blocks, the results are a list of objects from
    ## calc_errs
    p_sims <- p_sims_lst[[1]]
  } else {
    p_sims <- data.table::rbindlist(p_sims_lst)
  }

  return(p_sims)
}

#' Repeat experiment, reveal treatment effects from the potential outcomes, test within each block, summarize
#'
#' The function does hypothesis tests within blocks and then summarizes the results  of this testing across the blocks. It is designed for use by padj_test_fn.
#' @param idat Individual level data
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with outcome~treatment assignment  where treatment assignment is a factor
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param y1var Is the name of the potential outcome to treatment
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param p_adj_method Is an input to p.adjust like "fdr" or "holm" to use adjust_block_tests
#' @param afn NULL. Not used. Included here to enable same functions to call this function and reveal_po_and_test_siup
#' @param splitfn Must be null. Only exists so that we can use the same efffect creation and testing function for both approaches.
#' @param splitby Must be null. Only exists so that we can use the same efffect creation and testing function for both approaches.
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @param copydts TRUE or FALSE. TRUE if using [manytestsr::find_blocks] standalone. FALSE if copied objects are being sent to find_blocks from other functions.
#' @param stop_splitby_constant FALSE (not used here because this algorithm does not split) is the algorithm should stop splitting when the splitting criteria is constant within set/parent or whether it should continue but split randomly.
#' @param ncores The number of cores or threads to use for the test statistic creation and possible permutation testing
#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"), and a
#' node level dataset  ("detnodes"). Default here is FALSE. Only use TRUE when
#' not using simulations.
#' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
#' @export
reveal_po_and_test <- function(idat, bdat, blockid, trtid, fmla = NULL, ybase, y1var,
                               prop_blocks_0, tau_fn, tau_size, pfn, p_adj_method = "fdr",
                               afn = NULL, splitfn = NULL, splitby = NULL, thealpha = .05, copydts = FALSE,
                               stop_splitby_constant = FALSE, ncores = 1, return_details = FALSE) {
  stopifnot(is.null(splitfn) | splitfn == "NULL")
  # make no effects within block by shuffling treatment, this is the engine of variability in the sim
  idat[, newZ := sample(get(trtid)), by = blockid]
  # reveal relevant potential outcomes with possible known effect
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)]
  # the pvalue functions want a factor
  idat[, newZF := factor(newZ)]
  fmla <- Y ~ newZF
  # idat[, Y := get(y1var) * get(trtid) + get(ybase) * (1 - get(trtid))] # reveal relevant potential outcomes with possible known effect

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

  errs <- calc_errs(
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
#' Repeat experiment, reveal treatment effects from the potential outcomes, test within partitions, summarize
#'
#' The function does hypothesis tests within partitions of blocks (including individual blocks depending on the splitting algorithm) and then summarizes the results  of this testing across the blocks. It very much depends  on padj_test_fn.
#' @param idat Data at the unit level.
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with outcome~treatment assignment  | block where treatment assignment and block must be factors. (NOT USED HERE)
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param y1var Is the name of the potential outcome to treatment
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param afn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
#' @param p_adj_method Must be "split" here.
#' @param copydts TRUE or FALSE. TRUE if using find_blocks standalone. FALSE if copied objects are being sent to find_blocks from other functions.
#' @param splitfn A function to split the data into two pieces --- using bdat
#' @param splitby A string indicating which column in bdat contains a variable to guide splitting (for example, a column with block sizes or block harmonic mean weights or a column with a covariate (or a function of covariates))
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @param stop_splitby_constant TRUE is the algorithm should stop splitting when the splitting criteria is constant within set/parent or whether it should continue but split randomly.
#' @param ncores The number of cores or threads to use for the test statistic creation and possible permutation testing
#' @param blocksize The character name of the column that measures the block size.
#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"), and a
#' node level dataset  ("detnodes"). Default here is FALSE. Only use TRUE when
#' not using simulations.
#' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
#' @export
reveal_po_and_test_siup <- function(idat, bdat, blockid, trtid, fmla = Y ~ newZF | blockF, ybase, y1var,
                                    prop_blocks_0, tau_fn, tau_size, pfn, afn, p_adj_method = "split",
                                    copydts = FALSE, splitfn, splitby, thealpha = .05, stop_splitby_constant = TRUE, ncores = 1,
                                    blocksize = "hwt", return_details = FALSE) {
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
  # Then reveal an observed outcome that contains a treatment effect from y1var as a function of newZ. So the treatment effect is known.
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)]
  # bdat in find_blocks doesn't really  carry important information other than block id
  # setkeyv(bdat, blockid)
  idat[, newZF := factor(newZ)]
  fmla <- as.formula(paste("Y~newZF|", blockid, sep = ""))
  # idat[, Y := get(y1var) * get(trtid) + get(ybase) * (1 - get(trtid))] # reveal relevant potential outcomes with possible known effect

  res <- find_blocks(
    idat = idat, bdat = bdat, blockid = blockid, splitfn =
      splitfn, pfn = pfn, alphafn = afn, thealpha = thealpha, fmla = fmla,
    parallel = parallel, ncores = ncores, copydts = copydts, splitby = splitby,
    blocksize = blocksize,
    stop_splitby_constant = stop_splitby_constant
  )

  errs <- calc_errs(
    testobj = res,
    truevar_name = "trueate",
    trueeffect_tol = .Machine$double.eps,
    blockid = blockid,
    return_details = return_details
  )

  return(errs)
}


#' Calculate the error and success proportions of tests for a single iteration
#'
#' @description

#' This function takes output from [manytestsr::find_blocks] or an equivalent
#' bottom-up testing function such as `adjust_block_tests` and returns the
#' proportions of errors made. To use this function the input to find_blocks
#' must include a column containing a true block-level effect. Repeated uses of
#' this function allow us to assess false discovery rates and family wise error
#' rates among other metrics of testing success.
#'
#' @param testobj Is an object arising from [manytestsr::find_blocks] or
#' [adjust_block_tests]. It will contain block-level results.
#' @param truevar_name Is a string indicating the name of the variable
#' containing the true underlying causal effect (at the block level).
#' @param trueeffect_tol Is the smallest effect size below which we consider
#' the effect to be zero (by default is it floating point zero).
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param thealpha Is the error rate for a given test (for cases where alphafn
#' is NULL, or the starting alpha for alphafn not null)
#' @param fwer Indicates that we are trying to control FWER. Right now, we do
#' this by default in report_detections but this indicates when we should be
#' looking at FDR
#' @param return_details TRUE means that the function should return a list of
#' the original data ("detobj"), a summary of the results ("detresults"),  a
#' node level dataset  ("detnodes"), and a copy of the original object that was
#' provided as input. Default here is FALSE. Only use TRUE when not using
#' simulations.

#'
#' @returns

#' If `hit` means that \eqn{p \le \alpha} for a given block or group of blocks
#' where testing has stopped, `allnull` means that all of the effects in the
#' group of blocks (or the single block) are zero, `anynotnull` means that at
#' least one block has a non-zero effect, then this is the code that returns
#' the different basic descriptions of errors.
#' \preformatted{
#' deterrs <- detnodes[, .(
#'  nreject = sum(hit),
#'  naccept = sum(1 - hit),
#'  prop_reject = mean(hit),
#'  prop_accept = mean(1 - hit), # or 1-prop_reject
#'  # rejections of false null /detections of true non-null
#'  # (if any of the blocks have an effect, then we count this as a correct rejection or correct detection)
#'  true_pos_prop = mean(hit * anynotnull),
#'  # If we reject but *all* of the blocks have no effects, this is a false positive error
#'  # If we reject but only one of them has no effects but the other has effects,
#'  # then this is not an error --- but a correct detection
#'  false_pos_prop = mean(hit * allnull),
#'  # If we do not reject and all blocks are truly null, then we have no error.
#'  true_neg_prop = mean((1 - hit) * allnull),
#'  # If we do not reject/detect and at least one of the blocks actually has an effect, we have
#'  # a false negative error --- a failure to detect the truth
#'  false_neg_prop = mean((1 - hit) * anynotnull),
#'  # Now look at false and true discoveries: false rejections as a proportion of rejections
#'  false_disc_prop = sum(hit * allnull) / max(1, sum(hit)),
#'  true_disc_prop = sum(hit * anynotnull) / max(1, sum(hit)),
#'  true_nondisc_prop = sum((1 - hit) * allnull) / max(1, sum(1 - hit)),
#'  false_nondisc_prop = sum((1 - hit) * anynotnull) / max(1, sum(1 - hit)),
#'  meangrpsize = mean(grpsize),
#'  medgrpsize = median(grpsize)
#' )]
#' }
#' We also summarize the average treatment effects in the blocks (or groups of blocks) among blocks where the null hypothesis of no effects has been rejected and where it has not been rejected.
#' \itemize{
#' \item "nreject" Number of tests with \eqn{p \le \alpha}
#'  \item "naccept" Number of tests with \eqn{p > \alpha}
#'  \item "prop_reject" Proportion of tests with \eqn{p \le \alpha} out of all tests completed
#'  \item "prop_accept" Proportion of tests with \eqn{p > \alpha} out of all tests completed
#'  \item "true_pos_prop" Proportion of tests with \eqn{p \le \alpha} out of all tests of true null hypotheses (rejections of false null hypotheses)
#'  \item "false_pos_prop" Proportion of tests with \eqn{p \le \alpha} out of all tests of false null hypotheses (rejections of true null hypotheses)
#'  \item "true_neg_prop"
#'  \item "false_neg_prop"
#'  \item "false_disc_prop"
#'  \item "true_disc_prop"
#'  \item "true_nondisc_prop"
#'  \item "false_nondisc_prop"
#'  \item "meangrpsize"
#'  \item "medgrpsize"
#'  \item "hit1"
#'  \item "minate1"
#'  \item "meanate1"
#'  \item "medate1"
#'  \item "maxate1"
#'  \item "hit0"
#'  \item "minate0"
#'  \item "meanate0"
#'  \item "medate0"
#'  \item "maxate0"
#' }
#' @export
calc_errs <- function(testobj,
                      truevar_name,
                      trueeffect_tol = .Machine$double.eps,
                      blockid = "bF",
                      thealpha = .05,
                      fwer = FALSE,
                      return_details = FALSE) {
  simp_summary <- function(x) {
    ## x is the true effect size in the block (probably the true mean effect size).
    x <- abs(x) ## We want look at relative sizes of detected effects and don't care about large negative versus positive effects
    list(min(x, na.rm = TRUE), mean(x, na.rm = TRUE), median(x, na.rm = TRUE), max(x, na.rm = TRUE))
  }

  if (length(grep("biggrp", names(testobj))) > 0) {
    # this is for the top-down/split and test method
    detobj <- report_detections(testobj, fwer = fwer, alpha = thealpha, only_hits = FALSE)
    ## This records both group hits and hits in individual blocks
    detobj[, hit := as.numeric(hit)]
    ## But from the perspective of error rates, we only care about rejections or non-rejections at the individual block level
    detobj[, hitb := as.numeric(max_p <= max_alpha & blocksbygroup == 1)]
    detobj[, hitb2 := as.numeric(single_hit)]
    stopifnot(all.equal(detobj$hitb, detobj$hitb2))
    ## Coding whether the true effect is zero or not by block.
    detobj[, true0 := as.numeric(abs(get(truevar_name)) <= trueeffect_tol)]
    detobj[, truenot0 := as.numeric(abs(get(truevar_name)) > trueeffect_tol)]

    ## Accessing the results another way
    # thetree <- make_results_tree(testobj, blockid = blockid) %>%
    #  select(-label) %>%
    #  as.data.frame()
    # sigleaves <- thetree %>% filter(out_degree == 0 & hit == 1)
    ##  detobj[blockF %in% sigleaves$bF,.(blockF,max_p,single_hit,hitb,hitb2,fin_grp)]
    ## stopifnot(all.equal(sort(sigleaves$bF), sort(as.character(detobj[hitb == 1, get(blockid)]))))

    ## detnodes_effects <- detobj[, simp_summary(get(truevar_name)), by = list(hit, hit_grp)]
    ## setnames(detnodes_effects, c("hit", "hit_grp", "minate", "meanate", "medianate", "maxate"))
    ## setkey(detnodes_effects, hit_grp)
    ## detnodes <- detnodes[detnodes_effects, ]
  } else {
    # This is for the bottom-up/test every block method
    thetree <- NA

    detobj <- testobj[, .(
      blockid = get(blockid),
      p,
      max_p,
      hit = max_p <= thealpha,
      hit_grp = get(blockid),
      truevar_name = get(truevar_name),
      hit = as.numeric(hit),
      hitb = as.numeric(hit),
      true0 = as.numeric(abs(get(truevar_name)) <= trueeffect_tol),
      truenot0 = as.numeric(abs(get(truevar_name)) > trueeffect_tol),
      trueeffect = get(truevar_name),
      grpsize = 1
    )]
    setnames(detobj, "truevar_name", truevar_name)
  }

  # This is calculated at the level of blocks not nodes.
  deterrs <- detobj[, .(
    nreject = sum(hitb), # number of rejected blocks
    naccept = sum(1 - hitb),
    prop_reject = mean(hitb),
    prop_accept = mean(1 - hitb),
    tot_true0 = sum(true0), # number of true 0 blocks.
    tot_truenot0 = sum(truenot0),
    tot_reject_true0 = sum(hitb * true0),
    tot_not_reject_true0 = sum((1 - hitb) * true0),
    tot_reject_truenot0 = sum(hitb * truenot0),
    tot_not_reject_truenot0 = sum((1 - hitb) * truenot0),
    # or 1-prop_reject
    # Proportion of the total blocks that have an effect where we detect that effect:
    correct_pos_effect_prop = sum(hitb * truenot0) / max(1, sum(truenot0)),
    correct_pos_nulls_prop = sum(hitb * true0) / max(1, sum(true0)),
    ## if hitb=1 (reject) and truenot0=1 then this is a correct rejection. Use in power calculations.
    true_pos_prop = mean(hitb * truenot0),
    ## Proportion of rejections of a true null out of the total number of tests
    ##   if hitb=1 (meaning rejection) and true0=1 (meaning no true effect, true effect = 0) then a false positive
    false_pos_prop = mean(hitb * true0), # number blocks rejected when truly 0 / total number of blocks.
    # also sum(hitb*true0)/nrow(detobj)
    # If we do not reject and all blocks are truly null, then we have no error.
    prop_not_reject_true0 = mean((1 - hitb) * true0),
    # If we do not reject/detect and at least one of the blocks actually has an effect, we have
    # a false negative error --- a failure to detect the truth: proportion accept when true is not 0
    false_neg_prop = mean((1 - hitb) * truenot0),
    # Now look at false and true discoveries: false rejections as a proportion of rejections
    false_disc_prop = sum(hitb * true0) / max(1, sum(hitb)),
    true_disc_prop = sum(hitb * truenot0) / max(1, sum(hitb)),
    # Failure to reject when the null is not 0 out of the total failures to reject
    false_nondisc_prop = sum((1 - hitb) * truenot0) / max(1, sum(1 - hitb)),
    ## Failure to reject when the null is true (we want this) out of total accepts/failures to reject
    true_nondisc_prop = sum((1 - hitb) * true0) / max(1, sum(1 - hitb))
    # meangrpsize = mean(grpsize),
    # medgrpsize = median(grpsize)
  )]
  # One row of results
  #  detresults <- cbind(deterrs, detates)
  if (!return_details) {
    return(deterrs)
  } else {
    res <- list(
      detresults = deterrs,
      detobj = detobj,
      # detnodes = detnodes,
      testobj = testobj,
      tree = thetree
    )
    return(res)
  }
}
