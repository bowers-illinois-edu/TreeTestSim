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
#' We also summarize the average treatment effects in the blocks (or groups of
#' blocks) among blocks where the null hypothesis of no effects has been
#' rejected and where it has not been rejected.
#' @importFrom manytestsr make_results_tree report_detections
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
    ## We want look at relative sizes of detected effects and don't care about large negative versus positive effects
    x <- abs(x)
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

#' Calculate the error and success proportions of tests for a single iteration
#'
#' @description

#' This function takes output from [manytestsr::find_blocks] or an equivalent
#' bottom-up testing function such as [adjust_block_tests] and returns the
#' proportions of errors and discoveries made across the multiple tests. To use
#' this function the input to [manytestsr::find_blocks] must include a column
#' containing a true block-level effect in `truevar_name`. Repeated uses of
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
#' the original data, a summary of the results,  a node level dataset, and a
#' copy of the original object that was provided as input. Default here is
#' FALSE. Only use TRUE when not using simulations.

#'
#' @returns

#' Most basically this function returns a single vector (or data.table with 1
#' row) with a variety of summaries of the performance of the testing procedure
#' across the blocks (i.e. leaves) and nodes (i.e. groups of blocks) in the
#' design.

#' @importFrom manytestsr make_results_tree report_detections
#' @export
calc_errs_new <- function(testobj,
                          truevar_name,
                          trueeffect_tol = .Machine$double.eps,
                          blockid = "bF",
                          thealpha = .05,
                          fwer = FALSE,
                          return_details = FALSE) {
  simp_summary <- function(x) {
    ## x is the true effect size in the block (probably the true mean effect size).
    ## We want look at relative sizes of detected effects and don't care about large negative versus positive effects
    x <- abs(x)
    list(min(x, na.rm = TRUE), mean(x, na.rm = TRUE), median(x, na.rm = TRUE), max(x, na.rm = TRUE))
  }

  if (length(grep("biggrp", names(testobj))) > 0) {
    ## Do both node and leaf level calculations here
    res <- make_results_tree(testobj,
      block_id = blockid,
      truevar_name = truevar_name, return_what = c("nodes", "test_summary")
    )
    test_summary <- res$test_summary
  } else {
    ## This is for the bottom-up approach where testobj has a p-value for each block or leaf
    ## there are no node level calculations since all nodes are leaves
    res <- testobj[, .(
      blockid = get(blockid),
      a = thealpha,
      p,
      max_p,
      hit = max_p <= thealpha,
      hit_grp = get(blockid),
      truevar_name = get(truevar_name)
    )]

    res[, hit := as.numeric(hit)]
    res[, hitb := as.numeric(hit)]
    res[, nonnull := abs(as.numeric(truevar_name)) > trueeffect_tol]
    setnames(res, "truevar_name", truevar_name)
    # res[, nonnull := truenot0 == 1]

    num_nodes <- NA
    # the block level dataset has one row for each leaf or block
    num_leaves <- nrow(testobj)
    num_nodes_tested <- NA
    num_nonnull_nodes_tested <- NA
    node_rejections <- NA
    node_any_false_rejection <- NA
    node_false_rejection_prop <- NA
    node_num_false_rejections <- NA
    node_false_discovery_prop <- NA
    node_true_discoveries <- NA
    node_power <- NA

    num_leaves_tested <- sum(res[, !is.na(max_p)])
    num_nonnull_leaves_tested <- sum(!is.na(res$max_p) & res$nonnull)

    if (num_leaves_tested > 0) {
      leaf_power <- res[nonnull == TRUE & !is.na(max_p), mean(max_p <= a, na.rm = TRUE)]
      leaf_rejections <- res[!is.na(max_p), sum(max_p <= a, na.rm = TRUE)]
      leaf_true_discoveries <- res[nonnull == TRUE & !is.na(max_p), sum(max_p <= a, na.rm = TRUE)]
      leaf_any_false_rejection <- res[nonnull == FALSE & !is.na(max_p), any(max_p <= a)]
      leaf_false_rejection_prop <- res[nonnull == FALSE & !is.na(max_p), mean(max_p <= a)]
      leaf_false_discovery_prop <- res[nonnull == FALSE & !is.na(max_p), sum(max_p <= a) / max(1, leaf_rejections)]
    } else {
      leaf_power <- NA
      leaf_rejections <- NA
      leaf_true_discoveries <- NA
      leaf_any_false_rejection <- NA
      leaf_false_rejection_prop <- NA
      leaf_false_discovery_prop <- NA
    }

    test_summary <- data.table(
      num_nodes = num_nodes,
      num_leaves = num_leaves,
      num_nodes_tested = num_nodes_tested,
      num_nonnull_nodes_tested = num_nonnull_nodes_tested,
      node_rejections = node_rejections,
      node_any_false_rejection = node_any_false_rejection,
      node_false_rejection_prop = node_false_rejection_prop,
      node_num_false_rejections = node_num_false_rejections,
      node_false_discovery_prop = node_false_discovery_prop,
      node_true_discoveries = node_true_discoveries,
      node_power = node_power,
      num_leaves_tested = num_leaves_tested,
      num_nonnull_leaves_tested = num_nonnull_leaves_tested,
      leaf_power = leaf_power,
      leaf_rejections = leaf_rejections,
      leaf_true_discoveries = leaf_true_discoveries,
      leaf_any_false_rejection = leaf_any_false_rejection,
      leaf_false_rejection_prop = leaf_false_rejection_prop,
      leaf_false_discovery_prop = leaf_false_discovery_prop
    )
  }
  # One row of results
  #  detresults <- cbind(deterrs, detates)
  if (!return_details) {
    return(test_summary)
  } else {
    res <- list(
      detresults = deterrs,
      testobj = testobj,
      res = res
    )
    return(res)
  }
}
