# Test and develop functions to create tre,atment effects in simulations
# which blocks did it choose? Did it choose the right blocks?

# context("Simulation Functions")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(testthat)
  local_edition(3)
  library(here)
  library(data.table)
  library(dtplyr)
  library(dplyr)
  library(conflicted)
  conflicts_prefer(dplyr::filter)
  library(devtools)
  load_all() ## use  this during debugging
  load_all("~/repos/manytestsr")
  source(here("tests/testthat", "make_test_data.R"))
}

## For now bF is a character not a factor.
setDTthreads(1)
options(digits = 4)
#####
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]
## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]
bdat4[, bF := factor(bF)]
idat3[, bF := factor(bF)]
setkey(bdat4, bF)
setkey(bdat, bF)

set.seed(12355)
idat$y1test_null <-
  create_effects(
    idat = idat,
    ybase = "y0",
    blockid = "bF",
    tau_fn = tau_norm,
    tau_size = 0,
    prop_blocks_0 = .5
  )
idat$y1test_zeros <-
  create_effects(
    idat = idat,
    ybase = "y0",
    blockid = "bF",
    tau_fn = tau_norm,
    tau_size = 2,
    prop_blocks_0 = .5
  )
idat[, Y_zeros := y1test_zeros * Z + y0 * (1 - Z)]
idat[, Y_null := y1test_null * Z + y0 * (1 - Z)]

bdt1[, bF := as.factor(bF)]
idt[, bF := as.factor(bF)]
## This next is just a test of the reporting functions
set.seed(12345)
res_half <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = NULL,
  local_adj_p_fn = NULL, simthresh = 20, sims = 1000, maxtest = 2000,
  thealpha = .05, thew0 = 0.05 - 0.001,
  fmla = Y_half_tau1 ~ trtF | bF, splitby = "lvls_fac",
  blocksize = "nb", trace = TRUE, copydts = TRUE, stop_splitby_constant = TRUE,
  return_what = c("blocks", "nodes")
)

res_half_tree <- make_results_tree(res_half$bdat, block_id = "bF", truevar_name = "nonnull")
make_results_ggraph(res_half_tree$graph)
make_results_ggraph(res_half_tree$graph %>% filter(!is.na(p)))
res_half_results <- report_detections(res_half$bdat, blockid = "bF")

res_bdat4 <- find_blocks(
  fmla = Ynorm_dec ~ ZF | bF, idat = idat3,
  bdat = bdat4, pfn = pOneway, alphafn = NULL, blockid = "bF",
  splitfn = splitCluster, splitby = "hwt", local_adj_p_fn = NULL,
  thealpha = .05, blocksize = "hwt", trace = TRUE, copydts = TRUE,
  return_what = c("blocks", "nodes")
)

res_bdat4_tree <- make_results_tree(res_bdat4$bdat, block_id = "bF", truevar_name = "ate_norm_dec")
make_results_ggraph(res_bdat4_tree$graph)
## This next looks very weird because of all of the blocks that are not tested and so we don't exactly know where they are in the tree
## make_results_ggraph(res_bdat4_tree$graph, remove_na_p=FALSE)
res_bdat4_results <- report_detections(res_bdat4$bdat, blockid = "bF")

set.seed(12345)
res_half_investing <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = alpha_investing,
  local_adj_p_fn = NULL, simthresh = 20, sims = 1000, maxtest = 2000,
  thealpha = .05, thew0 = 0.05 - 0.001,
  fmla = Y_half_tau1 ~ trtF | bF, splitby = "lvls_fac",
  blocksize = "nb", trace = TRUE, copydts = TRUE, stop_splitby_constant = TRUE,
  return_what = c("blocks", "nodes")
)

res_half_investing_tree <- make_results_tree(res_half_investing$bdat, block_id = "bF", truevar_name = "nonnull")

res_half_investing_tree$nodes %>%
  mutate(across(where(is.numeric), zapsmall)) %>%
  select(-blocks)

test_that("Assess the error calculation function on the test in every block or bottom-up approach",{
res_half_bottom_up <- adjust_block_tests(idat = idt, bdat = bdt1, blockid = "bF", pfn = pOneway, p_adj_method = "hommel", fmla = Y_half_tau1 ~ trtF, copydts = TRUE)

res_half_bottom_up_errs <- calc_errs_new(res_half_bottom_up, truevar_name = "nonnull") %>% select(contains("lea"))
expect_equal(sum(res_half_bottom_up$nonnull), res_half_bottom_up_errs$num_nonnull_leaves_tested)
expect_equal(res_half_bottom_up_errs$leaf_rejections, sum(res_half_bottom_up[, max_p <= .05]))
expect_equal(res_half_bottom_up_errs$leaf_true_discoveries, sum(res_half_bottom_up[nonnull == TRUE, max_p <= .05]))
expect_equal(res_half_bottom_up_errs$leaf_power, mean(res_half_bottom_up[nonnull == TRUE, max_p <= .05]))
})

alpha_and_splits <- expand.grid(
  afn = c("alpha_investing", "alpha_saffron", "NULL"),
  sfn = c(
    "splitCluster",
    "splitEqualApprox",
    "splitLOO",
    "splitSpecifiedFactor",
    "splitSpecifiedFactorMulti"
  ),
  stringsAsFactors = FALSE
)
alpha_and_splits$splitby <- "hwt"
alpha_and_splits$splitby[grep("Specified", alpha_and_splits$sfn)] <- "lvs"

err_testing_fn <-
  function(afn,
           sfn,
           sby,
           local_adj_p_fn,
           fmla = Ynorm_dec ~ ZF | bF,
           idat = idat3,
           bdat = bdat4,
           truevar_name,
           blocksize = "hwt") {
    if (afn == "NULL") {
      theafn <- NULL
    } else {
      theafn <- get(afn)
    }

    thealpha <- .05
    trueeffect_tol <- .Machine$double.eps

    ## afn and sfn and sby are character names
    theres <- find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "bF",
      splitfn = get(sfn),
      local_adj_p_fn = get(local_adj_p_fn),
      pfn = pOneway,
      alphafn = theafn,
      thealpha = thealpha,
      fmla = fmla,
      parallel = "multicore",
      ncores = 2,
      copydts = TRUE,
      splitby = sby,
      blocksize = blocksize
    )
    detobj <-
      report_detections(theres$bdat, only_hits = FALSE, fwer = afn == "NULL", blockid = "bF")

    detobj[, hit := as.numeric(hit)]
    detobj[, hitb := as.numeric(max_p <= max_alpha & blocksbygroup == 1)]
    detobj[, hitb2 := as.numeric(single_hit)]
    stopifnot(all.equal(detobj$hitb, detobj$hitb2))
    ## Coding whether the true effect is zero or not by block.
    detobj[, true0 := as.numeric(abs(get(truevar_name)) <= trueeffect_tol)]
    detobj[, truenot0 := as.numeric(abs(get(truevar_name)) > trueeffect_tol)]
    ## Accessing the results another way
    thetree <- make_results_tree(theres$bdat, block_id = "bF", return_what = "all", truevar_name = truevar_name)

    ## Start recording key pieces of information:
    ### Number of blocks or possible leaves
    n_blocks <- nrow(theres$bdat)
    n_blocks2 <- nrow(detobj)
    expect_equal(n_blocks, n_blocks2)
    expect_equal(n_blocks, thetree$test_summary$num_leaves)

    ## Number of tests done total:
    thetree$test_summary$num_nodes_tested
    ## This next records all non-missing p-values (i.e. all nodes tested)
    expect_equal(
      thetree$test_summary$num_nodes_tested,
      nrow(theres$node_dat)
    )

    ## Number of leaves tested I think of a leaf as a single block --- a node
    ## with only one block This can happen at different depths in the tree
    ## depending on how we split the tree

    expect_equal(
      thetree$test_summary$num_leaves_tested,
      thetree$nodes %>%
        filter(num_leaves == 1 & !is.na(p)) %>%
        nrow()
    )
    expect_equal(theres$bdat %>% filter(blocksbygroup == 1) %>% nrow(), thetree$test_summary$num_leaves_tested)

    ## This next calculates errors and discoveries at the level of the block
    err_tab0 <- with(detobj, table(hitb, true0, exclude = c()))
    ## Make a table of rejections by hypotheses
    ##                      True 0            |  Not True 0 (actual effect)
    ## Not Reject   True Reject (not error)   |    False reject (low power error)
    ## Rejected        False positive (Error) |    Correct rejection (detection of the truth)
    if (!identical(dim(err_tab0), as.integer(c(2, 2)))) {
      blank_mat <-
        matrix(0, 2, 2, dimnames = list(c("0", "1"), c("0", "1")))
      blank_mat[rownames(err_tab0), colnames(err_tab0)] <- err_tab0
      err_tab <- blank_mat
    } else {
      err_tab <- err_tab0
    }
    errs <- calc_errs(theres$bdat,
      truevar_name = truevar_name,
      blockid = "bF",
      trueeffect_tol = .Machine$double.eps, fwer = afn == "NULL"
    )
    tottests <- sum(err_tab)
    tottests2 <- errs$nreject + errs$naccept
    expect_equal(tottests, tottests2)

    if (thetree$test_summary$num_leaves_tested > 0) {
      expect_equal(as.numeric(thetree$test_summary$leaf_rejections), errs$nreject)
    } else {
      expect_equal(thetree$test_summary$leaf_rejections==0, errs$nreject == 0)
    }

    expect_equal(errs$true_pos_prop, err_tab["1", "0"] / tottests)
    expect_equal(errs$prop_not_reject_true0, err_tab["0", "1"] / tottests)
    expect_equal(errs$false_pos_prop, err_tab["1", "1"] / tottests)
    expect_equal(errs$false_neg_prop, err_tab["0", "0"] / tottests)
    expect_equal(errs$true_disc_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))
    expect_equal(errs$false_disc_prop, err_tab["1", "1"] / max(1, sum(err_tab["1", ])))
    expect_equal(errs$true_nondisc_prop, err_tab["0", "1"] / max(1, sum(err_tab["0", ])))
    expect_equal(errs$false_nondisc_prop, err_tab["0", "0"] / max(1, sum(err_tab["0", ])))
  }


err_testing_bottom_up <-
  function(fmla = Ynorm_dec ~ ZF,
           idat = idat3,
           bdat = bdat4,
           truevar_name,
           thealpha = .05) {
    ## afn and sfn and sby are character names
    theres <- adjust_block_tests(
      idat = idat,
      bdat = bdat,
      blockid = "bF",
      p_adj_method = "BH", ## FDR adjustment
      pfn = pOneway,
      fmla = fmla,
      copydts = TRUE,
      parallel = "multicore",
      ncores = 2
    )

    theres[, hit := max_p < .05] ## doing FDR==.05 for now. All that really matters is some fixed number.

    blocks <- theres[, .(
      bF = bF,
      hit = as.numeric(hit),
      true0 = as.numeric(get(truevar_name) == 0),
      truenot0 = as.numeric(get(truevar_name) == 0),
      anynotnull = as.numeric(get(truevar_name) != 0)
    )]

    err_tab0 <- with(blocks, table(rejected = hit, non_zero_effect = anynotnull, exclude = c()))
    #err_tab0 <- with(blocks, table(rejected = hit, true0 = true0, exclude = c()))
    ## Make a table of rejections by hypotheses
    ##                      True 0      "0"      |  Not True 0 (actual effect) "1"
    ## Not Reject "0"  True Reject (not error)   |    False reject (low power error)
    ## Rejected "1"       False positive (Error) |    Correct rejection (detection of the truth)

    if (!identical(dim(err_tab0), as.integer(c(2, 2)))) {
      blank_mat <-
        matrix(0, 2, 2, dimnames = list(c("0", "1"), c("0", "1")))
      blank_mat[rownames(err_tab0), colnames(err_tab0)] <- err_tab0
      err_tab <- blank_mat
    } else {
      err_tab <- err_tab0
    }

## So:
  ## err_tab["1","1"] are the number of correct rejections
  ## err_tab["0","0"] are the number of correct not-rejections
  ## err_tab["1","0"] are the number of false positive rejections (rejecting a null effect)
  ## err_tab["0","1"] are the number of false negative rejections (failing to reject a null effect) (a measure of power)
  ## power would be 1 - (err_tab["0","1"]/sum(err_tab[,"1"])) --- 1 - proportion of non-null effects not rejected=proportion of non-effects rejected
  ## or err_tab["1","1"]/sum(err_tab[,"1"])

    errs <- calc_errs_new(
      testobj = theres,
      truevar_name = truevar_name,
      trueeffect_tol = .Machine$double.eps, fwer = FALSE
    )

    tottests <- sum(err_tab)
    tottests_from_calc_errs <- errs$num_leaves_tested
    expect_equal(tottests, tottests_from_calc_errs)
  ## This approach has to test in all leaves
    expect_equal(tottests, errs$num_leaves)

  expect_equal(errs$leaf_rejections,sum(err_tab["1",]))
    expect_equal(errs$leaf_true_discoveries, err_tab["1", "1"])
  if(!is.na(errs$leaf_power)){
    expect_equal(errs$leaf_power, err_tab["1", "1"] / max(1, sum(err_tab[,"1"])))
  }
    expect_equal(errs$leaf_false_discovery_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))

    if(!is.na(errs$leaf_false_rejection_prop)){
    expect_equal(errs$leaf_false_rejection_prop, err_tab["1", "0"] / tottests)
  }

  }

err_testing_fn(
  fmla = Ynull ~ ZF | bF, idat = idat3,
  bdat = bdat4, truevar_name = "ate_null", afn = "NULL",
  sfn = "splitCluster", sby = "hwt", local_adj_p_fn = "local_unadj_all_ps"
)

err_testing_fn(
  fmla = Ynorm_dec ~ ZF | bF, idat = idat3,
  bdat = bdat4, truevar_name = "ate_norm_dec", afn = "NULL",
  sfn = "splitCluster", sby = "hwt", local_adj_p_fn = "local_unadj_all_ps"
)

err_testing_fn(
  fmla = Y_half_tau1 ~ trtF | bF, idat = idt,
  bdat = bdt1, truevar_name = "nonnull", afn = "NULL", blocksize = "nb",
  sfn = "splitSpecifiedFactorMulti", sby = "lvls_fac", local_adj_p_fn = "local_unadj_all_ps"
)

err_testing_fn(
  afn = "alpha_investing",
  sfn = "splitCluster",
  sby = "hwt",
  idat = idat3,
  bdat = bdat4,
  fmla = Yhomog ~ ZF | bF,
  truevar_name = "ate_homog",
  local_adj_p_fn = "local_hommel_all_ps"
)

err_testing_bottom_up(fmla = Ynorm_dec ~ ZF, idat = idat3, bdat = bdat4, truevar_name = "ate_norm_dec")
err_testing_bottom_up(fmla = Ynorm_dec ~ ZF, idat = idat3, bdat = bdat4, truevar_name = "ate_tau")

resnms <- apply(alpha_and_splits, 1, function(x) {
  paste(x, collapse = "_", sep = "")
})


test_that("Error calculations for a given set of tests work: No effects at all", {
  ### No effects at all
  test_lst <- mapply(
    FUN = function(afn = afn,
                   sfn = sfn,
                   sby = sby,
                   truevar_name = "ate_null") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn,
        sfn = sfn,
        sby = sby,
        idat = idat3,
        bdat = bdat4,
        fmla = Ynull ~ ZF | bF,
        truevar_name = truevar_name,
        local_adj_p_fn = "local_unadj_all_ps"
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby,
    SIMPLIFY = FALSE
  )
})

test_that("Error calculations for a given set of tests work: large and homogenous effects", {
  test_lst <- mapply(
    FUN = function(afn = afn,
                   sfn = sfn,
                   sby = sby,
                   truevar_name = "ate_homog") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn,
        sfn = sfn,
        sby = sby,
        idat = idat3,
        bdat = bdat4,
        fmla = Yhomog ~ ZF | bF,
        truevar_name = truevar_name,
        local_adj_p_fn = "local_unadj_all_ps"
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby,
    SIMPLIFY = FALSE
  )
})


test_that(
  "Error calculations for a given set of tests work: individually heteogeneous effects and increase with block size. Also some completely null blocks.",
  {
    test_lst <- mapply(
      FUN = function(afn = afn,
                     sfn = sfn,
                     sby = sby,
                     truevar_name = "ate_norm_inc") {
        message(paste(afn, sfn, sby, collapse = ","))
        err_testing_fn(
          afn = afn,
          sfn = sfn,
          sby = sby,
          idat = idat3,
          bdat = bdat4,
          fmla = Ynorm_inc ~ ZF | bF,
          truevar_name = truevar_name,
          local_adj_p_fn = "local_unadj_all_ps"
        )
      },
      afn = alpha_and_splits$afn,
      sfn = alpha_and_splits$sfn,
      sby = alpha_and_splits$splitby,
      SIMPLIFY = FALSE
    )
  }
)

test_that(
  "Error calculations for a given set of tests work:individually heteogeneous effects and decrease with block size. Also some completely null blocks.",
  {
    test_lst <- mapply(
      FUN = function(afn = afn,
                     sfn = sfn,
                     sby = sby,
                     truevar_name = "ate_norm_dec") {
        message(paste(afn, sfn, sby, collapse = ","))
        err_testing_fn(
          afn = afn,
          sfn = sfn,
          sby = sby,
          idat = idat3,
          bdat = bdat4,
          fmla = Ynorm_dec ~ ZF | bF,
          truevar_name = truevar_name,
          local_adj_p_fn = "local_unadj_all_ps"
        )
      },
      afn = alpha_and_splits$afn,
      sfn = alpha_and_splits$sfn,
      sby = alpha_and_splits$splitby,
      SIMPLIFY = FALSE
    )
  }
)

test_that(
  "Error calculations for a given set of tests work:constant effect that cancel out at the high level (half large and positive, half large and negative).",
  {
    test_lst <- mapply(
      FUN = function(afn = afn,
                     sfn = sfn,
                     sby = sby,
                     truevar_name = "ate_tau") {
        message(paste(afn, sfn, sby, collapse = ","))
        err_testing_fn(
          afn = afn,
          sfn = sfn,
          sby = sby,
          idat = idat3,
          bdat = bdat4,
          fmla = Y ~ ZF | bF,
          truevar_name = truevar_name,
          local_adj_p_fn = "local_unadj_all_ps"
        )
      },
      afn = alpha_and_splits$afn,
      sfn = alpha_and_splits$sfn,
      sby = alpha_and_splits$splitby,
      SIMPLIFY = FALSE
    )
  }
)

test_that(
  "Error calculations for a given set of tests work:individually heterogeneous effects block-fixed effects and some completely null blocks.",
  {
    test_lst <- mapply(
      FUN = function(afn = afn,
                     sfn = sfn,
                     sby = sby,
                     truevar_name = "ate_tauv2") {
        message(paste(afn, sfn, sby, collapse = ","))
        err_testing_fn(
          afn = afn,
          sfn = sfn,
          sby = sby,
          idat = idat3,
          bdat = bdat4,
          fmla = Ytauv2 ~ ZF | bF,
          truevar_name = truevar_name,
          local_adj_p_fn = "local_unadj_all_ps"
        )
      },
      afn = alpha_and_splits$afn,
      sfn = alpha_and_splits$sfn,
      sby = alpha_and_splits$splitby,
      SIMPLIFY = FALSE
    )
  }
)

test_that(
  "Error calculations for a given set of tests work: Testing in every block and adjusting p-values via FDR/BH.",
  {
    truevar_names <- grep("^ate", names(bdat4), value = TRUE)
    outcome_names <- paste0("Y", gsub("ate_", "", truevar_names))
    outcome_names[1] <- "Y"
    stopifnot(all(outcome_names %in% names(idat3)))

    test_lst <- mapply(
      FUN = function(truevar_name, outcome_name) {
        fmla <- reformulate("ZF", response = outcome_name)
        message(paste(truevar_name, collapse = ","))
        err_testing_bottom_up(
          idat = idat3,
          bdat = bdat4,
          fmla = fmla,
          truevar_name = truevar_name
        )
      },
      truevar_name = truevar_names,
      outcome_name = outcome_names,
      SIMPLIFY = FALSE
    )
  }
)


test_that("Some simulation code runs without error",{
  skip()
  skip_on_ci()
  skip_on_cran()

### This next is less of a test with expected results and more to ensure that the code runs without error.
simparms <- cbind(alpha_and_splits, p_adj_method = rep("split", nrow(alpha_and_splits)))
simparms <- rbind(simparms, c("NULL", "NULL", "NULL", "fdr"))
simparms <- data.table(simparms)
simresnms <- apply(simparms, 1, function(x) {
  paste(x, collapse = "_", sep = "")
})

set.seed(12345)
p_sims_res <- lapply(
  seq_len(nrow(simparms)),
  FUN = function(i) {
    x <- simparms[i, ]
    xnm <- paste(x, collapse = "_")
    message(xnm)
    nsims <- 10 ## 100
    p_sims_tab <- padj_test_fn(
      idat = idat3,
      bdat = bdat4,
      blockid = "bF",
      trtid = "Z",
      fmla = Y ~ ZF | bF,
      ybase = "y0",
      prop_blocks_0 = .5,
      tau_fn = tau_norm_covariate_outliers,
      tau_size = .5,
      covariate = "v4",
      pfn = pOneway,
      nsims = nsims,
      afn = ifelse(x[["afn"]] != "NULL", getFromNamespace(x[["afn"]], ns = "manytestsr"), "NULL"),
      p_adj_method = x[["p_adj_method"]],
      splitfn = ifelse(x[["sfn"]] != "NULL", getFromNamespace(x[["sfn"]], ns = "manytestsr"), "NULL"),
      splitby = x[["splitby"]],
      local_adj_p_fn = local_unadj_all_ps,
      bottom_up_adj = "hommel",
      return_bottom_up = TRUE,
      return_details = FALSE,
      ncores = 2
    )
    return(p_sims_tab)
  }
)

names(p_sims_res) <- simresnms

lapply(p_sims_res, function(obj) {
  obj[, lapply(.SD, mean)]
})


p_sims_obj <- rbindlist(p_sims_res[1:15], idcol = TRUE)

})


# test_that(
#   "Clustering-based Splitters control FWER.",
#   {
#     simparms2 <- data.table(expand.grid(
#       splitby = c("hwt", "v4", "newcov"),
#       covariate = c("v4", "newcov"), stringsAsFactors = FALSE
#     ))
#     simparms2[, afn := "NULL"]
#     simparms2[, sfn := "splitCluster"]
#     simparms2[, p_adj_method := "split"]
#     skip_on_ci()
#     skip_on_cran()
#     set.seed(12345)
#     res2 <- lapply(
#       seq_len(nrow(simparms2)),
#       FUN = function(i) {
#         x <- simparms2[i, ]
#         xnm <- paste(x, collapse = "_")
#         message(xnm)
#         nsims <- 100
#         p_sims_tab <- padj_test_fn(
#           idat = idat3,
#           bdat = bdat4,
#           blockid = "bF",
#           trtid = "Z",
#           fmla = Y ~ ZF | blockF,
#           ybase = "y0",
#           prop_blocks_0 = 1,
#           tau_fn = tau_norm_covariate_outliers,
#           tau_size = 0,
#           covariate = x[["covariate"]],
#           pfn = pIndepDist,
#           nsims = nsims,
#           afn = "NULL",
#           p_adj_method = "split",
#           splitfn = getFromNamespace(x[["sfn"]], ns = "manytestsr"),
#           splitby = x[["splitby"]],
#           ncores = 1
#         )
#         return(p_sims_tab)
#       }
#     )
#     p_sims_obj2 <- rbindlist(res2, idcol = TRUE)
#     err_rates2 <- p_sims_obj2[, lapply(.SD, mean), .SDcols = c(
#       "true_pos_prop",
#       "false_pos_prop",
#       "prop_not_reject_true0",
#       "false_neg_prop",
#       "true_disc_prop",
#       "false_disc_prop",
#       "true_nondisc_prop",
#       "false_nondisc_prop"
#     ), by = .id]
#     err_rates2
#     expect_lte(max(err_rates2$false_pos_prop), .05 + 2 * (sqrt(.025 / 100)))
#   }
# )
