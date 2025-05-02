# Test and develop functions to create treatment effects in simulations
# which blocks did it choose? Did it choose the right blocks?

# context("Simulation Functions")


## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(dtplyr)
  library(dplyr)
  library(conflicted)
  conflicts_prefer(dplyr::filter)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  load_all() ## use  this during debugging
}


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

### Add testthat here for these two to assess the create_effects function
## Should not reject
oneway_test(Y_null ~ ZF | bF, data = idat)
## Should reject: half blocks are 0 but they are big enough and the effects are large enough

alpha_and_splits <- expand.grid(
  afn = c("alpha_investing", "alpha_saffron", "NULL"),
  sfn = c(
    "splitCluster",
    "splitEqualApprox",
    "splitLOO",
    "splitSpecifiedFactor"
  ),
  stringsAsFactors = FALSE
)
alpha_and_splits$splitby <- "hwt"
alpha_and_splits$splitby[grep("Specified", alpha_and_splits$sfn)] <-
  "lvs"
oneway_test(Y_zeros ~ ZF | bF, data = idat)
####

err_testing_fn <-
  function(afn,
           sfn,
           sby,
           fmla = Ytauv2 ~ ZF | bF,
           idat = idat3,
           bdat = bdat4,
           truevar_name) {
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
      pfn = pIndepDist,
      alphafn = theafn,
      thealpha = thealpha,
      fmla = fmla,
      parallel = "multicore",
      ncores = 2,
      copydts = TRUE,
      splitby = sby
    )
    detobj <-
      report_detections(theres, only_hits = FALSE, fwer = afn == "NULL", blockid = "bF")

    detobj[, hit := as.numeric(hit)]
    detobj[, hitb := as.numeric(max_p <= max_alpha & blocksbygroup == 1)]
    detobj[, hitb2 := as.numeric(single_hit)]
    stopifnot(all.equal(detobj$hitb, detobj$hitb2))
    ## Coding whether the true effect is zero or not by block.
    detobj[, true0 := as.numeric(abs(get(truevar_name)) <= trueeffect_tol)]
    detobj[, truenot0 := as.numeric(abs(get(truevar_name)) > trueeffect_tol)]
    ## Accessing the results another way
    thetree <- make_results_tree(theres, blockid = "bF") %>%
      select(-label) %>%
      as.data.frame()
    sigleaves <- thetree %>% filter(out_degree == 0 & hit == 1)
    ##  detobj[blockF %in% sigleaves$bF,.(blockF,max_p,single_hit,hitb,hitb2,fin_grp)]
    stopifnot(all.equal(sort(sigleaves$bF), sort(as.character(detobj[hitb == 1, bF]))))

    ## nodes <- detects[, .(
    ##   hit = as.numeric(unique(hit)),
    ##   numnull = sum(get(truevar_name) == 0),
    ##   anynull = as.numeric(any(get(truevar_name) == 0)),
    ##   allnull = as.numeric(all(get(truevar_name) == 0)),
    ##   numnotnull = sum(get(truevar_name) != 0),
    ##   anynotnull = as.numeric(any(get(truevar_name) != 0))
    ## ), by = hit_grp]
    # err_tab1 <- with(nodes1 , table(hit, anynotnull, exclude = c()))
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
    errs <- calc_errs(theres,
      truevar_name = truevar_name,
      trueeffect_tol = .Machine$double.eps, fwer = afn == "NULL"
    )
    tottests <- sum(err_tab)
    tottests2 <- errs$nreject + errs$naccept
    expect_equal(tottests, tottests2)
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
  function(fmla = Ytauv2 ~ ZF,
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
      pfn = pIndepDist,
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

    err_tab0 <- with(blocks, table(rejected = hit, true0 = true0, exclude = c()))
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

    errs <- calc_errs(
      testobj = theres,
      truevar_name = truevar_name,
      trueeffect_tol = .Machine$double.eps, fwer = FALSE
    )

    tottests <- sum(err_tab)
    tottests_from_calc_errs <- errs$nreject + errs$naccept
    expect_equal(tottests, tottests_from_calc_errs)

    expect_equal(errs$true_pos_prop, err_tab["1", "0"] / tottests)
    expect_equal(errs$tot_not_reject_true0 / tottests, err_tab["0", "1"] / tottests)
    expect_equal(errs$false_pos_prop, err_tab["1", "1"] / tottests)
    expect_equal(errs$false_neg_prop, err_tab["0", "0"] / tottests)
    expect_equal(errs$true_disc_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))
    expect_equal(errs$false_disc_prop, err_tab["1", "1"] / max(1, sum(err_tab["1", ])))
    expect_equal(errs$true_nondisc_prop, err_tab["0", "1"] / max(1, sum(err_tab["0", ])))
    expect_equal(errs$false_nondisc_prop, err_tab["0", "0"] / max(1, sum(err_tab["0", ])))
  }


err_testing_fn(
  fmla = Ynull ~ ZF | bF, idat = idat3,
  bdat = bdat4, truevar_name = "ate_null", afn = "NULL", sfn = "splitCluster", sby = "hwt"
)


err_testing_fn(
  afn = "alpha_saffron",
  sfn = "splitEqualApprox",
  sby = "hwt",
  idat = idat3,
  bdat = bdat4,
  fmla = Ynorm_dec ~ ZF | bF,
  truevar_name = "ate_norm_dec"
)

err_testing_bottom_up(fmla = Ytauv2 ~ ZF, idat = idat3, bdat = bdat4, truevar_name = "ate_tauv2")

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
        truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby,
    SIMPLIFY = FALSE
  )
  ## names(tau_null) <- resnms
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
        truevar_name = truevar_name
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
          truevar_name = truevar_name
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
          truevar_name = truevar_name
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
          truevar_name = truevar_name
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
          truevar_name = truevar_name
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


### This next is less of a test with expected results and more to ensure that the code runs without error.
simparms <- cbind(alpha_and_splits, p_adj_method = rep("split", nrow(alpha_and_splits)))
simparms <- rbind(simparms, c("NULL", "NULL", "NULL", "fdr"))
simparms <- data.table(simparms)
simresnms <- apply(simparms, 1, function(x) {
  paste(x, collapse = "_", sep = "")
})

set.seed(12345)
p_sims_res <- lapply(
  1:nrow(simparms),
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
      fmla = Y ~ ZF | blockF,
      ybase = "y0",
      prop_blocks_0 = .5,
      tau_fn = tau_norm_covariate_outliers,
      tau_size = .5,
      covariate = "v4",
      pfn = pIndepDist,
      nsims = nsims,
      afn = ifelse(x[["afn"]] != "NULL", getFromNamespace(x[["afn"]], ns = "manytestsr"), "NULL"),
      p_adj_method = x[["p_adj_method"]],
      splitfn = ifelse(x[["sfn"]] != "NULL", getFromNamespace(x[["sfn"]], ns = "manytestsr"), "NULL"),
      splitby = x[["splitby"]],
      ncores = 2
    )
    # p_sims_tab <- p_sims_tab[1,,]
    return(p_sims_tab)
  }
)

names(p_sims_res) <- simresnms

p_sims_obj <- rbindlist(p_sims_res, idcol = TRUE)

err_rates <- p_sims_obj[, lapply(.SD, mean), .SDcols = c(
  "true_pos_prop",
  "false_pos_prop",
  "correct_pos_effect_prop",
  "false_neg_prop",
  "true_disc_prop",
  "false_disc_prop",
  "true_nondisc_prop",
  "false_nondisc_prop"
), by = .id]

## This is just one run --- so these are not really error rates but error proportions
err_rates

test_that(
  "Clustering-based Splitters control FWER.",
  {
    simparms2 <- data.table(expand.grid(
      splitby = c("hwt", "v4", "newcov"),
      covariate = c("v4", "newcov"), stringsAsFactors = FALSE
    ))
    simparms2[, afn := "NULL"]
    simparms2[, sfn := "splitCluster"]
    simparms2[, p_adj_method := "split"]
    set.seed(12345)
    res2 <- lapply(
      seq_len(nrow(simparms2)),
      FUN = function(i) {
        x <- simparms2[i, ]
        xnm <- paste(x, collapse = "_")
        message(xnm)
        nsims <- 100
        p_sims_tab <- padj_test_fn(
          idat = idat3,
          bdat = bdat4,
          blockid = "bF",
          trtid = "Z",
          fmla = Y ~ ZF | blockF,
          ybase = "y0",
          prop_blocks_0 = 1,
          tau_fn = tau_norm_covariate_outliers,
          tau_size = 0,
          covariate = x[["covariate"]],
          pfn = pIndepDist,
          nsims = nsims,
          afn = "NULL",
          p_adj_method = "split",
          splitfn = getFromNamespace(x[["sfn"]], ns = "manytestsr"),
          splitby = x[["splitby"]],
          ncores = 1
        )
        return(p_sims_tab)
      }
    )
    p_sims_obj2 <- rbindlist(res2, idcol = TRUE)
    err_rates2 <- p_sims_obj2[, lapply(.SD, mean), .SDcols = c(
      "true_pos_prop",
      "false_pos_prop",
      "prop_not_reject_true0",
      "false_neg_prop",
      "true_disc_prop",
      "false_disc_prop",
      "true_nondisc_prop",
      "false_nondisc_prop"
    ), by = .id]
    err_rates2
    expect_lte(max(err_rates2$false_pos_prop), .05 + 2 * (sqrt(.025 / 100)))
  }
)
