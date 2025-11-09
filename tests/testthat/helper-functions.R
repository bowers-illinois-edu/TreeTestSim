
err_testing_fn <-
  function(afn, sfn, sby, local_adj_p_fn,
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
    thetree <- make_results_tree(theres, block_id = "bF", return_what = "all", truevar_name = truevar_name)

    ## Start recording key pieces of information:
    ### Number of blocks or possible leaves
    n_blocks <- nrow(theres$bdat)
    n_blocks2 <- nrow(detobj)
    expect_equal(n_blocks, n_blocks2)
    expect_equal(n_blocks, thetree$test_summary$num_blocks)

    ## Number of tests done total:
    thetree$test_summary$num_nodes_tested
    ## This next records all non-missing p-values (i.e. all nodes tested)
    expect_equal(
      thetree$test_summary$num_nodes_tested,
      nrow(theres$node_dat)
    )

    ## Number of leaves tested

  ## A leaf is technically a node with no children. Sometimes
  ## these are single blocks and sometimes they are collections
  ## of blocks.

    expect_equal(
      thetree$test_summary$num_leaves_tested,
      thetree$nodes %>%
        filter(node_type=="leaf" & !is.na(p)) %>%
        nrow()
    )
    expect_equal(theres$bdat %>% filter(blocksbygroup == 1) %>% nrow(), 
    sum(thetree$nodes$num_blocks==1,na.rm=TRUE))

    ## This next calculates errors and discoveries at the level of the block
    err_tab0 <- with(detobj, table(hitb, true0)) #, exclude = c()))
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

  errs2 <- calc_errs_new(theres,
      truevar_name = truevar_name,
      blockid = "bF",
      trueeffect_tol = .Machine$double.eps, fwer = afn == "NULL"
    )

  ## calc_errs is deprecated. Delete the below soon
   # errs <- calc_errs(theres$bdat,
   #   truevar_name = truevar_name,
   #   blockid = "bF",
   #   trueeffect_tol = .Machine$double.eps, fwer = afn == "NULL"
   # )
   ## tottests <- sum(err_tab)
   ## tottests2 <- errs$nreject + errs$naccept
   ## expect_equal(tottests, tottests2)

   ## if (thetree$test_summary$num_leaves_tested > 0) {
   ##   expect_equal(as.numeric(thetree$test_summary$leaf_rejections), errs$nreject)
   ## } else {
   ##   expect_equal(thetree$test_summary$leaf_rejections == 0, errs$nreject == 0)
   ## }

   ## expect_equal(errs$true_pos_prop, err_tab["1", "0"] / tottests)
   ## expect_equal(errs$prop_not_reject_true0, err_tab["0", "1"] / tottests)
   ## expect_equal(errs$false_pos_prop, err_tab["1", "1"] / tottests)
   ## expect_equal(errs$false_neg_prop, err_tab["0", "0"] / tottests)
   ## expect_equal(errs$true_disc_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))
   ## expect_equal(errs$false_disc_prop, err_tab["1", "1"] / max(1, sum(err_tab["1", ])))
   ## expect_equal(errs$true_nondisc_prop, err_tab["0", "1"] / max(1, sum(err_tab["0", ])))
   ## expect_equal(errs$false_nondisc_prop, err_tab["0", "0"] / max(1, sum(err_tab["0", ])))
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
    # err_tab0 <- with(blocks, table(rejected = hit, true0 = true0, exclude = c()))
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
    tottests_from_calc_errs <- errs$num_blocks
    expect_equal(tottests, tottests_from_calc_errs)

    expect_equal(errs$leaf_rejections, sum(err_tab["1", ]))
    expect_equal(errs$leaf_true_discoveries, err_tab["1", "1"])
    if (!is.na(errs$leaf_power)) {
      expect_equal(errs$leaf_power, err_tab["1", "1"] / max(1, sum(err_tab[, "1"])))
    }
    expect_equal(errs$leaf_false_discovery_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))

    if (!is.na(errs$leaf_false_rejection_prop)) {
      expect_equal(errs$leaf_false_rejection_prop, err_tab["1", "0"] / tottests)
    }
  }

