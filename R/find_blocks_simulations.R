# Functions to make learning about the properties of tests easier

#' Error assessment function
#'
#' This function produces p-values evaluating hypotheses of known forms to help
#' assess the operating characteristics of testing procedures. It creates causal effects of known forms and calls the
#' findBlocks function and passes objects and arguments onto it.
#'
#' @param idat Unit-level data. An object inheriting from class data.table
#' @param bdat Block-level data. An object inheriting from class data.table
#' @param pfn  A function to  produce p-values (see [pWilcox] for example).
#' @param splitfn A function for splitting the data, see [splitLOO] for example.
#' @param fmla A formula relating outcomes  to  treatment and blocks (to  be passed to  pfn).
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtvar Is the name of the treatment numeric, (0,1), variable
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param afn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
#' @param thealpha The alpha level of the test.
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param prop_blocks_0 The proportion of blocks having no treatment effect at all.
#' @param sims Is the number of simulations to run --- each simulation uses the same treatment effects be re-assigns treatment (re-shuffles treatment and re-reveals the observed outcomes as a function of the potential outcomes)
#' @return A data.table with final pvalues, the associated blocks,  the true treatment effects, and the order in which the tests were conducted.
#' @export
errsimfn <- function(idat, bdat, pfn, splitfn,
                     fmla = Y ~ newtrtF | blockF,
                     blockid,
                     trtvar,
                     ybase,
                     afn = NULL,
                     thealpha = 0.05, tau_fn, tau_size = 0, prop_blocks_0 = 0, sims) {
  if (!is.null(afn) & is.character(afn)) {
    if (afn == "NULL") {
      afn <- NULL
    } else {
      afn <- get(afn)
    }
  }
  # require(data.table) # not necessary now that we have  a package
  idatnew <- copy(idat) # to avoid updating the bdat and idat inputs outside of the function
  bdatnew <- copy(bdat)
  setkeyv(bdatnew, blockid)
  setkeyv(idatnew, blockid)

  idatnew$y1new <- create_effects(idat = idatnew, ybase = ybase, blockid = blockid, tau_fn = tau_fn, tau_size = tau_size, prop_blocks_0 = prop_blocks_0)
  idatnew[, Y := y1new * get(trtvar) + get(ybase) * (1 - get(trtvar))]
  idatnew[, truetaui := y1new - get(ybase)]
  tausb <- idatnew[, .(
    truetaub = mean(truetaui),
    mndiffb = mean(Y[get(trtvar) == "1"]) - mean(Y[get(trtvar) == "0"])
  ), by = blockid]
  stopifnot(key(tausb) == blockid)
  res <- findBlocks(
    idat = idatnew, bdat = bdatnew, pfn = pfn, alphafn = afn, splitfn = splitfn, thealpha = thealpha,
    fmla = fmla, parallel = "no", blockid = blockid, sims = sims
  )
  setkeyv(res, blockid)
  res <- merge(res, tausb)
  # return(unique(res$pfinalb))
  # stopifnot(all.equal(res$block,bdatnew$block))
  # return(res[, .(pfinalb, block = get(blockid), mndiffb, truetaub, biggrp)])
  return(res)
}


#' Simulate Treatment Effects
#'
#' This function creates treatment effects of known forms using existing potential outcomes.
#' We specify differences across blocks and within blocks
#'
#' @param idat Unit-level data. An object inheriting from class data.table
#' @param ybase Baseline potential outcome variable name, a string.
#' @param blockid Name of block variable (the blocking variable is a factor)
#' @param tau_fn A function taking y0 and tau_size to produce a treatment effect for each unit
#' @param tau_size The rough or mean size of the treatment effect (like the mean shift), in sds probably.
#' @param covariate Contains information about covariates currently a character name of a column in idat. It is NULL if not used.
#' @param prop_blocks_0 The proportion of blocks having zero treatment effect.
#' @param by_block TRUE if each block has a separate shift (like a mean shift) or to otherwise create the individual level treatment effect separately within block or FALSE to create the effect across all blocks (ignoring blocks).
#' @return A data.table with a single column containing potential outcome to treatment.
#' @export
create_effects <- function(idat, ybase, blockid, tau_fn, tau_size, covariate = NULL, prop_blocks_0 = 0, by_block = TRUE) {
  if (!is.null(covariate) & is.character(covariate)) {
    if (covariate == "" | covariate == "NULL") {
      covariate <- NULL
    }
  }
  stopifnot(is.null(covariate) | is.vector(covariate))
  # assume that ybase already has block-level shifts in it
  idatnew <- copy(idat) # to avoid updating the bdat and idat inputs outside of the function
  stopifnot(prop_blocks_0 >= 0 & prop_blocks_0 <= 1)
  setkeyv(idatnew, blockid)

  if (is.null(covariate)) {
    if (by_block) {
      idatnew[, y1sim := get(ybase) + tau_fn(get(ybase), tau_size), by = blockid]
    } else {
      idatnew[, y1sim := get(ybase) + tau_fn(get(ybase), tau_size)]
    }
  } else {
    if (by_block) {
      idatnew[, y1sim := get(ybase) + tau_fn(get(ybase), tau_size, get(covariate)), by = blockid]
    } else {
      idatnew[, y1sim := get(ybase) + tau_fn(get(ybase), tau_size, get(covariate))]
    }
  }

  ## Choose blocks at random to set to zero
  blocks <- sort(as.character(unique(idatnew[[blockid]])))
  num_blocks <- length(blocks)
  n_null_blocks <- round(num_blocks * prop_blocks_0)
  if (prop_blocks_0 > 0 & n_null_blocks >= 1) {
    null_blocks <- sort(as.character(sample(blocks, size = n_null_blocks)))
    idatnew[.(null_blocks), y1sim := get(ybase)]
  }
  return(idatnew$y1sim)
}


#' Causal Effect Functions
#'
#' These functions create individual level causal effects, tau_i, that can be combined with a potential outcome to control to create a potential outcome to treatment.
#' @param ybase Is a vector of potential outcome to control
#' @param tau_sds is the number of sds of ybase to shift the distribution
#' @param covariate Contains information about covariates currently a character name of a column in idat. It is NULL if not used. Mostly it is a vector the same length as ybase.
#' @return A vector of individual level causal effects (taus) that we will add to ybase (potential outcome to control) to get y1var or potential outcome to treatment.

#' @describeIn Tau_Functions A basic function with no outliers
#' @export
tau_norm <- function(ybase, tau_sds, covariate) {
  ## covariate ignored here.
  n <- length(ybase)
  rnorm(n, mean = sd(ybase) * tau_sds, sd = sd(ybase) / 2)
}

#' @describeIn Tau_Functions Draws from a Normal but also adds a few outliers.
#' @export
tau_norm_outliers <- function(ybase, tau_sds, covariate) {
  ## covariate ignored here.
  n <- length(ybase)
  num_outliers <- round(length(ybase) * .02)
  if (num_outliers < 1) {
    rnorm(n, mean = tau_sds * sd(ybase), sd = sd(ybase) / 2)
  } else {
    # 2 large outliers
    c(
      rnorm(n - num_outliers, mean = sd(ybase) * tau_sds, sd = sd(ybase) / 2),
      rnorm(num_outliers, mean = sd(ybase) * tau_sds * 4, sd = sd(ybase) * 4)
    )
  }
}

#' @describeIn Tau_Functions A basic function that specifies a tau_sds*2 size effect if cov>median(cov) and otherwise is a tau_sds/2 size effect. The idea is to keep the average individual effect size the same as other functions --- i.e. about tau_sds --- but to make a strong but simple relationship with a covariate.
#' @export
tau_norm_covariate <- function(ybase, tau_sds, covariate) {
  stopifnot(length(covariate) == length(ybase))
  n <- length(ybase)
  thetau <- rnorm(n, mean = sd(ybase) * tau_sds, sd = sd(ybase) / 2)
  ifelse(covariate > median(covariate), thetau * 2, thetau / 2)
}

#' @describeIn Tau_Functions A basic function that specifies a tau_sds size effect if cov>median(cov) and otherwise is a tau_sds/4 size effect
#' @export
tau_norm_covariate_outliers <- function(ybase, tau_sds, covariate) {
  stopifnot(length(covariate) == length(ybase))
  thetau <- tau_norm_outliers(ybase = ybase, tau_sds = tau_sds, covariate = covariate)
  ifelse(covariate > median(covariate), thetau * 2, thetau / 2)
}

#' @describeIn Tau_Functions A basic function that specifies a tau_sds size effect if cov>median(cov) and otherwise is a tau_sds/4 size effect
#' @export
tau_norm_covariate_cont <- function(ybase, tau_sds, covariate) {
  stopifnot(length(covariate) == length(ybase))
  n <- length(ybase)
  thetau <- covariate + rnorm(n, mean = sd(ybase) * tau_sds, sd = sd(ybase) / 2)
  return(thetau)
}

#' @describeIn Tau_Functions A basic function that specifies a tau_sds size effect that varies, randomly, by level of covariate. We envision this to be used with covariate that have relatively few levels into which the outcome values can be grouped.
#' @export
tau_norm_covariate_levels <- function(ybase, tau_sds, covariate) {
  stopifnot(length(covariate) == length(ybase))
  stopifnot(length(unique(covariate)) < length(unique(ybase)) / 10)
  n <- length(ybase)
  thetau <- unsplit(lapply(split(ybase, covariate), function(theys) {
    rnorm(length(theys),
      mean = (runif(1, min(theys), max(theys)) + sd(theys) * tau_sds),
      sd = sd(theys) / 2
    )
  }), covariate)
  return(thetau)
}
