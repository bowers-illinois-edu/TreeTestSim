## Make some testing data for use in the tests in this directory
message("test")
library(randomizr)
library(data.table)
library(here)
library(dtplyr)
library(dplyr)
# library(TreeTestsSim)
load_all("~/repos/TreeTestsSim")
library(conflicted)
conflicts_prefer(dplyr::filter)

setDTthreads(1)
set.seed(12345)

## For now, create a very simple canceling out arrangement and simple design setup
idat <- data.table(i = 1:1000, b = rep(c(1:10), length = 1000))
bdat <- data.table(b = 1:10)
setkey(bdat, b)
setkey(idat, b)
## Make vb1, vb2, vb3 all nest within each other like states, counties, districts (with the individual blocks nested within all three)
## 2 levels at highest level
bdat[, vb1 := rep(c(1, 2), length.out = .N)]
## 2 levels within each state
bdat[, vb2 := rep(1:2, length.out = .N), by = vb1]
bdat[, vbnest := interaction(vb1, vb2, lex.order = TRUE, drop = TRUE)]
ftable(bdat$vb1, bdat$vb2)
## Now merge the block data onto the individual level data
idat <- bdat[idat]
set.seed(12345)
idat[, vi1 := round(runif(.N, min = 0, max = 10)), by = b]
idat[, y0 := vi1 + vb1 + rnorm(100), by = b]
idat[, tau := ifelse(b <= 5, -5 * sd(y0), 5 * sd(y0))]
idat[, tauhomog := sd(y0) * 5]
idat[, taunormb := rnorm(.N, mean = (sd(y0) / b), sd = sd(y0) / 2), by = b]
idat[, y1 := y0 + tau]
idat[, y1homog := y0 + tauhomog]
idat[, y1normb := y0 + taunormb]
stopifnot(all.equal(mean(idat$y1 - idat$y0), 0))
idat[, Z := complete_ra(N = 100), by = b]
testra <- idat[, sum(Z), by = b]
stopifnot(all(testra$V1 == 50))
idat[, Y := Z * y1 + (1 - Z) * y0]
## Now ensure that the no effects case is really no effects
idat[, y0null := resid(lm(y0 ~ Z, data = .SD)), by = b]
idat[, y1null := y0null]
idat[, Ynull := Z * y1null + (1 - Z) * y0null]
stopifnot(all.equal(idat$Ynull, idat$y0null))
idat[, Yhomog := Z * y1homog + (1 - Z) * y0]
idat[, Ynormb := Z * y1normb + (1 - Z) * y0]
idat$ZF <- factor(idat$Z)
idat$bF <- as.character(idat$b)
idat$gF <- interaction(idat$bF, idat$ZF)
## Created aligned versions of the  variables
idat[, Ymd := Y - mean(Y), by = b]
idat[, Zmd := Z - mean(Z), by = b]

## Make some blocks have nothing going on so the splitting doesn't go all the way to the end
set.seed(12345)
idat[bF %in% c(7, 8), Ynormb := Ynull + rnorm(.N)]

## Observed Effects by Block
idat[, list(
  canceling = mean(Y[Z == 1] - Y[Z == 0]),
  null = mean(Ynull[Z == 1] - Ynull[Z == 0]),
  homog = mean(Yhomog[Z == 1] - Yhomog[Z == 0]),
  normb = mean(Ynormb[Z == 1] - Ynormb[Z == 0])
), by = b]

##  Make a smaller dataset
b1n6 <- droplevels(idat[b %in% c(1, 6), ])
setorder(b1n6, b)
b1n6
b1 <- droplevels(b1n6[b == 1, ])

bdat_tmp <- idat[, .(nb = .N, nt = sum(Z), nc = sum(1 - Z), pb = mean(Z)), by = b]
bdat_tmp[, hwt := (2 * (nc * nt) / (nc + nt))]
setkey(bdat_tmp, b)
bdat <- bdat_tmp[bdat]
bdat$bF <- factor(bdat$b)

setkey(idat, bF)
setkey(bdat, bF)

#### Now make data with unequal sized blocks, unequal taus, and more blocks
bdat2 <- bdat[, .(nb = round(seq(1, 20, length = 20), 2), bF2 = as.character(1:20))]
setkey(bdat2, "bF2")
idat[, b2 := {
  sel <- rep(c(0, 1), .N / 2)
  sel * b + (1 - sel) * (10 + b)
}, by = bF]
idat[, bF2 := factor(b2)]
setkey(idat, "bF2")
## Make a dataset with unequal sized blocks
idat2 <- bdat2[idat]
idat3 <- idat2[, sample(.N, size = unique((.N * nb) / 10), replace = TRUE), by = bF2]
## Now giving it the same blocking names as idat for ease with testing
idat3$bF <- idat3$bF2
setkey(idat3, bF)
idat3 <- idat3[, b := as.numeric(as.character(bF2))]
idat3[, c("v1", "v2", "v3", "v4") := lapply(1:4, function(i) {
  rnorm(.N, mean = .N / i, sd = 1 / i)
}), by = bF2]
idat3[, y0 := v1 + rnorm(.N), by = bF2]
idat3[, tau := ifelse(b <= 10, -5 * sd(y0), 5 * sd(y0))]
idat3[, tauhomog := 5 * sd(y0)]
## Taunorm_inc increases the effect with the block number and also size of the block.
## We should be increasingly likely to detect effects for blocks with higher numbers.
idat3[, taunorm_inc := rnorm(.N, mean = b / 10 * sd(y0), sd = .5), by = bF2]
idat3[, taunorm_inc := ifelse(taunorm_inc < 0, 0, taunorm_inc)]
## taunorm_dec has effect size decreasing in block size
idat3[, taunorm_dec := rnorm(.N, mean = 4 / b * sd(y0), sd = .5), by = bF2]
idat3[, taunorm_dec := ifelse(taunorm_dec < 0, 0, taunorm_dec)]
idat3[, sd(y0), by = b]
## tauv2 has no systematic relationship between block number or block size and effect size

idat3[, tauv2 := ifelse(v2 < quantile(v2, sample(c(.9, .1), size = 1)),
  sd(y0) * (taunorm_inc) + sd(y0) * runif(.N, min = .25, max = 2), 0
), by = bF2]

idat3[, lapply(.SD, mean), .SDcols = grep("tau|nb|y0", names(idat3), value = TRUE), by = b]
idat3[, y1 := y0 + tau]
idat3[, y1null := y0]
idat3[, y1homog := y0 + tauhomog]
set.seed(1235) ## set y1norm_inc=y0 for 1/4 of the blocks
nullbF <- sample(unique(idat3$bF), size = length(unique(idat3$bF)) / 4)
idat3[, y1norm_inc := ifelse(bF %in% nullbF, y0, y0 + taunorm_inc)]
idat3[, y1norm_dec := ifelse(bF %in% nullbF, y0, y0 + taunorm_dec)]
idat3[, y1tauv2 := y0 + tauv2]
idat3[, nb := .N, by = bF2]
set.seed(1235)
tmp <- idat3[, .(nb = .N), by = bF2]
mb <- round(tmp$nb * runif(20, min = .2, max = .8))
idat3[, Z := block_ra(blocks = bF2, block_m = mb)]
idat3[, Y := Z * y1 + (1 - Z) * y0]
idat3[, Ynull := Z * y1null + (1 - Z) * y0]
idat3[, Yhomog := Z * y1homog + (1 - Z) * y0]
idat3[, Ynorm_inc := Z * y1norm_inc + (1 - Z) * y0]
idat3[, Ynorm_dec := Z * y1norm_dec + (1 - Z) * y0]
idat3[, Ytauv2 := Z * y1tauv2 + (1 - Z) * y0]
idat3$ZF <- factor(idat3$Z)
idat3[, .(.N, mean(y1tauv2 - y0)), by = bF]
bdat3 <- idat3[, .(
  nb = .N, nt = sum(Z), nc = sum(1 - Z), pb = mean(Z),
  v1 = mean(v1), v2 = mean(v2), v3 = mean(v3), v4 = mean(v4),
  ate_tau = mean(y1 - y0),
  ate_null = mean(y1null - y0),
  ate_homog = mean(y1homog - y0),
  ate_norm_inc = mean(y1norm_inc - y0),
  ate_norm_dec = mean(y1norm_dec - y0),
  ate_tauv2 = mean(y1tauv2 - y0)
), by = bF]
bdat3[, hwt := (2 * (nc * nt) / (nc + nt))]
bdat3[, bF := as.character(bF)]
setkey(idat, bF)
setkey(bdat, bF)
setkey(idat3, bF)
setkey(bdat3, bF)

##############################
#### Now generate complete k-ary tree with known null and non-null blocks
### (blocks with known no-effects or known to have non-zero effects)

## We think this kind of structure will cause trouble to unadjusted approaches.
## A tree with 4 nodes per levels and 3 levels after the root where 50% of the
## leaves are non-null the tree is defined at the level of the block We are not
## going to pass this vector below when we assess weak FWER control

bdt <- generate_tree_DT(max_level = 3, k = 4, t = .5)
bdt[, leaf := as.numeric(level == max(level))]
table(bdt$leaf)
stopifnot(nrow(bdt) == sum(4^seq(0, 3)))

### This next converts the node level data into a block level dataset with a
### factor variable with "." indicating relationships between nodes (each "."
### is a level in the tree)

# build a named vector so parent_lookup["22"] == 6, etc.
parent_lookup <- setNames(bdt$parent, bdt$node)
## Compare to the mathematical way to calculate the parent
expect_equal(parent_lookup[22][[1]], floor((22 - 2) / 4) + 1)
expect_equal(parent_lookup[22][[1]], 6)

# for a given node, return the vector of its ancestors+itself, *excluding* the
# NA‐parent root. "node" is an integer that counts from 1 at the root upwards
# in a structured order (see the math above about how to calculate the node
# numbers of parents, similar calculations are used for children in
# generate_tree_DT)

get_path <- function(node) {
  path <- integer(0)
  current <- node
  # walk up until you hit an NA parent
  while (!is.na(parent_lookup[as.character(current)])) {
    path <- c(path, current) # record this node
    current <- parent_lookup[as.character(current)] # jump to parent
  }
  # reverse so it reads root→…→node
  rev(path)
}

#  apply get_path to every node, paste the names with ".", and coerce to factor since
#  find_blocks and splitSpecifiedFactorMulti requires such a factor to split
#  the data according to the pre-specified structure

bdt[, lvls_fac := factor(sapply(node, function(n) {
  paste(get_path(n), collapse = ".")
}))]
bdt[node == 6, ]
bdt[parent == 6, ]

## We are going to focus on the leaves which are the individual blocks.
## That is the leaf nodes are the block-level dataset.
## This is the kind of data that an experimenter would have access to directly.

bdt1 <- droplevels(bdt[leaf == 1, ])

## Assess the splitSpecifiedFactorMulti splitting function
## The splitSpecifiedFactorMulti function should recreate this tree.
bdt1[, split1 := splitSpecifiedFactorMulti(bid = node, x = lvls_fac)]
bdt1[, split2 := splitSpecifiedFactorMulti(node, lvls_fac), by = split1]
bdt1[, split3 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2)]

with(bdt1, table(split1, exclude = c()))
with(bdt1, ftable(split1 = split1, split2 = split2, split3 = split3))

## Everybody within the same parent should be in the same first level split
test_split_1 <- bdt1 %>%
  group_by(parent) %>%
  summarize(val = length(unique(split1)))
stopifnot(all(test_split_1$val == 1))

## Split2 is constant within split1 --- so there are blocks with split2==1 within each value of split1.
## i.e. first group within split1=1, split1=2, etc..
test_split_2 <- bdt1 %>%
  group_by(split1) %>%
  summarize(val = length(unique(split2)))
stopifnot(all(test_split_2$val == 4))
stopifnot(nrow(test_split_2) == 4)

bdt1 %>% filter(split1 == 1)
bdt1 %>% filter(split2 == 1)

## This next returns all 0s when splitting is no longer possible. So, consider
## watching out for this in the find_blocks code

bdt1[, split4 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3)]
bdt1[, split5 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3, split4)]

expect_equal(as.numeric(as.character(unique(bdt1$split4))), 0)
expect_equal(as.numeric(as.character(unique(bdt1$split5))), 0)

bdt1[, c("split1", "split2", "split3", "split4", "split5") := NULL]

#### Done assessing the splitSpecifiedFactorMulti function on the data created by generate_tree_DT

head(bdt1)

## Make individual level data with nb=50
## Idea here is to make something with high power within each block
bdt1[, nb := 50]
bdt1[, bary0 := 0]
bdt1[, bF := as.character(node)]

idt <- data.table(b = rep(bdt1$node, bdt1$nb))
idt[, bF := as.character(b)]
idt[, id := seq_len(nrow(idt))]
idt[, bary0 := rep(bdt1$bary0, bdt1$nb)]

## this next "by=b" is not necessary since it is the same y0 for all blocks.
## But it will not be that case in other situations

idt[, y0 := rnorm(.N, bary0, sd = 1), by = b]
idt <- merge(idt, bdt1[, .(bF, nonnull, lvls_fac)], by = "bF")

## Here we ignore the nonnull variable on bdt1 for the sake of creating a
## situation where there are no effects in any block or for any person.

idt[, y1 := create_effects(
  idat = idt, ybase = "y0",
  blockid = "bF", tau_fn = tau_null, tau_size = 0, prop_blocks_0 = 1, non_null_blocks = NULL
)]

## Treatment is completely randomized in each block
idt[, trt := sample(rep(c(0, 1), .N / 2)), by = bF]
## Reveal the observed Y as a function of the potential outcomes
idt[, Y := trt * y1 + (1 - trt) * y0]

## Some procedures need treatment to be a factor (esp those using the coin
## package)
idt[, trtF := factor(trt)]

## test using some non-null effect blocks so that the algorithm descends more into the tree
## This creates effects block-by-block.
idt[, y1_half_tau1 := create_effects(
  idat = idt, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 1, prop_blocks_0 = .5,
  non_null_blocks = "nonnull"
)]
## Check that the blocks with no effects actually have no effects
test_null_effects <- idt[nonnull == FALSE, sum(y1_half_tau1 - y0)]
expect_equal(unique(test_null_effects), 0)

test_not_null <- idt[nonnull == TRUE, sum(y1_half_tau1 - y0)]
expect_gt(abs(unique(test_not_null)), 0)
## Make an observed outcome
idt[, Y_half_tau1 := trt * y1_half_tau1 + (1 - trt) * y0]
