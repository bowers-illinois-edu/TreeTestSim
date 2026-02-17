# TreeTestSim 0.0.0.9303

## New features

* `simulate_test_DT()` gains `alpha_method = "adaptive_power_pruned"`.
  Like `"adaptive_power"`, this pre-computes a depth-indexed alpha schedule
  from estimated power decay. After each tree level, it counts how many
  nodes survived (rejected) and recomputes the schedule for deeper levels
  using the pruned node count instead of the full tree's `k^l`. When
  branches die early, surviving branches get more alpha budget. The FWER
  guarantee follows from the same telescoping-sum argument applied to the
  pruned subtree.

* `simulate_many_runs_DT()` passes `alpha_method = "adaptive_power_pruned"`
  through to `simulate_test_DT()` without modification.

# TreeTestSim 0.0.0.9302

## New features

* `padj_test_fn()` gains a `non_null_blocks` parameter. When provided (as the
  name of a logical column in `bdat`), the same blocks are treated as non-null
  across scenarios rather than being randomly assigned via `prop_blocks_0`.

## Breaking changes

* Removed the `tau` parameter from `simulate_test_DT()`.
  `manytestsr::compute_adaptive_alphas()` now uses the error load
  framework instead of `tau` — when the total error load is at most 1,
  nominal alpha is used everywhere (natural gating); when it exceeds 1,
  alpha is adjusted automatically. No tuning parameter is needed.

# TreeTestSim 0.0.0.9300

## Breaking changes

* `simulate_test_DT()` and `simulate_many_runs_DT()` now parameterize
  non-null p-value generation via `effect_size` (Cohen's d) instead of
  `beta_base`. The new parameter has a direct scientific interpretation: it
  is the standardized mean difference at each non-null leaf. Together with
  `N_total` and the tree structure, power at every level is derived
  automatically.

* Removed `beta_base`, `effN`, `adj_effN`, and `delta_hat` parameters from
  `simulate_test_DT()`. Removed `beta_base` and `adj_effN` from
  `simulate_many_runs_DT()`. The old `adj_effN` logical is replaced by
  `power_decay` (same default behavior: `TRUE`).

* Removed internal helper `derive_delta_hat()`. Its role is replaced by
  `effect_size_to_beta()`, which converts Cohen's d directly to a Beta shape
  parameter without the intermediate round-trip.

## New features

* `simulate_test_DT()` now returns `root_power` in `sim_res` — the
  theoretical power at the root level computed from `effect_size` and
  `N_total`.

* New internal helper `effect_size_to_beta(effect_size, N_level, alpha)`
  converts Cohen's d at a given sample size to the Beta(a, 1) shape
  parameter used for p-value generation.

# TreeTestSim 0.0.0.9201

## Bug fixes

* Fixed test helper `err_testing_fn()` using bare `filter()` which resolved to
  `stats::filter` instead of `dplyr::filter` during non-interactive test runs.
  Replaced with data.table subsetting to remove the implicit dplyr dependency.

# TreeTestSim 0.0.0.9200

## New features

* `simulate_test_DT()` gains `alpha_method = "adaptive_power"`, which uses
  `manytestsr::compute_adaptive_alphas()` to precompute a depth-indexed alpha
  schedule based on estimated power decay through the tree. New parameters
  `delta_hat` and `tau` control the schedule; when `delta_hat = NULL` it is
  derived automatically from `beta_base` and `N_total`.

* `simulate_test_DT()` now always computes **both** bottom-up hommel and
  bottom-up BH adjustments on leaf p-values, returning six new columns:
  `bu_hommel_false_error`, `bu_hommel_true_disc`, `bu_hommel_power`,
  `bu_bh_false_error`, `bu_bh_true_disc`, `bu_bh_power`. The existing
  `bottom_up_*` columns are preserved for backward compatibility.

* New function `generate_synthetic_experiment()` creates a synthetic
  block-randomized experiment organized as a complete k-ary tree. Returns
  individual-level and block-level data.tables with a `lvls_fac` factor
  suitable for `manytestsr::splitSpecifiedFactorMulti()`.

* New function `run_synthetic_comparison()` runs all six testing methods
  (four top-down, two bottom-up) on synthetic experiment data via
  `padj_test_fn()`. Returns a tidy data.table with one row per method.

* Internal helper `derive_delta_hat()` converts `beta_base` and `N_total` to
  the equivalent `delta_hat` for the normal power formula by matching
  root-level power.

## Dependencies

* Now requires `manytestsr (>= 0.0.4)` for `compute_adaptive_alphas()` and
  `alpha_adaptive()`.

# TreeTestSim 0.0.0.9100

* Initial development release with `simulate_test_DT()`,
  `simulate_many_runs_DT()`, `generate_tree_DT()`, local p-value adjustment
  functions, `padj_test_fn()`, `create_effects()`, tau functions, and error
  calculation utilities.
