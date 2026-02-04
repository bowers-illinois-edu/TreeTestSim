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
