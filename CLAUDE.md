# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

TreeTestSim is an R package for studying hierarchical hypothesis testing procedures when hypotheses are organized in tree structures. It supports both simulation-based assessment (using abstract beta-distributed p-values) and concrete block-randomized experiments.

## Related Projects

This package is part of a trio of related projects:

- **~/repos/manytests-paper**: The academic paper motivating and presenting the methods
- **~/repos/manytestsr**: The user-facing R package implementing the methods for practitioners
- **~/repos/TreeTestSim** (this repo): Simulation infrastructure for assessing method performance

TreeTestSim provides the simulation evidence (FWER control, power comparisons) that validates the methods described in the paper and implemented in manytestsr.

## Coding Preferences

See `CLAUDE_CODING.md` for detailed guidance on:
- Reading relevant files broadly before making changes
- Preferring boring code over clever code (but using vectorization over for-loops)
- Commenting **why** not **what**
- Writing tests that verify statistical principles before refactoring
- Pausing for review at checkpoints (after tests, after implementation, at design decisions)

## Development Commands

```bash
# Interactive R session with package loaded (enables help on package functions)
make interactive

# Run tests
make test

# Build documentation
make document

# Full package check
make check

# Install dependencies
make dependencies

# Build package
make build

# Clean generated files
make clean
```

Single test file: `R -e "devtools::test(filter = 'test_general')"`

## Architecture

The package has two parallel testing frameworks:

### 1. Abstract Tree Simulations (`R/simulate_test.R`, `R/simulate_many_runs.R`)

Uses beta-distributed p-values to simulate hierarchical testing without real data:

- `generate_tree_DT()` creates k-ary trees with known null/non-null status at leaves, propagated bottom-up
- `simulate_test_DT()` runs a single hierarchical test with configurable alpha allocation methods:
  - `"fixed"`: constant alpha at every node
  - `"spending"`: alpha reduced at each branch
  - `"investing"`: rejections earn bonus alpha
  - `"fixed_k_adj"` / `"adaptive_k_adj"`: conservative adjustments
  - `"adaptive_power"`: depth-indexed schedule via `manytestsr::compute_adaptive_alphas`
- `simulate_many_runs_DT()` averages results across simulation replicates

Key design: p-values are constrained to be monotonically non-decreasing down the tree (child p-values >= parent p-values).

### 2. Concrete Block-Randomized Experiments (`R/synthetic_data.R`, `R/padj_sim_fns.R`)

Works with actual individual-level data via `manytestsr::find_blocks`:

- `generate_synthetic_experiment()` creates individual/block data.tables for a k-ary tree experiment
- `padj_test_fn()` is the main simulation driver—creates treatment effects, re-randomizes, tests, summarizes errors
- Two testing modes controlled by `p_adj_method`:
  - `"split"`: top-down via `manytestsr::find_blocks` with splitting functions
  - `"fdr"/"holm"/etc.: bottom-up via `adjust_block_tests()` using `p.adjust`

### Supporting Modules

- `R/p_adjustment.R`: Local p-value adjustment functions (`local_simes`, `local_hommel_all_ps`, `local_bh_all_ps`, `local_unadj_all_ps`)
- `R/find_blocks_simulations.R`: Treatment effect creation (`create_effects()`, `tau_*` functions)
- `R/summarize_errors_discoveries.R`: Error rate calculations (`calc_errs`, `calc_errs_new`)

## Key Dependencies

- `manytestsr` (from GitHub: `bowers-illinois-edu/manytestsr`): provides `find_blocks`, `alpha_adaptive`, `compute_adaptive_alphas`, splitting functions
- `data.table`: all data manipulation uses data.table syntax
- `hommel`: Hommel p-value adjustment

## Testing Notes

Tests use `testthat` edition 3. Time-consuming FWER simulations are in `tests.notforCRAN.R` and skipped by default. The main test file `test_general.R` checks all combinations of `alpha_method`, `final_adj_method`, `local_adj_method`, and `adj_effN` parameters.

## Data Structures

Tree nodes are stored in data.tables with columns: `node`, `level` (0 = root), `parent`, `nonnull`. Simulation results include `p_val`, `alpha_alloc`, `p_sim`, and various adjusted p-value columns.

Block-level data uses columns: `bF`/`blockF` (block factor), `lvls_fac` (encodes tree ancestry as "parent.child.grandchild"), `nonnull`, `hwt` (harmonic mean weight).
