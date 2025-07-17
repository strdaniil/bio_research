## Overview

Given a phylogenetic tree and gene cluster assignments (`cls.` IDs) from a real dataset, we:

1. **Simulate genome evolution** from a root genome across a phylogenetic tree.
2. **Compute synteny blocks** for all pairs of leaf genomes.
3. **Compare distributions** of simulated vs. real synteny block lengths using Wasserstein distance.
4. **Optimize parameters** using Bayesian optimization to minimize this distance.

---

## Key Files

### 🔹 `simulation.py`

- Core simulation engine for genome evolution.
- Implements gain, loss, and inversion events along tree branches.
- Computes synteny blocks across all simulated leaf genome pairs.
- Provides `run_simulation(...)` to return a list of block lengths.

### 🔹 `sim_real_comparison.py`

- Loads the real phylogenetic tree and genome data.
- Loads true synteny block lengths (`real_lengths.txt`).
- Uses `skopt` to optimize gain/loss parameters via Bayesian optimization.
- Minimizes the Wasserstein distance between simulated and real block length distributions.

### 🔹 `calc_real_data.py`

- Parses real genome cluster assignments (`atgc.cc.csv`).
- Computes real synteny block lengths between leaf genomes using `findSyntenyReal`.
- Used to precompute and save the `real_lengths.txt` file.
