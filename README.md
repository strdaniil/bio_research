# Simulated Genome Evolution via Phylogenetic Tree Traversal

This project simulates the evolution of genomes along a phylogenetic tree using gene gain, loss, and inversion events. It then compares the synteny block distributions between simulated and real genome data to assess how well the simulation mimics biological patterns.

## Overview

- **Input**
  - A rooted Newick-format tree (e.g., `atgc.iq.r.tre`)
  - A `.cc.csv` file from the ATGC database containing gene orders and COG annotations
- **Output**
  - Simulated genomes along each branch of the tree
  - Synteny block comparisons:
    - Simulated vs simulated genomes
    - Real vs real genomes
    - Real vs corresponding simulated genome

## What the Script Does

1. **Loads real genome data** from an ATGC `.cc.csv` file and reconstructs gene orders per genome.
2. **Selects a root genome** as the ancestor for simulation.
3. **Simulates genome evolution** along a rooted tree by modeling:
   - Gene gain (insertion of novel genes)
   - Gene loss (deletion of existing genes)
   - Inversions (optional, disabled by default)
4. **Traverses the phylogenetic tree** and generates descendant genomes for each node using the specified mutation rates.
5. **Computes synteny blocks** — conserved contiguous gene segments — between selected genome pairs.
6. **Prints block length distributions** to compare the structure of real and simulated genome alignments.

