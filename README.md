Gene Gain–Loss & Synteny Simulation Toolkit

A lightweight Python workflow for simulating genome evolution on a phylogenetic tree and comparing simulated vs. observed synteny‐block patterns in microbial genomes (ATGC dataset).
1. Project Goals

    Simulate gene gains, losses, and inversions along every branch of a user-supplied Newick tree.

    Reconstruct extant (leaf) genomes from the simulation.

    Measure conserved gene-order (“synteny”) blocks between any pair of genomes.

    Benchmark simulated block-length distributions against real ATGC genomes.

2. Repository Layout

project-root/
├── ATGC0001/                # Example dataset from NCBI ATGC
│   ├── atgc.cc.csv          # Gene-content & order table
│   └── atgc.iq.r.tre        # Maximum-likelihood tree in Newick
├── simulate.py              # Main script (code shown below)
├── requirements.txt         # Exact Python package versions
└── README.md                # This file

Feel free to add multiple ATGCxxxx folders; pass the folder name to the script with --dataset.
3. Quick Start
3.1 Clone & install dependencies

git clone https://github.com/your-username/genome-synteny-sim.git
cd genome-synteny-sim

# create a fresh virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# install everything
pip install -r requirements.txt

<details><summary>requirements.txt (exact versions)</summary>

biopython==1.83
ete3==3.1.3
numpy==1.26.4
pandas==2.2.2

</details>
3.2 Run a simulation

python simulate.py \
    --dataset ATGC0001 \
    --gain 0.10 \
    --loss 0.10 \
    --inv  0.00

The script prints three sections:

    Simulated synteny blocks for two random leaves

    Real synteny blocks for the matching ATGC genomes

    Real vs. simulated block distribution on the same leaf

4. Script Overview (simulate.py)
Function	Purpose
traverse(clade)	Utility to print the tree topology with branch lengths.
evolve_genome_branch()	Applies Poisson-sampled gain/loss/inversion events to an input genome across one branch.
median_root_to_leaf_lengths()	Scaling factor so that per-gene rates map to branch-specific rates.
evolve_genome()	Recursively builds simulated genomes for every node in the tree.
synteny_blocks(g1, g2)	Returns a list of contiguous block lengths shared between two gene orders.
get_real_genomes_from_cc()	Loads ordered gene lists for all genomes in a CSV exported by ATGC.

Key command-line arguments:
Flag	Meaning	Default
--dataset	Folder name containing atgc.cc.csv & atgc.iq.r.tre	required
--gain	Per-gene gain rate	0.10
--loss	Per-gene loss rate	0.10
--inv	Per-gene inversion rate	0.00
5. Data Preparation

    Download an ATGC cluster from NCBI (or your own dataset) so that each folder contains

        atgc.cc.csv — annotated gene table with ordered atgc_cog_ID

        atgc.iq.r.tre — phylogenetic tree whose tip names match genome_ID in the CSV

    Place the folder inside project-root/.

    Reference it with --dataset.

6. Customising Simulations

    Rates Pass different --gain, --loss, --inv values to explore scenarios.

    Alternative root genome Edit the code block where the first genome column is chosen; pick any column or supply a curated root list.

    Different distance metrics Swap out synteny_blocks for another genome-comparison method if desired.

7. Example Output

Simulated synteny blocks between GCF_00000000.1 and GCF_00000001.1
Block length distribution: {3: 5, 2: 12, 1: 69}

Real synteny blocks between GCF_00000000.1 and GCF_00000001.1
Block length distribution: {10: 3, 4: 1, 12: 1, 69: 1}

Comparing Real Genome GCF_00000000.1 with Simulated GCF_00000000.1:
Synteny block lengths: [10, 3, 4, 12, 69]
Block length distribution: {3: 1, 4: 1, 69: 1, 10: 1, 12: 1}

Interpretation: Your simulation roughly reproduces the long conserved regions (~69 genes) but over-fragments shorter blocks—consider lowering the loss rate.
8. Troubleshooting
Symptom	Likely Cause	Fix
KeyError on genome name	Tip labels in tree ≠ genome_ID in CSV	Rename tips or update CSV.
Zero-length branch warnings	Tree contains polytomies or zero-length branches	Reroot/re-branch with IQ-TREE or RAxML.
All simulated genomes identical	Gain/loss/inv rates too low or branch_length==0	Inspect branch lengths; scale rates up.
9. References

    Wolf, Koonin et al. ATGC Database of Genomic Clusters. Nucleic Acids Res 2012.

    Biopython (Cock, P.J.A. et al.)

    Huerta-Cepas J., Dopazo J. ETE3: reconstruction, analysis and visualization of phylogenomic data sets.

10. License

MIT License – see LICENSE file for details.
