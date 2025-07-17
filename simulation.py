from Bio import Phylo
import pandas as pd
import numpy as np
import random
from itertools import combinations
from collections import Counter
import matplotlib.pyplot as plt
from synteny_tools import findSyntenyReal

def evolve_genome_branch(genome, branch_length, gain_rate, loss_rate, inv_rate, gain_genes):
    genome = genome.copy()
    gains = np.random.poisson(gain_rate * branch_length)
    losses = np.random.poisson(loss_rate * branch_length)
    inversions = np.random.poisson(inv_rate * branch_length)

    for _ in range(gains):
        new_gene = random.choice(gain_genes)
        idx = random.randrange(len(genome)+1)
        genome.insert(idx, new_gene)

    for _ in range(losses):
        if genome:
            idx = random.randrange(len(genome))
            genome[idx] = -1

    for _ in range(inversions):
        if len(genome) >= 2:
            i = random.randrange(len(genome))
            j = random.randrange(i, len(genome))
            genome[i:j+1] = list(reversed(genome[i:j+1]))
    return genome

def median_root_to_leaf_lengths(tree):
    return np.median([tree.distance(tree.root, leaf) for leaf in tree.get_terminals()])

def evolve_genome(tree, root_genome, per_gene_gain_rate, per_gene_loss_rate, per_gene_inv_rate, gain_genes):
    genomes = {}
    num_genes = len(root_genome)
    median_path_length = median_root_to_leaf_lengths(tree)

    branch_gain_rate = (num_genes * per_gene_gain_rate) / median_path_length
    branch_loss_rate = (num_genes * per_gene_loss_rate) / median_path_length
    branch_inv_rate = (num_genes * per_gene_inv_rate) / median_path_length

    genomes[tree.root] = root_genome.copy()
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            parent_genome = genomes[parent]
            child_genome = evolve_genome_branch(parent_genome, child.branch_length,
                                                branch_gain_rate, branch_loss_rate, branch_inv_rate, gain_genes)
            genomes[child] = child_genome
    return genomes

def synteny_blocks(genome1, genome2):
    pos_in_B = {gene: idx for idx, gene in enumerate(genome2)}
    blocks = []
    i = 0
    while i < len(genome1):
        gene = genome1[i]
        if gene in pos_in_B:
            j = pos_in_B[gene]
            block_len = 1
            while i + block_len < len(genome1):
                next_gene = genome1[i + block_len]
                if next_gene in pos_in_B and pos_in_B[next_gene] == j + block_len:
                    block_len += 1
                else:
                    break
            blocks.append(block_len)
            i += block_len
        else:
            i += 1
    return blocks

def get_real_genomes_from_cc(path):
    df = pd.read_csv(path, names=[
        "gene_ID", "genome_ID", "protein_ID", "protein_length",
        "atgc_cog_footprint", "atgc_cog_footprint_length",
        "atgc_cog_ID", "protein_cluster_ID", "match_class"
    ])
    df = df.dropna(subset=["atgc_cog_ID"])
    return df.groupby("genome_ID")["atgc_cog_ID"].apply(list).to_dict()

# Config
DATASET = "ATGC0070"
TREE_FILE = f"{DATASET}/atgc.iq.r.tre"
CC_FILE = f"{DATASET}/atgc.cc.csv"

# Load tree
tree = Phylo.read(TREE_FILE, "newick")




# Load root genome (just use first genome from file)
real_genomes = get_real_genomes_from_cc(CC_FILE)
with open("ATGC0070/atgc.info.tab", "r") as f:
    first_line = f.readline().strip()
    first_value = first_line.split()[0]
root_genome = real_genomes[first_value]

# for randomly inserting genes that exist in the tree not something else; IS THIS RIGHT?
all_cls_ids = list({gene for genes in real_genomes.values() for gene in genes})

# Simulate all genomes
simulated_genomes = evolve_genome(tree, root_genome,
                                   per_gene_gain_rate=0.1,
                                   per_gene_loss_rate=0.1,
                                   per_gene_inv_rate=0.0,
                                   gain_genes = all_cls_ids)

# Collect synteny block lengths for all leaf pairs
leaves = tree.get_terminals()
all_block_lengths = []

for leaf1, leaf2 in combinations(leaves, 2):
    genome1 = simulated_genomes[leaf1]
    genome2 = simulated_genomes[leaf2]
    blocks = findSyntenyReal(genome1, genome2)
    all_block_lengths.extend([abs(b[4]) for b in blocks if len(b) > 4])

# Summarize block length distribution
# print(len(all_block_lengths))

def run_simulation(tree, root_genome, gain, loss, inv_rate, gain_genes):
    genomes = evolve_genome(tree, root_genome, gain, loss, inv_rate, gain_genes)
    leaves = tree.get_terminals()
    all_block_lengths = []

    for leaf1, leaf2 in combinations(leaves, 2):
        g1 = genomes[leaf1].copy()
        g2 = genomes[leaf2].copy()
        blocks = findSyntenyReal(g1, g2)
        all_block_lengths.extend([abs(b[4]) for b in blocks if len(b) > 4])

    return all_block_lengths



# plt.hist(all_block_lengths, bins=range(1, max(all_block_lengths)+2), density=True, edgecolor='black')
# plt.xlabel("Synteny Block Length")
# plt.ylabel("Normalized Frequency")
# plt.title("Distribution of Synteny Blocks (All Leaf Pairs)")
# plt.tight_layout()
# plt.show()

# real_lengths = []
# for leaf1, leaf2 in combinations(leaves, 2):
#     g1 = real_genomes[leaf1.name]
#     g2 = real_genomes[leaf2.name]
#     blocks = findSyntenyReal(g1.copy(), g2.copy())
#     block_lengths = [block[4] for block in blocks]
#     real_lengths.extend(block_lengths)
# with open("real_lengths.txt", "w") as f:
#     f.write(" ".join(map(str, real_lengths)))
