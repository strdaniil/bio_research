from Bio import Phylo
import pandas as pd
import numpy as np
import random
from itertools import combinations
from collections import Counter
import matplotlib.pyplot as plt
from synteny_tools import findSyntenyReal
import random

def sample_truncated_powerlaw(L0, alpha=2.0):
    k_max = max(1, L0 // 2)
    ks = np.arange(1, k_max + 1)
    weights = ks ** (-alpha)
    weights /= weights.sum()
    return int(np.random.choice(ks, p=weights))


def evolve_genome_branch(genome, branch_length, gain_rate, loss_rate, inv_rate, trans_rate, gain_exp, loss_exp, inv_exp, trans_exp, next_gene_id_holder, L0_ancestral):
    genome = genome.copy()

    gains = np.random.poisson(gain_rate * branch_length)
    losses = np.random.poisson(loss_rate * branch_length)
    inversions = np.random.poisson(inv_rate * branch_length)
    translocations = np.random.poisson(trans_rate * branch_length)

    events = (["G"] * gains +
              ["L"] * losses +
              ["I"] * inversions +
              ["T"] * translocations)
    random.shuffle(events)

    for event in events:
        Lc = len(genome)
        loc_point = random.randrange(Lc)

        if event == "G":
            k = sample_truncated_powerlaw(L0_ancestral, gain_exp)
        elif event == "L":
            k = sample_truncated_powerlaw(L0_ancestral, loss_exp)
        elif event == "I":
            k = sample_truncated_powerlaw(L0_ancestral, inv_exp)
        else:
            k = sample_truncated_powerlaw(L0_ancestral, trans_exp)

        k = max(1, min(k, Lc))

        if event == "G":
            insert_pos = (loc_point + 1) % (Lc + 1)
            block = []
            for _ in range(k):
                new_gene = f"cog.{next_gene_id_holder[0]}"
                next_gene_id_holder[0] += 1
                block.append(new_gene)
            genome[insert_pos:insert_pos] = block

        elif event == "L":
            for t in range(k):
                pos = (loc_point + t) % len(genome)
                genome[pos] = -1

        elif event == "I":
            if k > 1 and len(genome) >= 2:
                idxs = [(loc_point + t) % len(genome) for t in range(k)]
                vals = [genome[p] for p in idxs][::-1]
                for p, v in zip(idxs, vals):
                    genome[p] = v

        elif event == "T":
            if len(genome) >= 2 and 0 < k < len(genome):
                Lc_now = len(genome)
                idxs = [(loc_point + t) % Lc_now for t in range(k)]
                block_vals = [genome[p] for p in idxs]
                j = random.randrange(Lc_now)
                if j in idxs:
                    continue
                removed = set(idxs)
                remaining = [g for t, g in enumerate(genome) if t not in removed]
                r_before = sum(1 for t in idxs if t < j)
                j_new = j - r_before
                remaining[j_new + 1:j_new + 1] = block_vals
                genome = remaining

    return genome


def median_root_to_leaf_lengths(tree):
    return np.median([tree.distance(tree.root, leaf) for leaf in tree.get_terminals()])

def evolve_genome(tree, root_genome, per_gene_gain_rate, per_gene_loss_rate, per_gene_inv_rate, per_gene_trans_rate, gain_exp, loss_exp, inv_exp, trans_exp, next_gene_id_holder):
    genomes = {}
    num_genes = len(root_genome)
    median_path_length = median_root_to_leaf_lengths(tree)

    branch_gain_rate = (num_genes * per_gene_gain_rate) / median_path_length
    branch_loss_rate = (num_genes * per_gene_loss_rate) / median_path_length
    branch_inv_rate = (num_genes * per_gene_inv_rate) / median_path_length
    branch_trans_rate = (num_genes * per_gene_trans_rate) / median_path_length

    L0_ancestral = len(root_genome)

    genomes[tree.root] = root_genome.copy()
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            parent_genome = genomes[parent]
            child_genome = evolve_genome_branch(
                parent_genome, child.branch_length,
                branch_gain_rate, branch_loss_rate, branch_inv_rate, branch_trans_rate,
                gain_exp, loss_exp, inv_exp, trans_exp,
                next_gene_id_holder,
                L0_ancestral
            )
            genomes[child] = child_genome
    return genomes


# def synteny_blocks(genome1, genome2):
#     pos_in_B = {gene: idx for idx, gene in enumerate(genome2)}
#     blocks = []
#     i = 0
#     while i < len(genome1):
#         gene = genome1[i]
#         if gene in pos_in_B:
#             j = pos_in_B[gene]
#             block_len = 1
#             while i + block_len < len(genome1):
#                 next_gene = genome1[i + block_len]
#                 if next_gene in pos_in_B and pos_in_B[next_gene] == j + block_len:
#                     block_len += 1
#                 else:
#                     break
#             blocks.append(block_len)
#             i += block_len
#         else:
#             i += 1
#     return blocks

def get_real_genomes_from_cc(path):
    df = pd.read_csv(path, names=[
        "gene_ID", "genome_ID", "protein_ID", "protein_length",
        "atgc_cog_footprint", "atgc_cog_footprint_length",
        "atgc_cog_ID", "protein_cluster_ID", "match_class"
    ])
    df = df.dropna(subset=["atgc_cog_ID"])
    return df.groupby("genome_ID")["atgc_cog_ID"].apply(list).to_dict()


DATASET = "ATGC0070"
TREE_FILE = f"{DATASET}/atgc.iq.r.tre"
CC_FILE = f"{DATASET}/atgc.cc.csv"

tree = Phylo.read(TREE_FILE, "newick")
real_genomes = get_real_genomes_from_cc(CC_FILE)

# Compute median genome size
median_length = int(np.median([len(g) for g in real_genomes.values()]))

# Build artificial root genome
root_genome = [f"cog.{i+1}" for i in range(median_length)]

next_gene_id_holder = [len(root_genome) + 1]

# --- FIX: add missing comma, and pass exponents expected by evolve_genome ---
simulated_genomes = evolve_genome(
    tree, root_genome,
    per_gene_gain_rate=0.1,
    per_gene_loss_rate=0.1,
    per_gene_inv_rate=0.0,
    per_gene_trans_rate=0.0,
    gain_exp=2.0, loss_exp=2.0, inv_exp=2.0, trans_exp=2.0,
    next_gene_id_holder=next_gene_id_holder
)

# Collect synteny block lengths per pair in a dictionary

# --- UPDATE: keep old call style working by giving defaults for new args ---
def run_simulation(
    tree, root_genome,
    per_gene_gain_rate, per_gene_loss_rate, per_gene_inv_rate,
    per_gene_trans_rate=0.0,
    gain_exp=2.0, loss_exp=2.0, inv_exp=2.0, trans_exp=2.0
):
    next_gene_id_holder = [len(root_genome) + 1]
    genomes = evolve_genome(
        tree, root_genome,
        per_gene_gain_rate, per_gene_loss_rate, per_gene_inv_rate, per_gene_trans_rate,
        gain_exp, loss_exp, inv_exp, trans_exp,
        next_gene_id_holder
    )
    leaves = tree.get_terminals()
    print(leaves)
    pairwise_blocks = {}

    for leaf1, leaf2 in combinations(leaves, 2):
        g1 = genomes[leaf1].copy()
        g2 = genomes[leaf2].copy()
        blocks = findSyntenyReal(g1, g2)
        block_lengths = [abs(b[4]) for b in blocks if len(b) > 4]
        pairwise_blocks[(leaf1.name, leaf2.name)] = block_lengths
    print(pairwise_blocks)

    return pairwise_blocks

# optional plot
sim_pairs = run_simulation(tree, root_genome, per_gene_gain_rate=0.1, per_gene_loss_rate=0.1, per_gene_inv_rate=0.0)
all_sim_lengths = [length for blocks in sim_pairs.values() for length in blocks]
plt.figure(figsize=(10, 5))
plt.hist(all_sim_lengths, bins=range(1, max(all_sim_lengths)+2), color="orange", edgecolor="black", alpha=0.7)
plt.title("Distribution of Simulated Synteny Block Lengths")
plt.xlabel("Synteny Block Length")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()
