from Bio import Phylo
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
from synteny_tools import findSyntenyReal


def get_real_genomes_from_cc(dataset):
    df = pd.read_csv(f"ATGC{dataset}/atgc.cc.csv", names=[
        "gene_ID", "genome_ID", "protein_ID", "protein_length",
        "atgc_cog_footprint", "atgc_cog_footprint_length",
        "cog_ID", "cls_ID", "match_class"
    ])

    # Filter rows where cls_ID exists (ignore missing annotations)
    df = df.dropna(subset=["cls_ID"])

    # Reconstruct gene order per genome
    return df.groupby("genome_ID")["cls_ID"].apply(list).to_dict()

def convert_to_numeric(cls_list):
    return [int(cls.replace("cls.", "")) for cls in cls_list]

tree = Phylo.read("ATGC0070/atgc.iq.r.tre", "newick")
genome_to_genes = get_real_genomes_from_cc("0070")
leaf_names = [leaf.name for leaf in tree.get_terminals()]

# List to store all block lengths from all genome pairs
all_block_lengths = []

# Compare all leaf-to-leaf genome pairs; doesnt store the syntenty blocks, but we can change this later
# not needed for my task as of now
for genome1, genome2 in combinations(leaf_names, 2):
    if genome1 not in genome_to_genes or genome2 not in genome_to_genes:
        continue  # skip missing genomes
    
    g1 = convert_to_numeric(genome_to_genes[genome1])
    g2 = convert_to_numeric(genome_to_genes[genome2])
    
    blocks = findSyntenyReal(g1.copy(), g2.copy())
    all_block_lengths.extend([b[4] for b in blocks])

# Plot the distribution
# plt.hist(all_block_lengths, bins=range(1, max(all_block_lengths)+2), density=True, edgecolor='black')
# plt.xlabel("Synteny Block Length")
# plt.ylabel("Normalized Frequency")
# plt.title("Distribution of Synteny Blocks (All Leaf Pairs)")
# plt.tight_layout()
# plt.show()


# with open("real_lengths.txt", "w") as f:
#     f.write(str(all_block_lengths))
