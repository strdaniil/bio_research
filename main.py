from Bio import Phylo
from io import StringIO
from ete3 import Tree
import random, numpy as np

# random example used earlier as test
# newick_str = "((A:0.3,B:0.3):0.2,C:0.5);"
# tree = Phylo.read(StringIO(newick_str), "newick")
# print(tree)
# Phylo.draw(tree)

# function to traverse tree, might be useful later
def traverse(clade, parent_name=""):
    name = clade.name or parent_name  
    print(f"Node {name} (branch length from parent = {clade.branch_length})")
    for child in clade.clades:
        traverse(child, parent_name=name)

# creates a random tree with 5 leaves since we have no data yet
t = Tree()
t.populate(5)  
newick_str = t.write(format=1)
tree = Phylo.read(StringIO(newick_str), "newick")


def evolve_genome(genome, branch_length, gain_rate, loss_rate, inv_rate):
    genome = genome.copy()  
    # total num of events
    total_rate = gain_rate + loss_rate + inv_rate
    N = np.random.poisson(total_rate * branch_length) 
    
    for _ in range(N):
        r = random.random() * total_rate
        if r < loss_rate:
            # gene loss
            if genome:
                idx = random.randrange(len(genome))
                genome.pop(idx)
        elif r < loss_rate + gain_rate:
            # gene gain
            new_gene = f"newGene{random.randrange(10000)}"
            idx = random.randrange(len(genome)+1)
            genome.insert(idx, new_gene)
        else:
            # inversion
            if len(genome) >= 2:
                i = random.randrange(len(genome))
                j = random.randrange(i, len(genome))
                segment = genome[i:j+1]
                segment.reverse()
                genome[i:j+1] = segment
    return genome

genomes = {} 

root_genome = [f"gene{i}" for i in range(1, 101)]
genomes[tree.root] = root_genome

for parent in tree.find_clades(order="level"):
    for child in parent.clades:
        parent_genome = genomes[parent]
        #these rates can be changed, for now just simple, inversion rate similar to the one you found in your previous research
        child_genome = evolve_genome(parent_genome, child.branch_length,
                                     gain_rate=0.4, loss_rate=0.4, inv_rate=0.1)
        genomes[child] = child_genome

def synteny_blocks(genome1, genome2):
    pos_in_B = {gene: idx for idx, gene in enumerate(genome2)}
    blocks = []
    i = 0
    while i < len(genome1):
        gene = genome1[i]
        if gene in pos_in_B:
            # potential start of a block
            j = pos_in_B[gene]
            block_len = 1
            # extend block as long if next genes match consecutively
            while i + block_len < len(genome1):
                next_gene = genome1[i + block_len]
                if next_gene in pos_in_B and pos_in_B[next_gene] == j + block_len:
                    block_len += 1
                else:
                    break
            blocks.append(block_len)
            i += block_len  # skip ahead past this block
        else:
            i += 1
    return blocks

# example of two leaves:
leaf1, leaf2 = tree.get_terminals()[0], tree.get_terminals()[1]
blocks_12 = synteny_blocks(genomes[leaf1], genomes[leaf2])
print("Blocks between", leaf1.name, "and", leaf2.name, ":", blocks_12)
print("Block length distribution:", {L: blocks_12.count(L) for L in set(blocks_12)})
