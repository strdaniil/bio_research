from scipy.stats import wasserstein_distance
from simulation import run_simulation, get_real_genomes_from_cc
from calc_real_data import get_pair_blocks
import numpy as np
from Bio import Phylo

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

leaves = tree.get_terminals()

best_distance = float('inf')
best_params = None

# for gain in np.linspace(0.05, 0.5, 10):
#     for loss in np.linspace(0.05, 0.5, 10):
#         sim_lengths = run_simulation(tree, root_genome, gain, loss, inv_rate=0.0, gain_genes=all_cls_ids)
#         dist = wasserstein_distance(real_lengths, sim_lengths)
#         print(f"Trying gain={gain:.3f}, loss={loss:.3f} -> Distance={dist:.4f}")
        
#         if dist < best_distance:
#             best_distance = dist
#             best_params = (gain, loss)

# print("Best parameters:", best_params, "Wasserstein distance:", best_distance)


real_pairs = get_pair_blocks()

sim_pairs = run_simulation(tree, root_genome, 0.1, 0.1, inv_rate=0.0)

#print(real_pairs, sim_pairs)


def find_error(realpairs, simpairs):
    total_distance = 0
    for pair in realpairs:
        if pair in simpairs:
            real_blocks = realpairs[pair]
            sim_blocks = simpairs[pair]
            dist = wasserstein_distance(real_blocks, sim_blocks)
            total_distance += dist
        else:
            print(f"Missing simulated data for pair: {pair}")
    return total_distance

print(f"Trying gain={0.1:.3f}, loss={0.1:.3f} -> Distance={find_error(real_pairs, sim_pairs):.4f}")


# real_blocks = real_pairs

# # Flatten all block lengths into single lists for plotting
# all_real_lengths = [length for blocks in real_blocks.values() for length in blocks]
# all_sim_lengths  = [length for blocks in sim_pairs.values() for length in blocks]

# import matplotlib.pyplot as plt

# # Plotting side-by-side histograms
# plt.figure(figsize=(12, 5))

# plt.subplot(1, 2, 1)
# plt.hist(all_real_lengths, bins=range(1, max(all_real_lengths)+2), color='blue', alpha=0.7, edgecolor='black')
# plt.title("Real Synteny Block Lengths")
# plt.xlabel("Block Length")
# plt.ylabel("Frequency")

# plt.subplot(1, 2, 2)
# plt.hist(all_sim_lengths, bins=range(1, max(all_sim_lengths)+2), color='orange', alpha=0.7, edgecolor='black')
# plt.title("Simulated Synteny Block Lengths")
# plt.xlabel("Block Length")
# plt.ylabel("Frequency")

# plt.tight_layout()
# plt.show()


#this is now the training part; will need to change, slow right now
from skopt import gp_minimize
from skopt.space import Real
from skopt.utils import use_named_args

# Define the search space
search_space = [
    Real(0.01, 0.5, name="gain"),
    Real(0.01, 0.5, name="loss")
]

best_distance = float("inf")
best_params = None

@use_named_args(search_space)
def objective(**params):
    gain = params["gain"]
    loss = params["loss"]
    
    sim_pairs = run_simulation(tree, root_genome, gain, loss, inv_rate=0.0)
    dist = find_error(sim_pairs, real_pairs)
    print(f"Trying gain={gain:.3f}, loss={loss:.3f} -> Distance={dist:.4f}")
    
    global best_distance, best_params
    if dist < best_distance:
        best_distance = dist
        best_params = (gain, loss)
    
    return dist

# Run Bayesian optimization
result = gp_minimize(
    func=objective,
    dimensions=search_space,
    acq_func="EI",  # Expected Improvement
    n_calls=30,
    n_initial_points=5,
    random_state=42
)

print("Best parameters:", best_params, "Wasserstein distance:", best_distance)


