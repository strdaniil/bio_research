from scipy.stats import wasserstein_distance
from simulation import run_simulation, get_real_genomes_from_cc
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

all_cls_ids = list({gene for genes in real_genomes.values() for gene in genes})

best_distance = float('inf')
best_params = None

import ast

with open("real_lengths.txt", "r") as f:
    real_lengths = ast.literal_eval(f.read())



# for gain in np.linspace(0.05, 0.5, 10):
#     for loss in np.linspace(0.05, 0.5, 10):
#         sim_lengths = run_simulation(tree, root_genome, gain, loss, inv_rate=0.0, gain_genes=all_cls_ids)
#         dist = wasserstein_distance(real_lengths, sim_lengths)
#         print(f"Trying gain={gain:.3f}, loss={loss:.3f} -> Distance={dist:.4f}")
        
#         if dist < best_distance:
#             best_distance = dist
#             best_params = (gain, loss)

# print("Best parameters:", best_params, "Wasserstein distance:", best_distance)

sim_lengths = run_simulation(tree, root_genome, 0.1, 0.1, inv_rate=0.0, gain_genes=all_cls_ids)
dist = wasserstein_distance(real_lengths, sim_lengths)
print(f"Trying gain={0.1:.3f}, loss={0.1:.3f} -> Distance={dist:.4f}")

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
    
    sim_lengths = run_simulation(tree, root_genome, gain, loss, inv_rate=0.0, gain_genes=all_cls_ids)
    dist = wasserstein_distance(real_lengths, sim_lengths)
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


