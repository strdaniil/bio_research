# --- sim_real_comparison.py additions ---------------------------------------
import os, csv, math, json, time
from datetime import datetime
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np
from Bio import Phylo
from scipy.stats import wasserstein_distance

# If you have a fast version, import/alias it here:
# from synteny_tools import findSyntenyReal
# or: from calc_real_data import findSyntenyReal_fast as findSyntenyReal
from synteny_tools import findSyntenyReal2

from simulation import run_simulation  # expects your signature
import pandas as pd


# ---------------------------
# Helpers (self-contained)
# ---------------------------
def _convert_to_numeric(cls_list):
    # Accepts like ["cls.123", "cls.45"] or ["cog.123", ...] -> ints
    out = []
    for x in cls_list:
        if isinstance(x, int):
            out.append(int(x))
        else:
            s = str(x)
            if "." in s:
                s = s.split(".", 1)[1]
            out.append(int(s))
    return out


def _load_real_genomes_from_cc(cc_path):
    """
    Reads ATGC cc file into dict: genome_id -> [cls/cog IDs as strings]
    Works for both 'cls_ID' or 'atgc_cog_ID' column names.
    """
    df = pd.read_csv(
        cc_path,
        names=[
            "gene_ID", "genome_ID", "protein_ID", "protein_length",
            "atgc_cog_footprint", "atgc_cog_footprint_length",
            "atgc_cog_ID", "protein_cluster_ID", "match_class"
        ],
        dtype=str,
    )
    # Back-compat if older schema uses "cls_ID"
    if "atgc_cog_ID" not in df.columns or df["atgc_cog_ID"].isnull().all():
        df = pd.read_csv(
            cc_path,
            names=[
                "gene_ID", "genome_ID", "protein_ID", "protein_length",
                "atgc_cog_footprint", "atgc_cog_footprint_length",
                "cog_ID", "cls_ID", "match_class"
            ],
            dtype=str,
        )
        key_col = "cls_ID"
    else:
        key_col = "atgc_cog_ID"

    df = df.dropna(subset=[key_col])
    return df.groupby("genome_ID")[key_col].apply(list).to_dict()


def _lengths_to_pmf(lengths):
    """
    list[int] -> (vals[np.int64], probs[np.float64]) for discrete pmf.
    Returns ([],[]) if empty.
    """
    c = Counter(lengths)
    tot = sum(c.values())
    if tot == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    vals = np.fromiter(sorted(c.keys()), dtype=int)
    probs = np.fromiter((c[v] / tot for v in vals), dtype=float)
    return vals, probs


def _median_root_to_leaf(tree):
    return float(np.median([tree.distance(tree.root, lf) for lf in tree.get_terminals()]))


def _build_real_pmfs(tree_path, cc_path, synteny_finder=findSyntenyReal2):
    """
    Compute real SBL pmfs for each cognate pair of leaves in the tree.
    Returns dict { (nameA,nameB): (vals, probs) } with sorted pair keys.
    """
    tree = Phylo.read(tree_path, "newick")
    leaf_names = [lf.name for lf in tree.get_terminals()]
    genomes = _load_real_genomes_from_cc(cc_path)

    real_pmfs = {}
    for a, b in combinations(leaf_names, 2):
        if a not in genomes or b not in genomes:
            continue
        g1 = _convert_to_numeric(genomes[a])
        g2 = _convert_to_numeric(genomes[b])
        blocks = synteny_finder(g1.copy(), g2.copy())
        lens = [int(abs(bk[4])) for bk in blocks if len(bk) >= 5]
        pair = tuple(sorted((a, b)))
        real_pmfs[pair] = _lengths_to_pmf(lens)
    return real_pmfs, tree


def _make_root_genome(root_mode, tree, cc_path, real_genomes=None):
    """
    root_mode: "median_synthetic" or "from_info_tab"
    Returns list[int] root genome (integer labels, 1..L0 or from real genome IDs).
    """
    if real_genomes is None:
        real_genomes = _load_real_genomes_from_cc(cc_path)

    if root_mode == "from_info_tab":
        info_tab = os.path.join(os.path.dirname(cc_path), "atgc.info.tab")
        if not os.path.exists(info_tab):
            raise FileNotFoundError(f"Root mode 'from_info_tab' but {info_tab} not found.")
        with open(info_tab, "r") as f:
            first_line = f.readline().strip()
            gid = first_line.split()[0]
        root_seq = _convert_to_numeric(real_genomes[gid])
        return root_seq

    # default: median_synthetic
    med_len = int(np.median([len(_convert_to_numeric(v)) for v in real_genomes.values()]))
    return list(range(1, med_len + 1))


# ---------------------------
# 1) Standalone simulation-evaluation entry point
# ---------------------------
def simulate_and_score(
    atgc_dir: str,
    tree_path: str,
    cc_path: str,
    root_mode: str,
    per_gene_gain: float,
    per_gene_loss: float,
    per_gene_inv: float,
    per_gene_trans: float,
    exp_gain: float,
    exp_loss: float,
    exp_inv: float,
    exp_trans: float,
    n_runs: int,
    seed: int | None = None,
    synteny_finder=findSyntenyReal2,
):
    """
    Runs n_runs simulations on the given tree/root genome, pools simulated SBL counts per pair,
    builds pmfs once, computes Wasserstein-1 vs real pmfs for each pair, and sums across pairs.

    Returns:
        dict with keys:
          'sum_w1', 'avg_w1', 'n_pairs', 'skipped_real', 'skipped_sim',
          'params' (echo), plus 'dataset', 'median_root_to_leaf', 'timestamp'
    """
    # Real pmfs (and the tree to get leaves & med path length)
    real_pmfs, tree = _build_real_pmfs(tree_path, cc_path, synteny_finder=synteny_finder)

    # Root genome
    real_genomes = _load_real_genomes_from_cc(cc_path)
    root_genome = _make_root_genome(root_mode, tree, cc_path, real_genomes=real_genomes)

    # Pool simulated counts across runs
    pooled_sim_counts: dict[tuple[str, str], Counter] = defaultdict(Counter)

    # Ensure deterministic-ish naming
    leaves = list(tree.get_terminals())
    leaf_names = [lf.name for lf in leaves]

    # Run simulations
    for r in range(n_runs):
        sim_pairs = run_simulation(
            tree,
            root_genome,
            per_gene_gain_rate=per_gene_gain,
            per_gene_loss_rate=per_gene_loss,
            per_gene_inv_rate=per_gene_inv,
            per_gene_trans_rate=per_gene_trans,
            gain_exp=exp_gain,
            loss_exp=exp_loss,
            inv_exp=exp_inv,
            trans_exp=exp_trans,
        )
        # Pool counts; standardize pair key ordering
        for (a, b), lens in sim_pairs.items():
            pair = tuple(sorted((a, b)))
            if lens:
                pooled_sim_counts[pair].update(int(x) for x in lens)

    # Compare
    total_w1 = 0.0
    compared = 0
    skipped_real = 0
    skipped_sim = 0

    for pair, (r_vals, r_probs) in real_pmfs.items():
        if r_vals.size == 0:
            skipped_real += 1
            continue
        sim_counts = pooled_sim_counts.get(pair)
        if not sim_counts or sum(sim_counts.values()) == 0:
            skipped_sim += 1
            continue
        s_vals = np.fromiter(sorted(sim_counts.keys()), dtype=int)
        s_probs = np.fromiter((sim_counts[v] / sum(sim_counts.values()) for v in s_vals), dtype=float)

        dist = wasserstein_distance(r_vals, s_vals, u_weights=r_probs, v_weights=s_probs)
        total_w1 += float(dist)
        compared += 1

    avg_w1 = (total_w1 / compared) if compared > 0 else float("nan")
    med_path = _median_root_to_leaf(tree)
    dataset = os.path.basename(os.path.normpath(atgc_dir))

    return {
        "sum_w1": total_w1,
        "avg_w1": avg_w1,
        "n_pairs": compared,
        "skipped_real": skipped_real,
        "skipped_sim": skipped_sim,
        "dataset": dataset,
        "median_root_to_leaf": med_path,
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "params": {
            "root_mode": root_mode,
            "per_gene_gain": per_gene_gain,
            "per_gene_loss": per_gene_loss,
            "per_gene_inv": per_gene_inv,
            "per_gene_trans": per_gene_trans,
            "exp_gain": exp_gain,
            "exp_loss": exp_loss,
            "exp_inv": exp_inv,
            "exp_trans": exp_trans,
            "n_runs": n_runs,
            "seed": seed,
            "tree_path": tree_path,
            "cc_path": cc_path,
        },
    }


# ---------------------------
# 2) 21×21 grid for the 2‑parameter (rf, rt) case
# ---------------------------
def run_grid_2d(
    atgc_dir: str,
    tree_filename: str = "atgc.iq.r.tre",
    cc_filename: str = "atgc.cc.csv",
    root_mode: str = "median_synthetic",
    rf_min: float = 1e-2,
    rf_max: float = 1e2,
    rt_min: float = 1e-3,
    rt_max: float = 1e1,
    points: int = 21,           # 21 points -> ~5 per decade over 4 decades
    n_runs: int = 4,            # start modest; you can raise later
    seed: int | None = 42,
    out_csv: str | None = None,
    huge_exp: float = 1e9,      # force k=1
    quiet: bool = False,
    synteny_finder=findSyntenyReal2,
):
    """
    Balanced flux: gain=loss=rf; inversion=0; translocation=rt; exponents huge (k=1).
    Writes one CSV row per point with results + params. Returns path to CSV.
    """
    atgc_dir = os.path.abspath(atgc_dir)
    tree_path = os.path.join(atgc_dir, tree_filename)
    cc_path = os.path.join(atgc_dir, cc_filename)
    if out_csv is None:
        out_csv = os.path.join(atgc_dir, f"{os.path.basename(atgc_dir)}_grid2d_results.csv")

    rf_vals = np.logspace(math.log10(rf_min), math.log10(rf_max), points)
    rt_vals = np.logspace(math.log10(rt_min), math.log10(rt_max), points)

    # Prepare CSV
    new_file = not os.path.exists(out_csv)
    with open(out_csv, "a", newline="") as fh:
        w = csv.writer(fh)
        if new_file:
            w.writerow([
                "dataset", "root_mode", "rf", "rt", "n_runs",
                "sum_w1", "avg_w1", "n_pairs", "skipped_real", "skipped_sim",
                "median_root_to_leaf", "seed", "tree_path", "cc_path", "timestamp"
            ])

        total = len(rf_vals) * len(rt_vals)
        done = 0
        t0 = time.time()

        for rf in rf_vals:
            for rt in rt_vals:
                res = simulate_and_score(
                    atgc_dir=atgc_dir,
                    tree_path=tree_path,
                    cc_path=cc_path,
                    root_mode=root_mode,
                    per_gene_gain=float(rf),
                    per_gene_loss=float(rf),
                    per_gene_inv=0.0,
                    per_gene_trans=float(rt),
                    exp_gain=huge_exp,
                    exp_loss=huge_exp,
                    exp_inv=huge_exp,
                    exp_trans=huge_exp,
                    n_runs=n_runs,
                    seed=seed,
                    synteny_finder=synteny_finder,
                )

                w.writerow([
                    res["dataset"], root_mode, f"{rf:.8g}", f"{rt:.8g}", n_runs,
                    f"{res['sum_w1']:.10g}", f"{res['avg_w1']:.10g}", res["n_pairs"],
                    res["skipped_real"], res["skipped_sim"],
                    f"{res['median_root_to_leaf']:.6g}", seed, tree_path, cc_path, res["timestamp"]
                ])
                fh.flush()

                done += 1
                if not quiet:
                    elapsed = time.time() - t0
                    print(f"[grid] rf={rf:.3g}, rt={rt:.3g}  -> sumW1={res['sum_w1']:.3f} "
                          f"({done}/{total}, {elapsed/60:.1f} min elapsed)")

    if not quiet:
        print(f"[grid] wrote results to: {out_csv}")
    return out_csv


# ---------------------------
# Optional CLI entry points
# ---------------------------
if __name__ == "__main__":
    # Example single-point call:
    # r = simulate_and_score(
    #     atgc_dir="ATGC0070",
    #     tree_path="ATGC0070/atgc.iq.r.tre",
    #     cc_path="ATGC0070/atgc.cc.csv",
    #     root_mode="median_synthetic",
    #     per_gene_gain=0.1, per_gene_loss=0.1, per_gene_inv=0.0, per_gene_trans=0.1,
    #     exp_gain=1e9, exp_loss=1e9, exp_inv=1e9, exp_trans=1e9,
    #     n_runs=4, seed=42,
    # )
    # print(json.dumps(r, indent=2))

    # Example 21×21 grid:
    run_grid_2d(atgc_dir="ATGC0070", n_runs=4, seed=42, quiet=False)
