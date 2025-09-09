# from Bio import Phylo
# import pandas as pd
# from itertools import combinations
# import matplotlib.pyplot as plt
# from synteny_tools import findSyntenyReal
# from synteny_tools import findSyntenyReal_fast


# def get_real_genomes_from_cc(dataset):
#     df = pd.read_csv(f"ATGC{dataset}/atgc.cc.csv", names=[
#         "gene_ID", "genome_ID", "protein_ID", "protein_length",
#         "atgc_cog_footprint", "atgc_cog_footprint_length",
#         "cog_ID", "cls_ID", "match_class"
#     ])

#     df = df.dropna(subset=["cls_ID"])

#     return df.groupby("genome_ID")["cls_ID"].apply(list).to_dict()

# def convert_to_numeric(cls_list):
#     return [int(cls.replace("cls.", "")) for cls in cls_list]

# tree = Phylo.read("ATGC0070/atgc.iq.r.tre", "newick")
# genome_to_genes = get_real_genomes_from_cc("0070")
# leaf_names = [leaf.name for leaf in tree.get_terminals()]

# # Dict to store block lengths per genome pair

# def get_pair_blocks():
#     real_pairwise_blocks = {}
#     for genome1, genome2 in combinations(leaf_names, 2):
#         if genome1 not in genome_to_genes or genome2 not in genome_to_genes:
#             continue  # skip missing genomes

#         g1 = convert_to_numeric(genome_to_genes[genome1])
#         g2 = convert_to_numeric(genome_to_genes[genome2])

#         blocks = findSyntenyReal(g1.copy(), g2.copy())
#         block_lengths = [b[4] for b in blocks]
#         real_pairwise_blocks[(genome1, genome2)] = block_lengths
#     return real_pairwise_blocks


# #optional plot
# real_blocks = get_pair_blocks()

# all_lengths = [length for blocklist in real_blocks.values() for length in blocklist]
# plt.figure(figsize=(10, 5))
# plt.hist(all_lengths, bins=range(1, max(all_lengths) + 2), color="skyblue", edgecolor="black")
# plt.title("Distribution of Synteny Block Lengths (Real Genome Pairs)")
# plt.xlabel("Synteny Block Length")
# plt.ylabel("Frequency")
# plt.tight_layout()
# plt.show()

# data_setup_from_pty.py
# calc_real_blocks.py


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Phylo
from itertools import combinations, product
import numpy as np
import csv
import os

# keep your finder where it is
from synteny_tools import findSyntenyReal2   # must return [s1,e1,s2,e2,L]  (1-based; L > 0)

# ----------------------------
# Config (adjust as needed)
# ----------------------------
DATASET = "0064"
LABEL   = "cls"   # "cls" (recommended) or "cog"

PTY_PATH  = f"ATGC{DATASET}/atgc.pty"
TREE_PATH = f"ATGC{DATASET}/atgc.iq.r.tre"
OUTDIR    = "output"
OUTFILE   = os.path.join(OUTDIR, f"ATGC{DATASET} - blocks_priority.tab")

# ----------------------------
# PTY parsing (partitions)
# ----------------------------
def parse_pty_partitions(pty_path: str, label: str = "cls"):
    """
    Returns:
      dict[str, list[dict(part_id:str, genes:list[int])]]
    Keeps chromosome/contig boundaries as separate partitions in order.
    PTY lines are assumed to have genome_id at col[3], partition_id at col[4],
    and label tokens like 'cls.000123' and/or 'cog.000456' near the end.
    """
    label = label.lower()
    prefix = label + "."
    partitions_by_genome = {}

    with open(pty_path, "r") as fh:
        current_gid = None
        current_pid = None
        for raw in fh:
            if not raw or raw[0] == "#":
                continue
            cols = raw.strip().split()
            if len(cols) < 6:
                continue

            gid = cols[3]              # e.g. GCF_...
            pid = cols[4]              # e.g. NC_...

            # find the right label token from the end (robust to extra cols)
            labtok = None
            for t in reversed(cols):
                tl = t.lower()
                if tl.startswith(prefix):
                    labtok = t
                    break
            if labtok is None:
                # missing this label on the line -> skip
                continue

            try:
                lab_id = int(labtok.split(".")[1])
            except Exception:
                continue

            if gid not in partitions_by_genome:
                partitions_by_genome[gid] = []
                partitions_by_genome[gid].append({"part_id": pid, "genes": [lab_id]})
                current_gid, current_pid = gid, pid
            else:
                # start a new partition when gid or pid changes
                if gid != current_gid or pid != current_pid:
                    partitions_by_genome[gid].append({"part_id": pid, "genes": []})
                    current_gid, current_pid = gid, pid
                partitions_by_genome[gid][-1]["genes"].append(lab_id)

    return partitions_by_genome

# ----------------------------
# Utilities for masking indices
# ----------------------------
def _forward_span_indices(start_1based: int, length: int, n: int):
    """Positions going forward from start, circular, 0-based indices."""
    return [((start_1based - 1) + k) % n for k in range(length)]

def _backward_span_indices(start_1based: int, length: int, n: int):
    """Positions going backward from start, circular, 0-based indices."""
    return [((start_1based - 1) - k) % n for k in range(length)]

def _orientation_from_coords(s: int, e: int, L: int, n: int):
    """Infer orientation in genome2 from (s,e,L) on a circular contig (1-based)."""
    # distance going forward from s to e (inclusive)
    fwd = ((e - s) % n) + 1
    if fwd == L:
        return "forward"
    # distance going backward from s to e (inclusive)
    bwd = ((s - e) % n) + 1
    if bwd == L:
        return "inverted"
    # fallback (rare)
    return "forward"

# ----------------------------
# Core: pairwise across partitions with priorities
# ----------------------------
def blocks_for_genome_pair_with_priorities(partsA, partsB):
    """
    partsA/B: list of {'part_id':str, 'genes':list[int]}
    Returns list of records:
      {
        'ia': idx in partsA, 'ib': idx in partsB,
        'genomeA_part_id': str, 'genomeB_part_id': str,
        'blocks': list[[s1,e1,s2,e2,L]]
      }
    Priority rules:
      - Within a partition pair: longer blocks win (handled by finder + sorting here).
      - Between partition pairs: process pairs by descending product of sizes.
      - No reuse of positions across different partition pairs for this genome pair.
    """
    masksA = [np.zeros(len(p["genes"]), dtype=bool) for p in partsA]
    masksB = [np.zeros(len(p["genes"]), dtype=bool) for p in partsB]

    # rank partition pairs by product of sizes (descending)
    pair_order = sorted(
        product(range(len(partsA)), range(len(partsB))),
        key=lambda ij: len(partsA[ij[0]]["genes"]) * len(partsB[ij[1]]["genes"]),
        reverse=True,
    )

    results = []

    for ia, ib in pair_order:
        genesA = partsA[ia]["genes"]
        genesB = partsB[ib]["genes"]
        if len(genesA) == 0 or len(genesB) == 0:
            continue

        # build working copies with already-used positions set to -1
        workA = [g if not masksA[ia][k] else -1 for k, g in enumerate(genesA)]
        workB = [g if not masksB[ib][k] else -1 for k, g in enumerate(genesB)]

        if all(x == -1 for x in workA) or all(x == -1 for x in workB):
            continue

        # find blocks (circular + inversion inside your function)
        blocks = findSyntenyReal2(workA.copy(), workB.copy())  # 1-based

        if not blocks:
            continue

        # longer-first within pair (ties keep upstream order)
        blocks.sort(key=lambda b: b[4], reverse=True)

        # persistently mark used positions to forbid reuse in later (smaller) pairs
        nA = len(genesA)
        nB = len(genesB)
        for (s1, e1, s2, e2, L) in blocks:
            idxA = _forward_span_indices(s1, L, nA)  # A is forward-coded in your routine
            orient = _orientation_from_coords(s2, e2, L, nB)
            idxB = (_forward_span_indices(s2, L, nB)
                    if orient == "forward"
                    else _backward_span_indices(s2, L, nB))
            masksA[ia][idxA] = True
            masksB[ib][idxB] = True

        results.append({
            "ia": ia,
            "ib": ib,
            "genomeA_part_id": partsA[ia]["part_id"],
            "genomeB_part_id": partsB[ib]["part_id"],
            "blocks": blocks,
        })

    return results

# ----------------------------
# Writer
# ----------------------------
def write_pair_results(writer, genomeA, genomeB, partsA, partsB, pair_records):
    """
    Writes:
      genomeA, genomeB, partA_id, partB_id, lenA, lenB
      <blocks...>
      <blank>
    """
    for rec in pair_records:
        ia, ib = rec["ia"], rec["ib"]
        pa, pb = partsA[ia], partsB[ib]
        writer.writerow([
            genomeA, genomeB,
            pa["part_id"], pb["part_id"],
            len(pa["genes"]), len(pb["genes"]),
        ])
        writer.writerows(rec["blocks"])
        writer.writerow([])

# ----------------------------
# High-level driver
# ----------------------------
def run_real_blocks(pty_path=PTY_PATH, tree_path=TREE_PATH, label=LABEL, out_path=OUTFILE):
    partitions = parse_pty_partitions(pty_path, label=label)

    # try to read the tree; if not present, just use all genomes from PTY
    try:
        tree = Phylo.read(tree_path, "newick")
        leaves = [t.name for t in tree.get_terminals()]
        present = [g for g in leaves if g in partitions]
        if not present:
            present = list(partitions.keys())
    except Exception:
        present = list(partitions.keys())

    present = sorted(present)

    # summary
    print(f"Genomes in PTY: {len(partitions)}")
    print(f"Used genomes  : {len(present)}\n")
    for gid in present:
        ps = partitions[gid]
        print(f"{gid}: {len(ps)} partitions")
        for i, p in enumerate(ps):
            print(f"  [{i}] {p['part_id']}  length={len(p['genes'])}")
        print()

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    pair_lengths = {}  # (gA,gB) -> list of block lengths
    with open(out_path, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        for gA, gB in combinations(present, 2):
            partsA = partitions[gA]
            partsB = partitions[gB]

            # compute blocks with priority + no-reuse across partition pairs
            chosen = blocks_for_genome_pair_with_priorities(partsA, partsB)

            # collect lengths for quick stats
            lengths = []
            for rec in chosen:
                lengths.extend([b[4] for b in rec["blocks"]])
            pair_lengths[(gA, gB)] = lengths

            # write to file
            write_pair_results(wr, gA, gB, partsA, partsB, chosen)

            nblocks = len(lengths)
            avg = (np.mean(lengths) if nblocks else 0.0)
            print(f"{gA} vs {gB}: {nblocks} blocks  (avg length={avg:.2f})")

    print(f"\nWrote: {out_path}")
    return pair_lengths

# ----------------------------
# Run
# ----------------------------

stats = run_real_blocks()
print(stats)

from collections import Counter

def print_length_distribution(pair_lengths, per_pair=False):
    """
    pair_lengths: dict[(genomeA, genomeB)] -> list[int]  (what run_real_blocks returns)

    If per_pair=False (default): prints a single global distribution like:
      1:345, 2:2305, 3:1120, ...
    If per_pair=True : prints one such line per genome pair.
    """
    if per_pair:
        for (gA, gB), lens in pair_lengths.items():
            c = Counter(lens)
            line = ", ".join(f"{L}:{c[L]}" for L in sorted(c))
            print(f"{gA} vs {gB} -> {line}")
    else:
        c = Counter()
        for lens in pair_lengths.values():
            c.update(lens)
        line = ", ".join(f"{L}:{c[L]}" for L in sorted(c))
        print(line)

# example usage:
# pair_lengths = run_real_blocks()
print_length_distribution(stats)          # global distribution
# print_length_distribution(pair_lengths, True)    # per-pair distribution

    # print(stats)  # uncomment if you want the dict of block-lengths per genome pair
