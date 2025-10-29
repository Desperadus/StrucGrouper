#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Build a UMAP embedding from Foldseek kNN results (no N^2 matrix).
# Input: knn_raw.m8 produced by "foldseek easy-search ... --format-output ..."

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap

def parse_args():
    p = argparse.ArgumentParser(description="UMAP from Foldseek kNN TSV")
    p.add_argument("--in", dest="in_tsv", required=True,
                   help="Foldseek TSV (knn_raw.m8)")
    p.add_argument("--out_prefix", default="umap",
                   help="Output prefix (default: umap)")
    p.add_argument("--k", type=int, default=50,
                   help="Neighbors per query to keep (must match or be ≤ --max-seqs)")
    p.add_argument("--umap_n_neighbors", type=int, default=30,
                   help="UMAP n_neighbors (≤ k+1 because we add self)")
    p.add_argument("--min_tm", type=float, default=0.0,
                   help="Optional TM threshold to keep edges (0–1)")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    return p.parse_args()

def main():
    args = parse_args()
    in_tsv = args.in_tsv
    K = int(args.k)
    UMAP_K = int(args.umap_n_neighbors)
    assert UMAP_K <= (K + 1), "UMAP n_neighbors must be ≤ k+1 (we add self as neighbor 0)"
    min_tm = float(args.min_tm)

    # Try to infer which columns you have based on our recommended --format-output
    # Preferred set: query,target,alnlen,qcov,tcov,evalue,score,alntmscore
    # Fallback set:  query,target,alnlen,qcov,tcov,evalue,score
    # Foldseek does not write a header; we assign names by position.
    base_cols_pref = ["query","target","alnlen","qcov","tcov","evalue","score","alntmscore"]
    base_cols_fbk  = ["query","target","alnlen","qcov","tcov","evalue","score"]

    # Read with more columns; then trim to available count.
    df = pd.read_csv(in_tsv, sep="\t", header=None, dtype=str, engine="c")
    if df.shape[1] >= len(base_cols_pref):
        df.columns = base_cols_pref + [f"extra_{i}" for i in range(df.shape[1]-len(base_cols_pref))]
    else:
        df.columns = base_cols_fbk[:df.shape[1]]

    # Coerce numeric fields that might exist
    for c in ["alnlen","qcov","tcov","evalue","score","alntmscore"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Compute TM-like similarity in [0,1]
    tm = None
    if "alntmscore" in df.columns:
        tm = df["alntmscore"]
    else:
        candidates = []
        if "score" in df.columns:
            # In TM-align mode: score ~ qTM * 100 (doc/issue notes)
            candidates.append(df["score"] / 100.0)
        if "evalue" in df.columns:
            # In TM-align mode some versions place TM or mean(TM) into evalue; keep only plausible 0–1 values
            ev = df["evalue"]
            ev_valid = ev.where((ev >= 0.0) & (ev <= 1.0))
            candidates.append(ev_valid)
        # rowwise max over available candidates
        tm = pd.concat(candidates, axis=1).max(axis=1)

    df["tm"] = tm
    df = df.dropna(subset=["tm"])

    # Remove self hits if present later, but we'll still need IDs; keep queries/targets union
    all_ids = pd.Index(sorted(set(df["query"]) | set(df["target"])))
    id2ix = {sid: i for i, sid in enumerate(all_ids)}
    N = len(all_ids)

    # Filter edges: remove self & low TM
    df = df[df["query"] != df["target"]]
    if min_tm > 0:
        df = df[df["tm"] >= min_tm]

    # Keep top-K per query by TM
    df = df.sort_values(["query", "tm"], ascending=[True, False])
    df = df.groupby("query", sort=False).head(K).reset_index(drop=True)

    # Build (indices, distances) arrays for UMAP precomputed_kNN
    # Shape: (N, K+1) because neighbor 0 must be self with distance 0
    indices = np.full((N, K + 1), -1, dtype=np.int32)
    dists   = np.full((N, K + 1), np.inf, dtype=np.float32)

    # Set self neighbors
    for i in range(N):
        indices[i, 0] = i
        dists[i, 0]   = 0.0

    # Fill neighbors from Foldseek edges
    for q, sub in df.groupby("query", sort=False):
        i = id2ix[q]
        js = sub["target"].map(id2ix).to_numpy(dtype=np.int32, copy=False)
        # Distance = 1 - TM (clamp to [0,1])
        ds = (1.0 - np.clip(sub["tm"].to_numpy(dtype=np.float32, copy=False), 0.0, 1.0)).astype(np.float32)
        n  = min(K, js.shape[0])
        indices[i, 1:1+n] = js[:n]
        dists[i,   1:1+n] = ds[:n]

    # Run UMAP with precomputed k-NN (no dense matrix)
    umap_k = min(UMAP_K, K + 1)
    reducer = umap.UMAP(
        n_neighbors=umap_k,
        min_dist=0.3,
        random_state=args.seed,
        precomputed_knn=(indices, dists),  # official fast path
        low_memory=True,
        verbose=True,
    )

    # Dummy X (not used for neighbor finding when precomputed_knn is supplied)
    X_dummy = np.zeros((N, 1), dtype=np.float32)
    Y = reducer.fit_transform(X_dummy)

    # Save outputs
    out_prefix = args.out_prefix
    coords_path = f"{out_prefix}_coords.tsv"
    fig_path    = f"{out_prefix}.png"

    pd.DataFrame({"id": all_ids, "umap_x": Y[:, 0], "umap_y": Y[:, 1]}) \
      .to_csv(coords_path, sep="\t", index=False)

    plt.figure(figsize=(8, 6))
    plt.scatter(Y[:, 0], Y[:, 1], s=2)
    plt.xlabel("UMAP-1")
    plt.ylabel("UMAP-2")
    plt.title(f"Foldseek kNN UMAP (N={N}, k={umap_k})")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=200)
    print(f"Wrote: {coords_path} and {fig_path}")

if __name__ == "__main__":
    main()

