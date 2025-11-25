#!/usr/bin/env python3
import argparse
import h5py
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# UMAP is in umap-learn
try:
    import umap
except ImportError as e:
    raise SystemExit(
        "umap-learn is required. Install with: pip install umap-learn"
    ) from e

import matplotlib.pyplot as plt
from pathlib import Path
import sys


def load_h5_embeddings(h5_path: str):
    ids = []
    vecs = []
    with h5py.File(h5_path, "r") as f:
        for key in f.keys():
            arr = np.array(f[key])
            if arr.ndim != 1:
                raise ValueError(
                    f"Dataset {key} is not 1D (got shape {arr.shape})")
            ids.append(key)
            vecs.append(arr)
    X = np.vstack(vecs)
    return np.array(ids), X


def fit_pca_2d(X: np.ndarray, whiten: bool = False, random_state: int = 42):
    # Standardize before PCA for stability
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xs = scaler.fit_transform(X)
    pca = PCA(n_components=2, random_state=random_state, whiten=whiten)
    coords = pca.fit_transform(Xs)
    return coords, pca, scaler


def fit_umap_2d(X: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.1, random_state: int = 42):
    # Standardize before UMAP
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xs = scaler.fit_transform(X)
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric="euclidean",
        random_state=random_state,
    )
    coords = reducer.fit_transform(Xs)
    return coords, reducer, scaler


def scatter_and_save(xy: np.ndarray, ids: np.ndarray, title: str, out_png: Path, s: float = 10.0):
    plt.figure(figsize=(6, 5), dpi=150)
    plt.scatter(xy[:, 0], xy[:, 1], s=s)
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description="Compute PCA and UMAP (2D) on UniProt embeddings in HDF5.")
    ap.add_argument("h5_file", help="Path to HDF5 file (no extension needed).")
    ap.add_argument("--outprefix", default="embeds",
                    help="Prefix for output files (default: embeds)")
    ap.add_argument("--neighbors", type=int, default=15,
                    help="UMAP n_neighbors (default: 15)")
    ap.add_argument("--min-dist", type=float, default=0.1,
                    help="UMAP min_dist (default: 0.1)")
    ap.add_argument("--whiten", action="store_true",
                    help="Whiten in PCA (optional)")
    args = ap.parse_args()

    h5_path = args.h5_file
    if not Path(h5_path).exists():
        # try opening without extension if user omitted it
        if Path(h5_path + ".h5").exists():
            h5_path = h5_path + ".h5"
        else:
            print(f"File not found: {args.h5_file}", file=sys.stderr)
            sys.exit(1)

    ids, X = load_h5_embeddings(h5_path)
    if len(ids) != X.shape[0]:
        raise RuntimeError(
            "Mismatch between number of IDs and embedding rows.")

    # PCA 2D
    pca_xy, pca_model, pca_scaler = fit_pca_2d(X, whiten=args.whiten)
    pca_df = pd.DataFrame(
        {"id": ids, "umap_x": pca_xy[:, 0], "umap_y": pca_xy[:, 1]})
    pca_tsv = Path(f"{args.outprefix}_pca.tsv")
    pca_df.to_csv(pca_tsv, sep="\t", index=False)

    # UMAP 2D
    umap_xy, umap_model, umap_scaler = fit_umap_2d(
        X, n_neighbors=args.neighbors, min_dist=args.min_dist
    )
    umap_df = pd.DataFrame(
        {"id": ids, "umap_x": umap_xy[:, 0], "umap_y": umap_xy[:, 1]})
    umap_tsv = Path(f"{args.outprefix}_umap.tsv")
    umap_df.to_csv(umap_tsv, sep="\t", index=False)

    # Plots
    scatter_and_save(pca_xy, ids, "PCA (2D) of embeddings",
                     Path(f"{args.outprefix}_pca.png"))
    scatter_and_save(umap_xy, ids, "UMAP (2D) of embeddings",
                     Path(f"{args.outprefix}_umap.png"))

    # Also create a quick combined TSV in the format you asked for (UMAP)
    # Already saved as <prefix>_umap.tsv with: id, umap_x, umap_y

    print(f"Saved: {pca_tsv}, {umap_tsv}")
    print(f"Saved: {args.outprefix}_pca.png, {args.outprefix}_umap.png")


if __name__ == "__main__":
    main()
