# StrucGrouper

Run Foldseek-based structural nearest-neighbour searches on AlphaFold (or other) models, collapse the results into a 2D UMAP embedding, and explore the proteins interactively in a browser with live structure previews and UniProt lookups.

![UMAP overview](umap.png)

## Highlights
- Reproducible Foldseek pipeline for building top-*k* structural neighbourhood graphs from a directory of PDB/mmCIF files (`create_foldseek_alignments.sh`).
- Fast UMAP projection that consumes foldseek `.m8` output without materialising dense distance matrices (`foldseek_umap.py`).
- Embedding + plotting utilities for HDF5 feature vectors such as UniProt/AlphaFold embeddings (`embed_reduce.py`).
- Browser-based UMAP explorer (`index.html`) with Plotly + 3Dmol.js: click any point to stream the corresponding AlphaFold model and fetch UniProt metadata.
- Helper scripts for fixing AlphaFold identifiers and expanding UniProt taxonomic lineage metadata (`tools/af2uniprot.sh`, `tools/expand_taxy.py`).
- Sample data (`data/`) covering lipocalin neighbourhoods and UniProt GO:0005549 proteins to test-drive the workflow end-to-end.

## Repository Layout
- `create_foldseek_alignments.sh` – build a Foldseek database and emit per-query top-*k* neighbours.
- `foldseek_umap.py` – convert Foldseek `.m8` results into a 2D UMAP embedding (with optional TM-score filtering).
- `embed_reduce.py` – PCA + UMAP projections for HDF5 embeddings downloaded from Uniprot, including static scatter plots.
- `index.html` – standalone AlphaFold UMAP Explorer (Plotly scatter + 3Dmol.js structure viewer).
- `tools/` – small CLI helpers (`af2uniprot.sh`, `expand_taxy.py`).

## Getting Started

### 1. Create a Conda environment
```bash
conda env create -f environment_raw.yml
conda activate strucgrouper
```
*Minimal stack*: `environment_raw.yml`.

### 2. Prepare structural neighbours with Foldseek
```bash
./create_foldseek_alignments.sh \
  -i afdb_mouse_lipocalin_pdb \        # directory with PDB/mmCIF files
  -o knn_raw_mouse_lipocalin.m8 \      # output edge list
  -T foldseek_tmp \                    # scratch directory (auto-created)
  -k 50                                # neighbours per query
```
- By default the query set equals the database (`-q`), so the command performs an all-vs-all search.
- Uses TM-align mode (`--alignment-type 1`) and requests `alntmscore`; older Foldseek builds fall back to `score`.
- Add `-f` to rebuild an existing Foldseek database.

### 3. Build a 2D embedding from the Foldseek output
```bash
python foldseek_umap.py \
  --in knn_raw_mouse_lipocalin.m8 \
  --out_prefix data/lipocalins \
  --k 50 \
  --umap_n_neighbors 30 \
  --min_tm 0.4
```
Outputs:
- `<prefix>_coords.tsv` (`id`, `umap_x`, `umap_y`) for plotting.
- `<prefix>.png` quick-look scatter plot.

Identifiers in Foldseek output are preserved. `--min_tm` lets you drop low-confidence edges (0–1).

### 4. Optional: reduce precomputed embeddings
For UniProt HDF5 embeddings (one dataset per accession):
```bash
python embed_reduce.py uniprot_embeddings.h5 \
  --outprefix data/uniprot \
  --neighbors 30 \
  --min-dist 0.1
```
Produces both PCA (`*_pca.tsv`/`.png`) and UMAP (`*_umap.tsv`/`.png`) summaries.

### 5. Tidy identifiers and metadata
- `tools/af2uniprot.sh` rewrites AlphaFold IDs (e.g. `AF-A0A0A6YW05-F1-model_v6`) back to UniProt accessions, which helps when joining with UniProt metadata TSVs.
  ```bash
  cat data/protein_umap_coords.tsv | tools/af2uniprot.sh > data/uniprot_umap_coords.tsv
  ```
- `tools/expand_taxy.py` normalises the “Taxonomic lineage” column in UniProt exports into one column per lineage type:
  ```bash
  python tools/expand_taxy.py data/uniprotkb_GO_0005549_2025_10_30.tsv
  # -> data/uniprotkb_GO_0005549_2025_10_30_extended.tsv
  ```

## AlphaFold UMAP Explorer
1. Serve the repository (browser security blocks `file://` AJAX):
   ```bash
   python -m http.server 8000
   ```
   Visit `http://localhost:8000/index.html`.
2. Upload a coordinates TSV (`id`, `umap_x`, `umap_y`). Example: `data/umap_coords_lipocalins.tsv`.
3. (Optional) Upload metadata TSV where the first column matches the identifiers in the coordinate file. Example: `data/uniprotkb_GO_0005549_2025_10_30_extended.tsv`.
4. Click points to stream AlphaFold models (via 3Dmol.js) and fetch UniProt metadata.

Tips
- Hold Shift to box/lasso select clusters in Plotly.
- Use the “Colour by” dropdown once metadata is loaded; each column becomes an option.
- If direct downloads from `alphafold.ebi.ac.uk` are blocked, proxy requests via `/fetch-structure?...` (the viewer already tries those endpoints if a suitable backend is present).

## License
This project is released under the MIT License (see `LICENSE`).
