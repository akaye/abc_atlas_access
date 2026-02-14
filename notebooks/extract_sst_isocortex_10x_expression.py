#!/usr/bin/env python3
"""
Extract receptor gene expression data for Sst neurons across isocortex from 10x data.

This script processes the large Isocortex h5ad files one at a time,
extracting only the 28 receptor genes for Sst cells, then saves the
result as a small CSV file for the notebook to load.

It handles disk space constraints by deleting each ~8-12 GB h5ad file
after extracting the subset (~57k cells × 28 genes).
"""
import os
import gc
import time
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from pathlib import Path

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# --- Configuration ---
download_base = Path('../../data/abc_atlas')
output_csv = Path('sst_isocortex_10x_expression.csv')
metadata_csv = Path('sst_isocortex_10x_metadata.csv')

# Receptor genes (mouse nomenclature)
serotonin_receptors = [
    'Htr1a', 'Htr1b', 'Htr1d', 'Htr1f',
    'Htr2a', 'Htr2b', 'Htr2c',
    'Htr3a', 'Htr3b',
    'Htr4', 'Htr5a', 'Htr5b', 'Htr6', 'Htr7'
]
norepinephrine_receptors = [
    'Adra1a', 'Adra1b', 'Adra1d',
    'Adra2a', 'Adra2b', 'Adra2c',
    'Adrb1', 'Adrb2', 'Adrb3'
]
dopamine_receptors = ['Drd1', 'Drd2', 'Drd3', 'Drd4', 'Drd5']
all_receptors = serotonin_receptors + norepinephrine_receptors + dopamine_receptors

# Sst subclasses
sst_subclasses = ['053 Sst Gaba', '056 Sst Chodl Gaba']

# Isocortex dissection regions
isocortex_rois = [
    'VIS', 'MOp', 'AI', 'PL-ILA-ORB', 'ACA',
    'TEa-PERI-ECT', 'SS-GU-VISC', 'RSP', 'AUD',
    'MO-FRP', 'VIS-PTLp', 'SSp', 'AUD-TEa-PERI-ECT'
]

# --- Initialize ---
print("Initializing ABC Atlas cache...")
abc_cache = AbcProjectCache.from_s3_cache(download_base)

# --- Load metadata ---
print("Loading cell metadata...")
cell = abc_cache.get_metadata_dataframe(
    directory='WMB-10X', file_name='cell_metadata', dtype={'cell_label': str}
)
cell.set_index('cell_label', inplace=True)
print(f"  Total cells: {len(cell):,}")

print("Loading gene metadata...")
gene = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')
gene.set_index('gene_identifier', inplace=True)

print("Loading taxonomy...")
cluster_details = abc_cache.get_metadata_dataframe(
    directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_pivoted',
    keep_default_na=False
)
cluster_details.set_index('cluster_alias', inplace=True)
cell_extended = cell.join(cluster_details, on='cluster_alias')

# --- Identify Sst cells in isocortex ---
sst_cells = cell_extended[
    (cell_extended['subclass'].isin(sst_subclasses)) &
    (cell_extended['region_of_interest_acronym'].isin(isocortex_rois))
].copy()
print(f"\nSst cells in isocortex: {len(sst_cells):,}")
print(f"\nBy dissection region:")
for roi, cnt in sst_cells.groupby('region_of_interest_acronym', observed=True).size().sort_values(ascending=False).items():
    print(f"  {roi}: {cnt:,}")
print(f"\nBy subclass:")
for sc, cnt in sst_cells.groupby('subclass', observed=True).size().sort_values(ascending=False).items():
    print(f"  {sc}: {cnt:,}")

# Save metadata
meta_cols = ['subclass', 'supertype', 'class', 'neurotransmitter',
             'region_of_interest_acronym', 'feature_matrix_label', 'dataset_label']
sst_cells[meta_cols].to_csv(metadata_csv)
print(f"\nSaved metadata: {metadata_csv}")

# --- Identify receptor genes ---
available_genes = gene[gene['gene_symbol'].isin(all_receptors)]
receptor_genes = [g for g in all_receptors if g in set(available_genes['gene_symbol'])]
gene_ensembl_ids = set(available_genes.index.tolist())
ensembl_to_symbol = dict(zip(available_genes.index, available_genes['gene_symbol']))
print(f"\nReceptor genes found: {len(receptor_genes)}")

# --- Extract expression data ---
sst_matrices = sst_cells.groupby('feature_matrix_label').size()
print(f"\nExpression matrices to process: {len(sst_matrices)}")
for mat, count in sst_matrices.items():
    print(f"  {mat}: {count:,} cells")

expression_frames = []

for i, matrix_label in enumerate(sst_matrices.index):
    dataset_label = sst_cells[
        sst_cells['feature_matrix_label'] == matrix_label
    ]['dataset_label'].iloc[0]
    file_name = f"{matrix_label}/log2"

    print(f"\n[{i+1}/{len(sst_matrices)}] Processing {file_name}...")
    t0 = time.time()

    file_path = abc_cache.get_file_path(directory=dataset_label, file_name=file_name)
    file_path_str = str(file_path)

    with h5py.File(file_path_str, 'r') as f:
        obs_cell_labels = f['obs']['cell_label'][:].astype(str)
        var_gene_ids = f['var']['gene_identifier'][:].astype(str)
        var_gene_symbols = f['var']['gene_symbol'][:].astype(str)

        n_cells_total = len(obs_cell_labels)
        n_genes_total = len(var_gene_ids)

        sst_cell_label_set = set(
            sst_cells[sst_cells['feature_matrix_label'] == matrix_label].index
        )
        cell_idx = np.array([j for j, name in enumerate(obs_cell_labels)
                             if name in sst_cell_label_set])
        gene_idx = np.array([j for j, gid in enumerate(var_gene_ids)
                             if gid in gene_ensembl_ids])

        print(f"  Total cells in matrix: {n_cells_total:,}")
        print(f"  Sst cells found: {len(cell_idx):,}")
        print(f"  Receptor genes found: {len(gene_idx)}")

        if len(cell_idx) == 0:
            print(f"  Skipping (no matching cells)")
        else:
            print(f"  Loading sparse matrix...")
            data = f['X']['data'][:]
            indices = f['X']['indices'][:]
            indptr = f['X']['indptr'][:]
            X_sparse = sparse.csr_matrix((data, indices, indptr),
                                         shape=(n_cells_total, n_genes_total))

            cell_idx_sorted = np.sort(cell_idx)
            gene_idx_sorted = np.sort(gene_idx)
            subset = X_sparse[cell_idx_sorted][:, gene_idx_sorted].toarray()

            cell_labels = obs_cell_labels[cell_idx_sorted]
            gene_symbols = var_gene_symbols[gene_idx_sorted]

            expr_df = pd.DataFrame(subset, index=cell_labels, columns=gene_symbols)
            expression_frames.append(expr_df)

            elapsed = time.time() - t0
            print(f"  Done in {elapsed:.1f}s")

            del X_sparse, data, indices, indptr, subset
            gc.collect()

    gc.collect()

    # Delete Isocortex h5ad files to free disk space
    if 'Isocortex' in file_path_str and os.path.exists(file_path_str):
        fsize = os.path.getsize(file_path_str) / 1e9
        os.remove(file_path_str)
        print(f"  Deleted {file_path_str} ({fsize:.1f} GB freed)")

# --- Combine and save ---
expression_data = pd.concat(expression_frames)
# Reorder columns to match receptor_genes order
expression_data = expression_data[receptor_genes]
expression_data.to_csv(output_csv)
print(f"\nSaved expression data: {expression_data.shape[0]:,} cells x {expression_data.shape[1]} genes")
print(f"Output: {output_csv}")
