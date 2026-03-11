#!/usr/bin/env python3
"""
Extract neuromodulator receptor expression for locus coeruleus & peri-LC
cell types from ABC Atlas 10x snRNA-seq data.

Target region: Pontine dissection regions (P) that contain LC cell types.
Key cell type: NTS Dbh Glut (noradrenergic neurons of the locus coeruleus).

Gene panel: 28 neuromodulator receptors
  - Serotonin (14): Htr1a-Htr7
  - Norepinephrine (9): Adra1a-Adrb3
  - Dopamine (5): Drd1-Drd5
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
output_csv = Path('lc_10x_expression.csv')
metadata_csv = Path('lc_10x_metadata.csv')

# Receptor genes (full panel available in 10x)
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

# Pontine dissection regions that contain LC and peri-LC cell types
pons_rois = ['P']

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

# --- Identify pontine cells ---
pons_cells = cell_extended[
    cell_extended['region_of_interest_acronym'].isin(pons_rois)
].copy()
print(f"\nPontine (P) cells: {len(pons_cells):,}")

print(f"\nBy class:")
for cls, cnt in pons_cells.groupby('class', observed=True).size().sort_values(ascending=False).items():
    print(f"  {cls}: {cnt:,}")

print(f"\nBy subclass (>= 10 cells):")
sc_counts = pons_cells.groupby('subclass', observed=True).size().sort_values(ascending=False)
for sc_name, cnt in sc_counts.items():
    if cnt >= 10:
        print(f"  {sc_name}: {cnt:,}")

# Check for Dbh+ (noradrenergic) neurons
dbh_mask = pons_cells['subclass'].str.contains('Dbh', na=False)
print(f"\nDbh+ (noradrenergic) neurons: {dbh_mask.sum():,}")

# Save metadata
meta_cols = ['subclass', 'supertype', 'class', 'neurotransmitter',
             'region_of_interest_acronym', 'feature_matrix_label', 'dataset_label']
pons_cells[meta_cols].to_csv(metadata_csv)
print(f"Saved metadata: {metadata_csv}")

# --- Identify receptor genes ---
available_genes = gene[gene['gene_symbol'].isin(all_receptors)]
receptor_genes = [g for g in all_receptors if g in set(available_genes['gene_symbol'])]
gene_ensembl_ids = set(available_genes.index.tolist())
print(f"\nReceptor genes found: {len(receptor_genes)} / {len(all_receptors)}")

# --- Extract expression data ---
pons_matrices = pons_cells.groupby('feature_matrix_label').size()
print(f"\nExpression matrices to process: {len(pons_matrices)}")
for mat, count in pons_matrices.items():
    print(f"  {mat}: {count:,} cells")

expression_frames = []

for i, matrix_label in enumerate(pons_matrices.index):
    dataset_label = pons_cells[
        pons_cells['feature_matrix_label'] == matrix_label
    ]['dataset_label'].iloc[0]
    file_name = f"{matrix_label}/log2"

    print(f"\n[{i+1}/{len(pons_matrices)}] Processing {file_name}...")
    t0 = time.time()

    file_path = abc_cache.get_file_path(directory=dataset_label, file_name=file_name)
    file_path_str = str(file_path)

    with h5py.File(file_path_str, 'r') as f:
        obs_cell_labels = f['obs']['cell_label'][:].astype(str)
        var_gene_ids = f['var']['gene_identifier'][:].astype(str)
        var_gene_symbols = f['var']['gene_symbol'][:].astype(str)

        n_cells_total = len(obs_cell_labels)
        n_genes_total = len(var_gene_ids)

        pons_cell_label_set = set(
            pons_cells[pons_cells['feature_matrix_label'] == matrix_label].index
        )
        cell_idx = np.array([j for j, name in enumerate(obs_cell_labels)
                             if name in pons_cell_label_set])
        gene_idx = np.array([j for j, gid in enumerate(var_gene_ids)
                             if gid in gene_ensembl_ids])

        print(f"  Total cells in matrix: {n_cells_total:,}")
        print(f"  Pons cells found: {len(cell_idx):,}")
        print(f"  Receptor genes found: {len(gene_idx)}")

        if len(cell_idx) == 0:
            print(f"  Skipping (no matching cells)")
        else:
            # Memory-efficient: read indptr to extract rows one chunk at a time
            print(f"  Reading sparse matrix row-by-row (memory-efficient)...")
            indptr = f['X']['indptr'][:]
            gene_idx_sorted = np.sort(gene_idx)
            gene_idx_set = set(gene_idx_sorted.tolist())
            cell_idx_sorted = np.sort(cell_idx)

            # Pre-allocate output array
            n_target = len(cell_idx_sorted)
            n_genes_target = len(gene_idx_sorted)
            result = np.zeros((n_target, n_genes_target), dtype=np.float32)

            # Map gene indices to output columns
            gene_to_col = {g: c for c, g in enumerate(gene_idx_sorted)}

            # Process in chunks
            chunk_size = 5000
            for chunk_start in range(0, n_target, chunk_size):
                chunk_end = min(chunk_start + chunk_size, n_target)
                chunk_cells = cell_idx_sorted[chunk_start:chunk_end]

                for local_i, global_i in enumerate(chunk_cells):
                    row_start = indptr[global_i]
                    row_end = indptr[global_i + 1]
                    if row_start == row_end:
                        continue
                    row_indices = f['X']['indices'][row_start:row_end]
                    row_data = f['X']['data'][row_start:row_end]
                    for idx, val in zip(row_indices, row_data):
                        if idx in gene_to_col:
                            result[chunk_start + local_i, gene_to_col[idx]] = val

                if (chunk_end % 10000 == 0) or chunk_end == n_target:
                    print(f"    Processed {chunk_end:,}/{n_target:,} cells...")

            cell_labels = obs_cell_labels[cell_idx_sorted]
            gene_symbols = var_gene_symbols[gene_idx_sorted]

            expr_df = pd.DataFrame(result, index=cell_labels, columns=gene_symbols)
            expression_frames.append(expr_df)

            elapsed = time.time() - t0
            print(f"  Done in {elapsed:.1f}s")

            del result
            gc.collect()

    gc.collect()

    # Delete h5ad to free disk space
    if os.path.exists(file_path_str):
        fsize = os.path.getsize(file_path_str) / 1e9
        os.remove(file_path_str)
        print(f"  Deleted {file_path_str} ({fsize:.1f} GB freed)")

# --- Combine and save ---
expression_data = pd.concat(expression_frames)
expression_data = expression_data[receptor_genes]
expression_data.to_csv(output_csv)
print(f"\nSaved expression data: {expression_data.shape[0]:,} cells x {expression_data.shape[1]} genes")
print(f"Output: {output_csv}")
