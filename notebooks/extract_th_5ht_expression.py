#!/usr/bin/env python3
"""
Extract serotonin receptor + marker expression for thalamus (TH) from ABC Atlas 10x.

We extract the full TH region, then filter to lateral habenula (LH) cells
in the notebook.

Target genes:
  - Htr1a (5-HT1A), Htr2a (5-HT2A), Htr2c (5-HT2C): serotonin receptors
  - Slc17a7 (Vglut1): glutamatergic marker
  - Gad1, Gad2: GABAergic markers
  - Sst, Pvalb, Vip, Lamp5: interneuron subclass markers
  - Slc6a4 (SERT): serotonin transporter

Region: TH (thalamus — includes lateral habenula LH subclass)
"""
import os
import gc
import time
import h5py
import numpy as np
import pandas as pd
from pathlib import Path

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# --- Configuration ---
download_base = Path('../../data/abc_atlas')
region_acronym = 'TH'
output_csv = Path('th_5ht_expression.csv')
output_meta = Path('th_5ht_metadata.csv')

target_genes = [
    'Htr1a', 'Htr2a', 'Htr2c',
    'Slc17a7', 'Gad1', 'Gad2',
    'Sst', 'Pvalb', 'Vip', 'Lamp5',
    'Slc6a4',
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

# --- Identify region cells ---
region_cells = cell_extended[cell_extended['region_of_interest_acronym'] == region_acronym].copy()
print(f"\n{region_acronym} cells: {len(region_cells):,}")

region_cells[['subclass', 'supertype', 'class', 'neurotransmitter',
            'region_of_interest_acronym']].to_csv(output_meta)
print(f"Saved metadata: {output_meta}")

# --- Identify target genes ---
available_genes = gene[gene['gene_symbol'].isin(target_genes)]
found_genes = [g for g in target_genes if g in set(available_genes['gene_symbol'])]
gene_ensembl_ids = set(available_genes.index.tolist())
print(f"\nTarget genes found: {len(found_genes)} / {len(target_genes)}")
print(f"  {found_genes}")

# --- Extract expression data ---
region_matrices = region_cells.groupby('feature_matrix_label').size()
print(f"\nExpression matrices to process: {len(region_matrices)}")
for mat, count in region_matrices.items():
    print(f"  {mat}: {count:,} cells")

expression_frames = []

for i, matrix_label in enumerate(region_matrices.index):
    dataset_label = region_cells[
        region_cells['feature_matrix_label'] == matrix_label
    ]['dataset_label'].iloc[0]
    file_name = f"{matrix_label}/log2"

    print(f"\n[{i+1}/{len(region_matrices)}] Processing {file_name}...")
    t0 = time.time()

    file_path = abc_cache.get_file_path(directory=dataset_label, file_name=file_name)
    file_path_str = str(file_path)

    with h5py.File(file_path_str, 'r') as f:
        obs_cell_labels = f['obs']['cell_label'][:].astype(str)
        var_gene_ids = f['var']['gene_identifier'][:].astype(str)
        var_gene_symbols = f['var']['gene_symbol'][:].astype(str)

        n_cells_total = len(obs_cell_labels)

        region_cell_label_set = set(
            region_cells[region_cells['feature_matrix_label'] == matrix_label].index
        )
        cell_idx = np.array([j for j, name in enumerate(obs_cell_labels)
                             if name in region_cell_label_set])
        gene_idx = np.array([j for j, gid in enumerate(var_gene_ids)
                             if gid in gene_ensembl_ids])

        print(f"  Total cells in matrix: {n_cells_total:,}")
        print(f"  {region_acronym} cells found: {len(cell_idx):,}")
        print(f"  Target genes found: {len(gene_idx)}")

        if len(cell_idx) == 0:
            print(f"  Skipping (no matching cells)")
        else:
            print(f"  Reading sparse matrix row-by-row (memory-efficient)...")
            indptr = f['X']['indptr'][:]
            gene_idx_sorted = np.sort(gene_idx)
            cell_idx_sorted = np.sort(cell_idx)

            n_target = len(cell_idx_sorted)
            n_genes_target = len(gene_idx_sorted)
            result = np.zeros((n_target, n_genes_target), dtype=np.float32)

            gene_to_col = {g: c for c, g in enumerate(gene_idx_sorted)}

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

    if os.path.exists(file_path_str):
        fsize = os.path.getsize(file_path_str) / 1e9
        os.remove(file_path_str)
        print(f"  Deleted {file_path_str} ({fsize:.1f} GB freed)")

# --- Combine and save ---
expression_data = pd.concat(expression_frames)
expression_data = expression_data[found_genes]
expression_data.to_csv(output_csv)
print(f"\nSaved expression data: {expression_data.shape[0]:,} cells x {expression_data.shape[1]} genes")
print(f"Output: {output_csv}")
