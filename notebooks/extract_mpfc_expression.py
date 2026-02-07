#!/usr/bin/env python3
"""
Extract receptor gene expression data for mPFC cells from ABC Atlas.

This script processes the large Isocortex h5ad files one at a time,
extracting only the 28 receptor genes for mPFC cells, then saves the
result as a small CSV file for the notebook to load.

It handles disk space constraints by deleting each ~8 GB h5ad file
after extracting the tiny subset (~165k cells × 28 genes).
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
output_csv = Path('mpfc_receptor_expression.csv')

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

# mPFC dissection regions
mpfc_rois = ['PL-ILA-ORB', 'ACA']

# Excitatory subclasses
excitatory_subclasses = [
    '007 L2/3 IT CTX Glut', '006 L4/5 IT CTX Glut', '005 L5 IT CTX Glut',
    '022 L5 ET CTX Glut', '032 L5 NP CTX Glut', '004 L6 IT CTX Glut',
    '030 L6 CT CTX Glut', '029 L6b CTX Glut', '021 L4 RSP-ACA Glut',
    '020 L2/3 IT RSP Glut',
]
# Interneuron subclasses
interneuron_subclasses = [
    '052 Pvalb Gaba', '051 Pvalb chandelier Gaba',
    '053 Sst Gaba', '056 Sst Chodl Gaba',
    '046 Vip Gaba', '049 Lamp5 Gaba',
    '050 Lamp5 Lhx6 Gaba', '047 Sncg Gaba',
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

# --- Identify mPFC cells ---
mpfc_cells = cell_extended[cell_extended['region_of_interest_acronym'].isin(mpfc_rois)]
mpfc_subclass_set = set(mpfc_cells['subclass'].unique())
excitatory_subclasses = [s for s in excitatory_subclasses if s in mpfc_subclass_set]
interneuron_subclasses = [s for s in interneuron_subclasses if s in mpfc_subclass_set]
selected_subclasses = sorted(set(excitatory_subclasses + interneuron_subclasses))
mpfc_selected = mpfc_cells[mpfc_cells['subclass'].isin(selected_subclasses)].copy()
print(f"mPFC selected cells: {len(mpfc_selected):,}")

# --- Identify receptor genes ---
available_genes = gene[gene['gene_symbol'].isin(all_receptors)]
receptor_genes = [g for g in all_receptors if g in set(available_genes['gene_symbol'])]
gene_ensembl_ids = set(available_genes.index.tolist())
# Map ensembl -> symbol for column naming
ensembl_to_symbol = dict(zip(available_genes.index, available_genes['gene_symbol']))
print(f"Receptor genes found: {len(receptor_genes)}")

# --- Extract expression data ---
mpfc_matrices = mpfc_selected.groupby('feature_matrix_label').size()
print(f"\nExpression matrices to process: {len(mpfc_matrices)}")
for mat, count in mpfc_matrices.items():
    print(f"  {mat}: {count:,} cells")

expression_frames = []

for i, matrix_label in enumerate(mpfc_matrices.index):
    dataset_label = mpfc_selected[
        mpfc_selected['feature_matrix_label'] == matrix_label
    ]['dataset_label'].iloc[0]
    file_name = f"{matrix_label}/log2"

    print(f"\n[{i+1}/{len(mpfc_matrices)}] Processing {file_name}...")
    t0 = time.time()

    file_path = abc_cache.get_file_path(directory=dataset_label, file_name=file_name)
    file_path_str = str(file_path)

    # Use h5py directly for much faster access than backed anndata
    with h5py.File(file_path_str, 'r') as f:
        # Read cell labels and gene identifiers
        obs_cell_labels = f['obs']['cell_label'][:].astype(str)
        var_gene_ids = f['var']['gene_identifier'][:].astype(str)
        var_gene_symbols = f['var']['gene_symbol'][:].astype(str)

        n_cells_total = len(obs_cell_labels)
        n_genes_total = len(var_gene_ids)

        # Find indices for our target cells and genes
        mpfc_cell_label_set = set(
            mpfc_selected[mpfc_selected['feature_matrix_label'] == matrix_label].index
        )
        cell_idx = np.array([j for j, name in enumerate(obs_cell_labels)
                             if name in mpfc_cell_label_set])
        gene_idx = np.array([j for j, gid in enumerate(var_gene_ids)
                             if gid in gene_ensembl_ids])

        print(f"  Total cells in matrix: {n_cells_total:,}")
        print(f"  mPFC cells found: {len(cell_idx):,}")
        print(f"  Receptor genes found: {len(gene_idx)}")

        if len(cell_idx) == 0:
            print(f"  Skipping (no matching cells)")
        else:
            # Load the sparse CSR matrix
            print(f"  Loading sparse matrix...")
            data = f['X']['data'][:]
            indices = f['X']['indices'][:]
            indptr = f['X']['indptr'][:]
            X_sparse = sparse.csr_matrix((data, indices, indptr),
                                         shape=(n_cells_total, n_genes_total))

            # Extract subset: selected cells × selected genes
            cell_idx_sorted = np.sort(cell_idx)
            gene_idx_sorted = np.sort(gene_idx)
            subset = X_sparse[cell_idx_sorted][:, gene_idx_sorted].toarray()

            # Build DataFrame
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
