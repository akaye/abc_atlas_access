#!/usr/bin/env python3
"""
Extract phosphodiesterase (PDE) gene expression data for mPFC cells from ABC Atlas.

Same approach as extract_mpfc_expression.py but for PDE genes instead of
neuromodulator receptors. Processes the large Isocortex h5ad files one at
a time, extracting 17 neurally-relevant PDE genes for mPFC cells.
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
output_csv = Path('mpfc_pde_expression.csv')
metadata_csv = Path('mpfc_pde_metadata.csv')

# Phosphodiesterase genes (neurally relevant)
pde1_family = ['Pde1a', 'Pde1b', 'Pde1c']       # Ca2+/CaM-stimulated, cAMP+cGMP
pde2_family = ['Pde2a']                           # cGMP-stimulated, cAMP+cGMP
pde3_family = ['Pde3a', 'Pde3b']                  # cGMP-inhibited, cAMP
pde4_family = ['Pde4a', 'Pde4b', 'Pde4c', 'Pde4d']  # cAMP-specific (major targets)
pde5_family = ['Pde5a']                            # cGMP-specific
pde7_family = ['Pde7a', 'Pde7b']                  # cAMP-specific
pde8_family = ['Pde8a', 'Pde8b']                  # cAMP-specific
pde9_family = ['Pde9a']                            # cGMP-specific
pde10_family = ['Pde10a']                          # dual cAMP+cGMP (striatal)
pde11_family = ['Pde11a']                          # dual cAMP+cGMP (hippocampal)

all_pdes = (pde1_family + pde2_family + pde3_family + pde4_family +
            pde5_family + pde7_family + pde8_family + pde9_family +
            pde10_family + pde11_family)

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

# Save metadata
meta_cols = ['subclass', 'supertype', 'class', 'neurotransmitter',
             'region_of_interest_acronym', 'feature_matrix_label', 'dataset_label']
mpfc_selected[meta_cols].to_csv(metadata_csv)
print(f"Saved metadata: {metadata_csv}")

# --- Identify PDE genes ---
available_genes = gene[gene['gene_symbol'].isin(all_pdes)]
pde_genes = [g for g in all_pdes if g in set(available_genes['gene_symbol'])]
gene_ensembl_ids = set(available_genes.index.tolist())
print(f"PDE genes found: {len(pde_genes)} — {pde_genes}")

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

    with h5py.File(file_path_str, 'r') as f:
        obs_cell_labels = f['obs']['cell_label'][:].astype(str)
        var_gene_ids = f['var']['gene_identifier'][:].astype(str)
        var_gene_symbols = f['var']['gene_symbol'][:].astype(str)

        n_cells_total = len(obs_cell_labels)
        n_genes_total = len(var_gene_ids)

        mpfc_cell_label_set = set(
            mpfc_selected[mpfc_selected['feature_matrix_label'] == matrix_label].index
        )
        cell_idx = np.array([j for j, name in enumerate(obs_cell_labels)
                             if name in mpfc_cell_label_set])
        gene_idx = np.array([j for j, gid in enumerate(var_gene_ids)
                             if gid in gene_ensembl_ids])

        print(f"  Total cells in matrix: {n_cells_total:,}")
        print(f"  mPFC cells found: {len(cell_idx):,}")
        print(f"  PDE genes found: {len(gene_idx)}")

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
expression_data = expression_data[pde_genes]
expression_data.to_csv(output_csv)
print(f"\nSaved expression data: {expression_data.shape[0]:,} cells x {expression_data.shape[1]} genes")
print(f"Output: {output_csv}")
