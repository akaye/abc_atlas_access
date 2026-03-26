"""Extract ORB expression using h5py, commit+push CSVs after each file."""
import pandas as pd
import numpy as np
import h5py
import os
import gc
import json
import subprocess
import urllib.request
from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

download_base = Path('../../data/abc_atlas')
abc_cache = AbcProjectCache.from_s3_cache(download_base)

# Load metadata
cell = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata', dtype={'cell_label': str})
cell.set_index('cell_label', inplace=True)
gene = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')
gene.set_index('gene_identifier', inplace=True)
cluster_details = abc_cache.get_metadata_dataframe(
    directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_pivoted',
    keep_default_na=False
)
cluster_details.set_index('cluster_alias', inplace=True)
cell_extended = cell.join(cluster_details, on='cluster_alias')
orb_cells = cell_extended[cell_extended['region_of_interest_acronym'] == 'PL-ILA-ORB'].copy()
print(f"ORB cells: {len(orb_cells):,}")

# Save metadata immediately
outdir = Path('scube1_orb')
outdir.mkdir(exist_ok=True)
orb_cells.to_csv(outdir / 'orb_metadata.csv')
print(f"Saved orb_metadata.csv ({len(orb_cells):,} rows)")

# Genes
gene_symbols = ['Scube1', 'Drd1', 'Drd2', 'Foxp2', 'Chat', 'Pvalb', 'Sst', 'Th', 'Aqp4', 'Mog']
gene_ids = gene[gene['gene_symbol'].isin(gene_symbols)].index.tolist()
gene_sym_map = dict(zip(gene.loc[gene_ids].index, gene.loc[gene_ids]['gene_symbol']))

# Get isocortex file URLs
manifest_path = download_base / abc_cache.current_manifest
with open(manifest_path) as f:
    manifest = json.load(f)

iso_files = {}
for key, val in manifest['file_listing']['WMB-10Xv2']['expression_matrices'].items():
    if 'Isocortex' in key:
        h5ad_info = val.get('log2', {}).get('files', {}).get('h5ad', None)
        if h5ad_info:
            iso_files[key] = h5ad_info

print(f"Files to process: {sorted(iso_files.keys())}")

orb_cell_labels = set(orb_cells.index)
all_expression = []

def git_commit_push(message):
    """Stage, commit, and push scube1_orb directory."""
    subprocess.run(['git', 'add', 'notebooks/scube1_orb/'], cwd='/home/user/abc_atlas_access', check=True)
    subprocess.run(['git', 'commit', '-m', message], cwd='/home/user/abc_atlas_access', check=True)
    for attempt in range(4):
        result = subprocess.run(
            ['git', 'push', '-u', 'origin', 'claude/regenerate-scube1-plots-8Bdh0'],
            cwd='/home/user/abc_atlas_access', capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"  Pushed successfully")
            return
        import time
        time.sleep(2 ** (attempt + 1))
    print(f"  Push failed after retries: {result.stderr}")

# Commit metadata first
git_commit_push("Add ORB cell metadata CSV for Scube1 dotplot")

for file_key in sorted(iso_files.keys()):
    info = iso_files[file_key]
    local_path = download_base / info['relative_path']

    if not local_path.exists():
        print(f"\nDownloading {file_key}...")
        local_path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(info['url'], local_path)

    print(f"\nProcessing {file_key} ({local_path.stat().st_size / 1e9:.1f} GB)...")

    with h5py.File(local_path, 'r') as f:
        obs_names = f['obs']['cell_label'][:].astype(str)
        var_names = f['var']['gene_identifier'][:].astype(str)

        cell_mask = np.isin(obs_names, list(orb_cell_labels))
        n_match = cell_mask.sum()
        print(f"  {n_match:,} ORB cells in this file")

        if n_match > 0:
            cell_indices = np.where(cell_mask)[0]
            matching_labels = obs_names[cell_indices]

            valid_gene_indices = []
            valid_gene_names = []
            for gid in gene_ids:
                idx = np.where(var_names == gid)[0]
                if len(idx) > 0:
                    valid_gene_indices.append(idx[0])
                    valid_gene_names.append(gene_sym_map[gid])

            expr_data = np.zeros((n_match, len(valid_gene_indices)), dtype=np.float32)

            X = f['X']
            if isinstance(X, h5py.Group):
                # CSR sparse matrix
                data = X['data']
                indices = X['indices']
                indptr = X['indptr']
                for row_idx, ci in enumerate(cell_indices):
                    start, end = indptr[ci], indptr[ci + 1]
                    row_indices = indices[start:end]
                    row_data = data[start:end]
                    for j, gi in enumerate(valid_gene_indices):
                        pos = np.searchsorted(row_indices, gi)
                        if pos < len(row_indices) and row_indices[pos] == gi:
                            expr_data[row_idx, j] = row_data[pos]
            else:
                # Dense
                for j, gi in enumerate(valid_gene_indices):
                    expr_data[:, j] = X[cell_indices, gi]

            df = pd.DataFrame(expr_data, index=matching_labels, columns=valid_gene_names)
            all_expression.append(df)
            print(f"  Extracted {df.shape}")

            # Save incremental CSV and commit
            partial = pd.concat(all_expression)
            partial.to_csv(outdir / 'orb_expression.csv')
            print(f"  Saved incremental CSV ({partial.shape[0]:,} cells total)")
            git_commit_push(f"Update ORB expression CSV: add {n_match:,} cells from {file_key} ({partial.shape[0]:,} total)")

    # Delete h5ad to free space
    print(f"  Deleting {local_path.name}...")
    os.remove(local_path)
    gc.collect()

print(f"\nDone! Final expression: {pd.concat(all_expression).shape}")
