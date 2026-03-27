"""Extract STR expression using h5py, save CSV, commit+push."""
import pandas as pd, numpy as np, h5py, os, gc, json, subprocess, urllib.request, re
from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

download_base = Path('../../data/abc_atlas')
abc_cache = AbcProjectCache.from_s3_cache(download_base)

cell = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata', dtype={'cell_label': str})
cell.set_index('cell_label', inplace=True)
gene = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')
gene.set_index('gene_identifier', inplace=True)
cluster_details = abc_cache.get_metadata_dataframe(directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_pivoted', keep_default_na=False)
cluster_details.set_index('cluster_alias', inplace=True)
cell_extended = cell.join(cluster_details, on='cluster_alias')

str_cells = cell_extended[cell_extended['region_of_interest_acronym'].isin(['STRd', 'STRv'])].copy()
print(f"STR cells: {len(str_cells):,}")

outdir = Path('scube1')
outdir.mkdir(exist_ok=True)
str_cells.to_csv(outdir / 'str_metadata.csv')
print("Saved str_metadata.csv")

gene_symbols = ['Scube1', 'Drd1', 'Drd2', 'Foxp2', 'Chat', 'Pvalb', 'Sst', 'Th', 'Aqp4', 'Mog']
gene_ids = gene[gene['gene_symbol'].isin(gene_symbols)].index.tolist()
gene_sym_map = dict(zip(gene.loc[gene_ids].index, gene.loc[gene_ids]['gene_symbol']))

manifest_path = download_base / abc_cache.current_manifest
with open(manifest_path) as f:
    manifest = json.load(f)

h5ad_info = manifest['file_listing']['WMB-10Xv3']['expression_matrices']['WMB-10Xv3-STR']['log2']['files']['h5ad']
local_path = download_base / h5ad_info['relative_path']

if not local_path.exists():
    print(f"Downloading WMB-10Xv3-STR...")
    local_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(h5ad_info['url'], local_path)

print(f"Processing WMB-10Xv3-STR ({local_path.stat().st_size / 1e9:.1f} GB)...")
str_cell_labels = set(str_cells.index)

with h5py.File(local_path, 'r') as f:
    obs_names = f['obs']['cell_label'][:].astype(str)
    var_names = f['var']['gene_identifier'][:].astype(str)
    cell_mask = np.isin(obs_names, list(str_cell_labels))
    n_match = cell_mask.sum()
    print(f"  {n_match:,} STR cells in file")

    cell_indices = np.where(cell_mask)[0]
    matching_labels = obs_names[cell_indices]
    valid_gene_indices, valid_gene_names = [], []
    for gid in gene_ids:
        idx = np.where(var_names == gid)[0]
        if len(idx) > 0:
            valid_gene_indices.append(idx[0])
            valid_gene_names.append(gene_sym_map[gid])

    expr_data = np.zeros((n_match, len(valid_gene_indices)), dtype=np.float32)
    X = f['X']
    data, indices, indptr = X['data'], X['indices'], X['indptr']
    for row_idx, ci in enumerate(cell_indices):
        start, end = indptr[ci], indptr[ci + 1]
        ri = indices[start:end]
        rd = data[start:end]
        for j, gi in enumerate(valid_gene_indices):
            pos = np.searchsorted(ri, gi)
            if pos < len(ri) and ri[pos] == gi:
                expr_data[row_idx, j] = rd[pos]

    df = pd.DataFrame(expr_data, index=matching_labels, columns=valid_gene_names)
    df.to_csv(outdir / 'str_expression.csv')
    print(f"  Saved str_expression.csv ({df.shape})")

print("Deleting h5ad...")
os.remove(local_path)

# Commit and push
subprocess.run(['git', 'add', 'notebooks/scube1/str_expression.csv', 'notebooks/scube1/str_metadata.csv'],
               cwd='/home/user/abc_atlas_access', check=True)
subprocess.run(['git', 'commit', '-m', 'Add STR expression and metadata CSVs for Scube1 dotplot'],
               cwd='/home/user/abc_atlas_access', check=True)
for attempt in range(4):
    r = subprocess.run(['git', 'push', '-u', 'origin', 'claude/regenerate-scube1-plots-8Bdh0'],
                       cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
    if r.returncode == 0:
        print("Pushed successfully")
        break
    import time; time.sleep(2 ** (attempt + 1))

print("Done!")
