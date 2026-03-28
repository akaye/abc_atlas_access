"""Extract TH expression using h5py, save CSV, commit+push after each file."""
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

th_cells = cell_extended[cell_extended['region_of_interest_acronym'] == 'TH'].copy()
print(f"TH cells: {len(th_cells):,}")

outdir = Path('scube1_th')
outdir.mkdir(exist_ok=True)
th_cells.to_csv(outdir / 'th_metadata.csv')
print("Saved th_metadata.csv")

gene_symbols = ['Scube1', 'Drd1', 'Drd2', 'Foxp2', 'Chat', 'Pvalb', 'Sst', 'Th', 'Aqp4', 'Mog']
gene_ids = gene[gene['gene_symbol'].isin(gene_symbols)].index.tolist()
gene_sym_map = dict(zip(gene.loc[gene_ids].index, gene.loc[gene_ids]['gene_symbol']))

manifest_path = download_base / abc_cache.current_manifest
with open(manifest_path) as f:
    manifest = json.load(f)

def git_commit_push(message):
    subprocess.run(['git', 'add', 'notebooks/scube1_th/'], cwd='/home/user/abc_atlas_access', check=True)
    r = subprocess.run(['git', 'commit', '-m', message], cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  Nothing to commit (already up to date)")
        return
    for attempt in range(4):
        r = subprocess.run(['git', 'push', '-u', 'origin', 'claude/regenerate-scube1-plots-8Bdh0'],
                           cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
        if r.returncode == 0: print("  Pushed"); return
        import time; time.sleep(2 ** (attempt + 1))
    print(f"  Push failed: {r.stderr}")

git_commit_push("Add TH cell metadata CSV for Scube1 dotplot")

th_cell_labels = set(th_cells.index)
all_expression = []

for ver in ['WMB-10Xv2', 'WMB-10Xv3']:
    key = f'{ver}-TH'
    if key not in manifest['file_listing'][ver].get('expression_matrices', {}):
        continue
    h5ad_info = manifest['file_listing'][ver]['expression_matrices'][key]['log2']['files']['h5ad']
    local_path = download_base / h5ad_info['relative_path']

    if not local_path.exists():
        print(f"\nDownloading {key}...")
        local_path.parent.mkdir(parents=True, exist_ok=True)
        import time as _time
        for dl_attempt in range(5):
            try:
                urllib.request.urlretrieve(h5ad_info['url'], local_path)
                break
            except Exception as e:
                print(f"  Download attempt {dl_attempt+1} failed: {e}")
                if local_path.exists(): os.remove(local_path)
                _time.sleep(2 ** (dl_attempt + 1))
        else:
            raise RuntimeError(f"Failed to download {key} after 5 attempts")

    print(f"\nProcessing {key} ({local_path.stat().st_size / 1e9:.1f} GB)...")

    with h5py.File(local_path, 'r') as f:
        obs_names = f['obs']['cell_label'][:].astype(str)
        var_names = f['var']['gene_identifier'][:].astype(str)
        cell_mask = np.isin(obs_names, list(th_cell_labels))
        n_match = cell_mask.sum()
        print(f"  {n_match:,} TH cells in file")

        if n_match > 0:
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
            all_expression.append(df)
            print(f"  Extracted {df.shape}")

            partial = pd.concat(all_expression)
            partial.to_csv(outdir / 'th_expression.csv')
            print(f"  Saved incremental CSV ({partial.shape[0]:,} cells total)")
            git_commit_push(f"Update TH expression CSV: add {n_match:,} cells from {key} ({partial.shape[0]:,} total)")

    print(f"  Deleting {local_path.name}...")
    os.remove(local_path)
    gc.collect()

print(f"\nDone! Final: {pd.concat(all_expression).shape}")
