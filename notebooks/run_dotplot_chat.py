#!/usr/bin/env python3
"""Generate striatum dotplots with Chat+ manual-select group."""

import pandas as pd
import numpy as np
import re
import anndata
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
from abc_atlas_access.abc_atlas_cache.anndata_utils import get_gene_data

# ---------- 1. Load data ----------
download_base = Path('../../data/abc_atlas')
abc_cache = AbcProjectCache.from_s3_cache(download_base)

cell = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata', dtype={'cell_label': str})
cell.set_index('cell_label', inplace=True)
print(f"Total cells: {len(cell):,}")

gene = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')
gene.set_index('gene_identifier', inplace=True)

cluster_details = abc_cache.get_metadata_dataframe(
    directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_pivoted',
    keep_default_na=False
)
cluster_details.set_index('cluster_alias', inplace=True)
cell_extended = cell.join(cluster_details, on='cluster_alias')

# ---------- 2. Filter to striatum ----------
str_cells = cell_extended[cell_extended['region_of_interest_acronym'].isin(['STRd', 'STRv'])].copy()

cortex_patterns = ['CTX', 'L2/3', 'L4/5', 'L5 ', 'L6 ', 'L6b', 'RSP']
cortical = [s for s in str_cells['subclass'].unique() if any(p in s for p in cortex_patterns)]
str_cells = str_cells[~str_cells['subclass'].isin(cortical)].copy()
print(f"Striatum cells (after cortex exclusion): {len(str_cells):,}")

# ---------- 3. Gene list ----------
primary_gene = ['Scube1']
msn_markers = ['Drd1', 'Drd2', 'Foxp2']
interneuron_markers = ['Chat', 'Pvalb', 'Sst', 'Th']
glial_markers = ['Aqp4', 'Mog']
all_genes = primary_gene + msn_markers + interneuron_markers + glial_markers

found_symbols = set(gene[gene['gene_symbol'].isin(all_genes)]['gene_symbol'])
gene_list = [g for g in all_genes if g in found_symbols]
print(f"Genes: {gene_list}")

# ---------- 4. Load expression ----------
print("Loading expression data (this may take a few minutes)...")
expression_data = get_gene_data(
    abc_atlas_cache=abc_cache,
    all_cells=str_cells,
    all_genes=gene,
    selected_genes=gene_list,
    data_type='log2'
)
expression_data = expression_data.dropna(how='all')
print(f"Expression: {expression_data.shape[0]:,} cells x {expression_data.shape[1]} genes")

# ---------- 5. Build AnnData ----------
common_cells = str_cells.index.intersection(expression_data.index)
expression_data = expression_data.loc[common_cells, gene_list]
str_cells_aligned = str_cells.loc[common_cells]

region_map = {'STRd': 'Dorsal striatum', 'STRv': 'Ventral striatum'}
str_cells_aligned['region_label'] = str_cells_aligned['region_of_interest_acronym'].map(region_map)

adata = anndata.AnnData(
    X=expression_data.values.astype(np.float32),
    obs=str_cells_aligned[['subclass', 'supertype', 'class', 'region_of_interest_acronym', 'region_label']].copy(),
    var=pd.DataFrame(index=gene_list)
)
adata.obs['subclass_short'] = adata.obs['subclass'].apply(lambda x: re.sub(r'^\d+\s+', '', x))
adata.obs['subclass_short'] = pd.Categorical(adata.obs['subclass_short'])
adata.obs['subclass'] = pd.Categorical(adata.obs['subclass'])

# ---------- 6. Add Chat+ manual select group ----------
chat_col_idx = list(adata.var_names).index('Chat')
chat_expr = adata.X[:, chat_col_idx]
chat_mask = chat_expr > 0
n_chat_pos = int(chat_mask.sum())
print(f"\nChat+ cells (expression > 0): {n_chat_pos:,} out of {adata.n_obs:,}")

adata_chat = adata[chat_mask].copy()
adata_chat.obs['subclass_short'] = 'Chat+ (manual select)'
adata_chat.obs['subclass'] = 'Chat+ (manual select)'

adata_with_chat = anndata.concat([adata, adata_chat])

# Order: insert Chat+ right after PAL-STR Gaba-Chol
original_cats = list(adata.obs['subclass_short'].cat.categories)
insert_after = 'PAL-STR Gaba-Chol'
if insert_after in original_cats:
    idx = original_cats.index(insert_after) + 1
    ordered_cats = original_cats[:idx] + ['Chat+ (manual select)'] + original_cats[idx:]
else:
    ordered_cats = original_cats + ['Chat+ (manual select)']

adata_with_chat.obs['subclass_short'] = pd.Categorical(
    adata_with_chat.obs['subclass_short'], categories=ordered_cats, ordered=True
)

# Print cell counts
print(f"\n{'='*60}")
print(f"Cell counts per group (all striatum):")
print(f"{'='*60}")
group_counts = adata_with_chat.obs.groupby('subclass_short', observed=True).size()
for name in ordered_cats:
    if name in group_counts.index:
        print(f"  {name}: n={group_counts[name]:,}")
n_groups = len([c for c in ordered_cats if c in group_counts.index])
print(f"\nTotal groups: {n_groups}")

# ---------- 7. Dotplots ----------
gene_groups = {
    'Primary': [g for g in primary_gene if g in gene_list],
    'MSN markers': [g for g in msn_markers if g in gene_list],
    'Interneuron': [g for g in interneuron_markers if g in gene_list],
    'Glia': [g for g in glial_markers if g in gene_list],
}

outdir = Path('scube1')
outdir.mkdir(exist_ok=True)

# All striatum
dp = sc.pl.dotplot(
    adata_with_chat, var_names=gene_groups, groupby='subclass_short',
    standard_scale='var', cmap='Reds',
    figsize=(12, max(6, n_groups * 0.5)),
    show=False, return_fig=True
)
dp.style(dot_edge_color='black', dot_edge_lw=0.5)
dp.savefig(outdir / 'dotplot_striatum_Scube1_by_subclass.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved: {outdir / 'dotplot_striatum_Scube1_by_subclass.png'}")

# Per region
for region in ['STRd', 'STRv']:
    adata_region = adata_with_chat[adata_with_chat.obs['region_of_interest_acronym'] == region].copy()
    region_name = region_map[region]

    min_cells = 10
    sub_counts = adata_region.obs.groupby('subclass_short', observed=True).size()
    valid_subs = sub_counts[sub_counts >= min_cells].index
    adata_region = adata_region[adata_region.obs['subclass_short'].isin(valid_subs)].copy()

    region_cats = [c for c in ordered_cats if c in valid_subs]
    adata_region.obs['subclass_short'] = pd.Categorical(
        adata_region.obs['subclass_short'], categories=region_cats, ordered=True
    )
    n_sub = len(region_cats)

    print(f"\n{'='*60}")
    print(f"{region_name}: {adata_region.n_obs:,} cells, {n_sub} groups")
    print(f"{'='*60}")
    reg_counts = adata_region.obs.groupby('subclass_short', observed=True).size()
    for name in region_cats:
        if name in reg_counts.index:
            print(f"  {name}: n={reg_counts[name]:,}")

    dp = sc.pl.dotplot(
        adata_region, var_names=gene_groups, groupby='subclass_short',
        standard_scale='var', cmap='Reds',
        figsize=(12, max(5, n_sub * 0.5)),
        title=f'Scube1 & markers — {region_name}',
        show=False, return_fig=True
    )
    dp.style(dot_edge_color='black', dot_edge_lw=0.5)
    fname = outdir / f'dotplot_{region}_Scube1_by_subclass.png'
    dp.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fname}")

print("\nDone!")
