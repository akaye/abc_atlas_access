#!/usr/bin/env python3
"""Dotplot of Scube1 in thalamus with focus on medial & lateral habenula."""

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

# ---------- 2. Filter to thalamus ----------
th_cells = cell_extended[cell_extended['region_of_interest_acronym'] == 'TH'].copy()
print(f"Thalamus (TH) cells: {len(th_cells):,}")

# ---------- 3. Gene list ----------
primary_gene = ['Scube1']
habenula_markers = ['Tac1', 'Tac2', 'Chat', 'Pou4f1']
thalamus_markers = ['Slc17a6', 'Gad1', 'Prkcd']
glial_markers = ['Aqp4', 'Mog']

all_genes = primary_gene + habenula_markers + thalamus_markers + glial_markers
found_symbols = set(gene[gene['gene_symbol'].isin(all_genes)]['gene_symbol'])
gene_list = [g for g in all_genes if g in found_symbols]
missing = [g for g in all_genes if g not in found_symbols]
if missing:
    print(f"WARNING \u2014 genes not found: {missing}")
print(f"Genes: {gene_list}")

# ---------- 4. Load expression ----------
print("Loading expression data...")
expression_data = get_gene_data(
    abc_atlas_cache=abc_cache,
    all_cells=th_cells,
    all_genes=gene,
    selected_genes=gene_list,
    data_type='log2'
)
expression_data = expression_data.dropna(how='all')
print(f"Expression: {expression_data.shape[0]:,} cells x {expression_data.shape[1]} genes")

# ---------- 5. Build AnnData ----------
common_cells = th_cells.index.intersection(expression_data.index)
expression_data = expression_data.loc[common_cells, gene_list]
th_cells_aligned = th_cells.loc[common_cells]

adata = anndata.AnnData(
    X=expression_data.values.astype(np.float32),
    obs=th_cells_aligned[['subclass', 'supertype', 'class', 'region_of_interest_acronym']].copy(),
    var=pd.DataFrame(index=gene_list)
)
adata.obs['subclass_short'] = adata.obs['subclass'].apply(lambda x: re.sub(r'^\d+\s+', '', x))
adata.obs['supertype_short'] = adata.obs['supertype'].apply(lambda x: re.sub(r'^\d+\s+', '', x))

# ---------- 6. Dotplot: All TH subclasses (\u226510 cells) ----------
gene_groups = {
    'Primary': [g for g in primary_gene if g in gene_list],
    'Habenula': [g for g in habenula_markers if g in gene_list],
    'Thalamus': [g for g in thalamus_markers if g in gene_list],
    'Glia': [g for g in glial_markers if g in gene_list],
}

outdir = Path('scube1')
outdir.mkdir(exist_ok=True)

min_cells = 10
sub_counts = adata.obs.groupby('subclass_short', observed=True).size()
valid_subs = sub_counts[sub_counts >= min_cells].index
adata_filt = adata[adata.obs['subclass_short'].isin(valid_subs)].copy()
adata_filt.obs['subclass_short'] = pd.Categorical(adata_filt.obs['subclass_short'])
n_sub = adata_filt.obs['subclass_short'].nunique()

print(f"\n{'='*60}")
print(f"Thalamus subclasses (n\u2265{min_cells}): {n_sub}")
print(f"{'='*60}")
filt_counts = adata_filt.obs.groupby('subclass_short', observed=True).size().sort_values(ascending=False)
for name, n in filt_counts.items():
    print(f"  {name}: n={n:,}")

dp = sc.pl.dotplot(
    adata_filt, var_names=gene_groups, groupby='subclass_short',
    standard_scale='var', cmap='Reds',
    figsize=(12, max(6, n_sub * 0.5)),
    title='Scube1 & markers \u2014 Thalamus (all subclasses)',
    show=False, return_fig=True
)
dp.style(dot_edge_color='black', dot_edge_lw=0.5)
dp.savefig(outdir / 'dotplot_TH_Scube1_by_subclass.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {outdir / 'dotplot_TH_Scube1_by_subclass.png'}")

# ---------- 7. Dotplot: Habenula supertypes only ----------
hab_subclasses = ['MH Tac2 Glut', 'LH Pou4f1 Sox1 Glut']
adata_hab = adata[adata.obs['subclass_short'].isin(hab_subclasses)].copy()
adata_hab.obs['supertype_short'] = pd.Categorical(adata_hab.obs['supertype_short'])
n_st = adata_hab.obs['supertype_short'].nunique()

print(f"\n{'='*60}")
print(f"Habenula supertypes: {n_st}")
print(f"{'='*60}")
st_counts = adata_hab.obs.groupby('supertype_short', observed=True).size().sort_values(ascending=False)
for name, n in st_counts.items():
    print(f"  {name}: n={n:,}")

dp = sc.pl.dotplot(
    adata_hab, var_names=gene_groups, groupby='supertype_short',
    standard_scale='var', cmap='Reds',
    figsize=(12, max(4, n_st * 0.6)),
    title='Scube1 & markers \u2014 Habenula supertypes',
    show=False, return_fig=True
)
dp.style(dot_edge_color='black', dot_edge_lw=0.5)
dp.savefig(outdir / 'dotplot_habenula_Scube1_by_supertype.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {outdir / 'dotplot_habenula_Scube1_by_supertype.png'}")

print("\nDone!")
