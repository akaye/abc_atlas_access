#!/usr/bin/env python3
"""
Regenerate the mPFC NE/5-HT/DA receptor dot plot, simplified to:
  - 5 cortical subclasses: L2/3 IT, L5 IT, L5 PT (= L5 ET), Pvalb, Sst
  - 5-HT2A and 5-HT2C plus all Gs-coupled and Gi-coupled monoamine receptors

Loads from existing pre-extracted CSVs (mpfc_receptor_expression.csv +
mpfc_10x_metadata.csv) so no data needs to be re-downloaded.
"""
import re
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

NOTEBOOK_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = NOTEBOOK_DIR.parent / 'outputs'
OUTPUT_DIR.mkdir(exist_ok=True)

EXPR_CSV = NOTEBOOK_DIR / 'mpfc_receptor_expression.csv'
META_CSV = NOTEBOOK_DIR / 'mpfc_10x_metadata.csv'

REGION = 'PL-ILA-ORB'
TARGET_SUBCLASSES = [
    'L2/3 IT CTX Glut',
    'L5 IT CTX Glut',
    'L5 ET CTX Glut',  # L5 ET = L5 PT (pyramidal tract)
    'Pvalb Gaba',
    'Sst Gaba',
]

# 5-HT2A + 5-HT2C explicitly requested (Gq); plus Gs- and Gi-coupled receptors
# across 5-HT, NE, and DA systems. Drops Gq adrenergics (Adra1*), Htr2b,
# and ionotropic Htr3a/b that aren't of interest here.
GENE_GROUPS = {
    '5-HT Gi': ['Htr1a', 'Htr1b', 'Htr1d', 'Htr1f', 'Htr5a', 'Htr5b'],
    '5-HT Gq': ['Htr2a', 'Htr2c'],
    '5-HT Gs': ['Htr4', 'Htr6', 'Htr7'],
    'NE Gi':   ['Adra2a', 'Adra2b', 'Adra2c'],
    'NE Gs':   ['Adrb1', 'Adrb2', 'Adrb3'],
    'DA Gs':   ['Drd1', 'Drd5'],
    'DA Gi':   ['Drd2', 'Drd3', 'Drd4'],
}


def main():
    print(f'Loading {EXPR_CSV.name} ...')
    expr = pd.read_csv(EXPR_CSV, index_col=0)
    print(f'  {expr.shape[0]:,} cells x {expr.shape[1]} genes')

    print(f'Loading {META_CSV.name} ...')
    meta = pd.read_csv(META_CSV, index_col=0)
    meta = meta.loc[meta.index.intersection(expr.index)].copy()
    meta['subclass_short'] = meta['subclass'].apply(
        lambda x: re.sub(r'^\d+\s+', '', str(x))
    )

    mask = (
        (meta['region_of_interest_acronym'] == REGION)
        & (meta['subclass_short'].isin(TARGET_SUBCLASSES))
    )
    meta_sel = meta.loc[mask]
    expr_sel = expr.loc[meta_sel.index]
    print(f'\nFiltered to {REGION} + 5 subclasses: {len(meta_sel):,} cells')
    for sc_name in TARGET_SUBCLASSES:
        n = (meta_sel['subclass_short'] == sc_name).sum()
        print(f'  {sc_name}: {n:,}')

    # Flatten + verify gene availability
    available = set(expr_sel.columns)
    gene_groups = {
        grp: [g for g in genes if g in available]
        for grp, genes in GENE_GROUPS.items()
    }
    gene_groups = {k: v for k, v in gene_groups.items() if v}
    flat_genes = [g for genes in gene_groups.values() for g in genes]
    print(f'\nReceptors plotted ({len(flat_genes)}):')
    for grp, genes in gene_groups.items():
        print(f'  {grp}: {genes}')

    expr_sel = expr_sel[flat_genes]

    # Build AnnData with categorical subclass in requested order
    adata = anndata.AnnData(
        X=expr_sel.values,
        obs=meta_sel[['subclass_short']].copy(),
        var=pd.DataFrame(index=flat_genes),
    )
    adata.obs['subclass_short'] = pd.Categorical(
        adata.obs['subclass_short'],
        categories=TARGET_SUBCLASSES,
        ordered=True,
    )

    n_sc = len(TARGET_SUBCLASSES)
    fig_w = max(8, 0.55 * len(flat_genes) + 3)
    fig_h = max(3.5, 0.55 * n_sc + 2)

    dp = sc.pl.dotplot(
        adata,
        var_names=gene_groups,
        groupby='subclass_short',
        standard_scale='var',
        cmap='Reds',
        figsize=(fig_w, fig_h),
        show=False,
        return_fig=True,
    )
    dp.style(dot_edge_color='black', dot_edge_lw=0.5)

    out = OUTPUT_DIR / 'dotplot_mPFC_receptors_simplified.png'
    dp.savefig(out, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f'\nSaved: {out}')


if __name__ == '__main__':
    main()
