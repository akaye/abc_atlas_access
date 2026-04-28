#!/usr/bin/env python3
"""
Regenerate the mPFC NE/5-HT/DA receptor dot plot, simplified to:
  - 5 cortical subclasses: L2/3 IT, L5 IT, L5 PT (= L5 ET), Pvalb, Sst
  - 5-HT2A and 5-HT2C plus all Gs-coupled and Gi-coupled monoamine receptors

Plots raw log2(CPM+1) expression (the WMB-10X "log2" matrices already store
log2 CPM+1), drops receptors that are essentially silent across all 5
subclasses (max fraction expressing < 5%), and uses 3x-larger fonts.

Loads from existing pre-extracted CSVs (mpfc_receptor_expression.csv +
mpfc_10x_metadata.csv) so no data needs to be re-downloaded.
"""
import re
from pathlib import Path

import anndata
import matplotlib as mpl
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

# Drop receptors expressed in fewer than this fraction of cells in EVERY
# selected subclass (i.e. essentially absent from this circuit).
MIN_MAX_FRAC = 0.05
FONT_SCALE = 3.0


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

    # Drop receptors that are essentially absent (max fraction expressing
    # across the 5 subclasses < MIN_MAX_FRAC).
    expr_with_sc = expr_sel.copy()
    expr_with_sc['_sc'] = meta_sel['subclass_short'].values
    frac = expr_with_sc.groupby('_sc', observed=True).apply(
        lambda x: (x.drop(columns='_sc') > 0).mean()
    )
    max_frac = frac.max(axis=0)
    dropped = sorted(g for g, v in max_frac.items() if v < MIN_MAX_FRAC)
    if dropped:
        print(f'\nDropped (max fraction expressing < {MIN_MAX_FRAC:.0%}):')
        for g in dropped:
            print(f'  {g}: {max_frac[g]:.3f}')

    gene_groups = {
        grp: [g for g in genes if g not in dropped]
        for grp, genes in gene_groups.items()
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

    # Bump every font 3x via matplotlib rcParams (default font.size = 10).
    base = 10.0
    new_size = base * FONT_SCALE
    mpl.rcParams.update({
        'font.size': new_size,
        'axes.titlesize': new_size,
        'axes.labelsize': new_size,
        'xtick.labelsize': new_size,
        'ytick.labelsize': new_size,
        'legend.fontsize': new_size,
        'legend.title_fontsize': new_size,
    })

    n_sc = len(TARGET_SUBCLASSES)
    # Scale figure modestly so the larger fonts have room to breathe but
    # don't dwarf the dots; scale `largest_dot` to keep dots visually
    # proportional to the bigger labels.
    fig_w = (0.85 * len(flat_genes) + 6) * 1.4
    fig_h = (0.85 * n_sc + 3) * 1.4

    # No standard_scale: plot raw log2(CPM+1) so absolute expression is
    # comparable across genes and subclasses.
    dp = sc.pl.dotplot(
        adata,
        var_names=gene_groups,
        groupby='subclass_short',
        cmap='Reds',
        figsize=(fig_w, fig_h),
        dot_min=0.25,
        dot_max=0.75,
        colorbar_title='Mean\nlog2(CPM+1)',
        size_title='% expressing',
        show=False,
        return_fig=True,
    )
    dp.style(dot_edge_color='black', dot_edge_lw=0.5,
             largest_dot=200 * (FONT_SCALE ** 1.5))
    # Give the legends room to breathe at 3x font size.
    dp.legends_width = 2.5 * FONT_SCALE

    out = OUTPUT_DIR / 'dotplot_mPFC_receptors_simplified.png'
    dp.savefig(out, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f'\nSaved: {out}')


if __name__ == '__main__':
    main()
