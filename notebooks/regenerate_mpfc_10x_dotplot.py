#!/usr/bin/env python3
"""
Companion to regenerate_mpfc_simplified_dotplot.py.

Generates a broader mPFC NE/5-HT/DA receptor dot plot using the WMB-10X
single-cell data:
  - 10 major cell classes: L2/3 IT, L5 IT, L5 ET (=PT), L6 IT, L6 CT, L6b,
    Pvalb, Sst, Vip, Lamp5
  - All receptors profiled in the 10x panel that aren't completely silent
    across these classes -- so the "extra" receptors absent from MERFISH
    (Adra1a/b/d, Htr1d, Htr2b, Htr3a/b, Adrb2, Drd3/4/5, etc.) are included.

Loads from the existing pre-extracted CSVs (mpfc_receptor_expression.csv +
mpfc_10x_metadata.csv) so nothing is re-downloaded.
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

# 10 subclasses: deep + superficial cortical glutamatergic incl. L5 PT (=ET)
# and layer-6 types, plus the 4 main isocortical interneuron classes.
TARGET_SUBCLASSES = [
    'L2/3 IT CTX Glut',
    'L5 IT CTX Glut',
    'L5 ET CTX Glut',   # L5 ET = L5 PT (pyramidal tract)
    'L6 IT CTX Glut',
    'L6 CT CTX Glut',
    'L6b CTX Glut',
    'Pvalb Gaba',
    'Sst Gaba',
    'Vip Gaba',
    'Lamp5 Gaba',
]

# Full 10x receptor panel, grouped by neuromodulator + canonical G-protein
# coupling. Htr3a/b are ionotropic but kept under "5-HT iono" so the plot
# still shows the difference vs the metabotropic groups.
GENE_GROUPS = {
    '5-HT Gi':   ['Htr1a', 'Htr1b', 'Htr1d', 'Htr1f', 'Htr5a', 'Htr5b'],
    '5-HT Gq':   ['Htr2a', 'Htr2b', 'Htr2c'],
    '5-HT Gs':   ['Htr4', 'Htr6', 'Htr7'],
    '5-HT iono': ['Htr3a', 'Htr3b'],
    'NE Gq':     ['Adra1a', 'Adra1b', 'Adra1d'],
    'NE Gi':     ['Adra2a', 'Adra2b', 'Adra2c'],
    'NE Gs':     ['Adrb1', 'Adrb2', 'Adrb3'],
    'DA Gs':     ['Drd1', 'Drd5'],
    'DA Gi':     ['Drd2', 'Drd3', 'Drd4'],
}

# Looser threshold than the simplified plot so the "extra" 10x receptors
# (Adra1*, Htr1d, Htr2b, Htr3a/b, Adrb2, Drd3/4/5) still show up: drop
# only receptors that are essentially absent everywhere in the selected
# subclasses.
MIN_MAX_FRAC = 0.01
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
    print(f'\nFiltered to {REGION} + {len(TARGET_SUBCLASSES)} subclasses: '
          f'{len(meta_sel):,} cells')
    for sc_name in TARGET_SUBCLASSES:
        n = (meta_sel['subclass_short'] == sc_name).sum()
        print(f'  {sc_name}: {n:,}')

    available = set(expr_sel.columns)
    gene_groups = {
        grp: [g for g in genes if g in available]
        for grp, genes in GENE_GROUPS.items()
    }

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
    fig_w = (0.85 * len(flat_genes) + 6) * 1.4
    fig_h = (0.85 * n_sc + 3) * 1.4

    dp = sc.pl.dotplot(
        adata,
        var_names=gene_groups,
        groupby='subclass_short',
        cmap='Reds',
        figsize=(fig_w, fig_h),
        dot_min=0.25,
        dot_max=0.75,
        colorbar_title='Mean\nlog2(CPM+1)',
        size_title='% expressing\n\n',
        show=False,
        return_fig=True,
    )
    dp.style(dot_edge_color='black', dot_edge_lw=0.5,
             largest_dot=120 * (FONT_SCALE ** 1.5))
    dp.legends_width = 4.0 * FONT_SCALE


    out = OUTPUT_DIR / 'dotplot_mPFC_10x_receptors.png'
    dp.savefig(out, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f'\nSaved: {out}')


if __name__ == '__main__':
    main()
