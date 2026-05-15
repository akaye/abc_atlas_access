#!/usr/bin/env python3
"""
BLA companion to the Xiao-style mPFC dot plot.

Uses cached BLA 10x extraction (region_of_interest_acronym == CTXsp): major
BLA glutamatergic subclasses + the same interneurons (Pvalb, Sst, Vip) and
the same 5-HT receptors (Htr1a, Htr1b, Htr2a, Htr2b, Htr2c, Htr4, Htr7).

Run: python notebooks/regenerate_bla_xiao_style_dotplot.py
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

NOTEBOOK_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = NOTEBOOK_DIR.parent / 'outputs'
OUTPUT_DIR.mkdir(exist_ok=True)

EXPR_CSV = NOTEBOOK_DIR / 'bla_10x_neuronal_expression.csv'
META_CSV = NOTEBOOK_DIR / 'bla_10x_metadata.csv'

TARGET_SUBCLASSES = [
    'LA-BLA-BMA-PA Glut',
    'MEA-COA-BMA Ccdc42 Glut',
    'Pvalb Gaba',
    'Sst Gaba',
    'Vip Gaba',
]
DISPLAY_LABELS = {
    'LA-BLA-BMA-PA Glut':       'LA/BLA/BMA/PA Glut',
    'MEA-COA-BMA Ccdc42 Glut':  'MEA/COA/BMA Glut',
    'Pvalb Gaba':               'Pvalb',
    'Sst Gaba':                 'Sst',
    'Vip Gaba':                 'Vip',
}

GENES = ['Htr1a', 'Htr1b', 'Htr2a', 'Htr2b', 'Htr2c', 'Htr4', 'Htr7']

MIN_FRAC_DISPLAY = 0.20
VMAX = 6.0


def main():
    print(f'Loading {EXPR_CSV.name} ...')
    expr = pd.read_csv(EXPR_CSV, index_col=0)
    meta = pd.read_csv(META_CSV, index_col=0)
    meta = meta.loc[meta.index.intersection(expr.index)].copy()
    meta['subclass_short'] = meta['subclass'].apply(
        lambda x: re.sub(r'^\d+\s+', '', str(x))
    )

    mask = meta['subclass_short'].isin(TARGET_SUBCLASSES)
    meta_sel = meta.loc[mask]
    expr_sel = expr.loc[meta_sel.index, GENES]
    print(f'Cells used (CTXsp): {len(meta_sel):,}')
    for sc_name in TARGET_SUBCLASSES:
        n = (meta_sel['subclass_short'] == sc_name).sum()
        print(f'  {DISPLAY_LABELS[sc_name]:<22} ({sc_name}): {n:,}')

    df = expr_sel.copy()
    df['_sc'] = meta_sel['subclass_short'].values
    mean_expr = df.groupby('_sc', observed=True)[GENES].mean()
    frac_expr = df.groupby('_sc', observed=True)[GENES].apply(
        lambda x: (x > 0).mean()
    )
    mean_expr = mean_expr.loc[TARGET_SUBCLASSES]
    frac_expr = frac_expr.loc[TARGET_SUBCLASSES]

    print('\nMean log2(CPM+1):')
    print(mean_expr.round(2).to_string())
    print('\nFraction expressing:')
    print(frac_expr.round(2).to_string())

    n_rows = len(TARGET_SUBCLASSES)
    n_cols = len(GENES)

    fig = plt.figure(figsize=(0.95 * n_cols + 5.5, 0.55 * n_rows + 2.5))
    ax = fig.add_axes([0.18, 0.15, 0.50, 0.72])

    cmap = LinearSegmentedColormap.from_list(
        'greens',
        ['#f7fcf5', '#c7e9c0', '#74c476', '#238b45', '#00441b'],
    )
    norm = Normalize(vmin=0, vmax=VMAX)

    SIZE_MIN_PCT = MIN_FRAC_DISPLAY * 100
    SIZE_MAX_PCT = 100.0
    AREA_MIN = 30
    AREA_MAX = 600

    def size_for(p):
        pct = p * 100
        if pct < SIZE_MIN_PCT:
            return None
        t = (pct - SIZE_MIN_PCT) / (SIZE_MAX_PCT - SIZE_MIN_PCT)
        t = float(np.clip(t, 0, 1))
        return AREA_MIN + t * (AREA_MAX - AREA_MIN)

    for i, sc_name in enumerate(TARGET_SUBCLASSES):
        y = n_rows - 1 - i
        for j, gene in enumerate(GENES):
            m = mean_expr.loc[sc_name, gene]
            f = frac_expr.loc[sc_name, gene]
            s = size_for(f)
            if s is None:
                continue
            ax.scatter(
                j, y, s=s,
                c=[cmap(norm(min(m, VMAX)))],
                edgecolors='none',
            )

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(GENES, style='italic')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [DISPLAY_LABELS[s] for s in TARGET_SUBCLASSES][::-1]
    )
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-0.7, n_rows - 0.3)
    ax.set_title('BLA', fontsize=11)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)

    size_ax = fig.add_axes([0.73, 0.50, 0.10, 0.40])
    size_ax.set_xlim(0, 1)
    size_ax.set_ylim(0, 1)
    size_ax.axis('off')
    size_ax.text(0.5, 1.02, 'Percentage\nof cells',
                 ha='center', va='bottom', fontsize=9)
    for k, pct in enumerate([20, 40, 60, 80, 100]):
        y = 0.85 - k * 0.18
        size_ax.scatter(0.25, y, s=size_for(pct / 100), c='black')
        size_ax.text(0.55, y, f'{pct}', va='center', fontsize=8)

    cbar_ax = fig.add_axes([0.90, 0.50, 0.025, 0.40])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Gene expression\nlog$_2$(CPM + 1)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    out = OUTPUT_DIR / 'dotplot_BLA_xiao_style.png'
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved: {out}')


if __name__ == '__main__':
    main()
