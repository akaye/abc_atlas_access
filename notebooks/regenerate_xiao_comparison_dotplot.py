#!/usr/bin/env python3
"""
Render an mPFC Htr dot plot that mimics the layout of Xiao et al. (Kwan,
Nature 2025) Fig. 1a to make our data directly comparable.

Xiao panel: regions ACA + ALM + ORB + PL-ILA; subclasses L2/3 IT, L4/5 IT,
L5 IT, L6 IT, L5 PT, L6 CT; genes Htr1a, Htr2a, Htr2b, Htr2c (+ Slc17a7,
Gad1 as sanity markers, omitted here -- not in cache).

We use the cached extraction (PL-ILA-ORB + ACA only -- ALM is not pulled to
save time). Dot size is %-expressing clipped to [20, 100]; color is mean
log2(CPM+1) clipped at >=6 to match Xiao's scale.

Run: python notebooks/regenerate_xiao_comparison_dotplot.py
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

EXPR_CSV = NOTEBOOK_DIR / 'mpfc_receptor_expression.csv'
META_CSV = NOTEBOOK_DIR / 'mpfc_10x_metadata.csv'

REGIONS = ['PL-ILA-ORB', 'ACA']

# Order matches Xiao panel (top to bottom).
TARGET_SUBCLASSES = [
    'L2/3 IT CTX Glut',
    'L4/5 IT CTX Glut',
    'L5 IT CTX Glut',
    'L6 IT CTX Glut',
    'L5 ET CTX Glut',   # = L5 PT
    'L6 CT CTX Glut',
]
DISPLAY_LABELS = {
    'L2/3 IT CTX Glut': 'L2/3 IT',
    'L4/5 IT CTX Glut': 'L4/5 IT',
    'L5 IT CTX Glut':   'L5 IT',
    'L6 IT CTX Glut':   'L6 IT',
    'L5 ET CTX Glut':   'L5 PT',
    'L6 CT CTX Glut':   'L6 CT',
}

GENES = ['Htr1a', 'Htr2a', 'Htr2b', 'Htr2c']

# Xiao clipping
MIN_FRAC_DISPLAY = 0.20   # dots below 20% are not drawn (their legend starts at 20)
VMAX = 6.0                 # color clipped at >=6


def main():
    print(f'Loading {EXPR_CSV.name} ...')
    expr = pd.read_csv(EXPR_CSV, index_col=0)
    meta = pd.read_csv(META_CSV, index_col=0)
    meta = meta.loc[meta.index.intersection(expr.index)].copy()
    meta['subclass_short'] = meta['subclass'].apply(
        lambda x: re.sub(r'^\d+\s+', '', str(x))
    )

    mask = (
        meta['region_of_interest_acronym'].isin(REGIONS)
        & meta['subclass_short'].isin(TARGET_SUBCLASSES)
    )
    meta_sel = meta.loc[mask]
    expr_sel = expr.loc[meta_sel.index, GENES]
    print(f'Cells used (ACA + PL-ILA-ORB, 6 subclasses): {len(meta_sel):,}')
    for sc_name in TARGET_SUBCLASSES:
        n = (meta_sel['subclass_short'] == sc_name).sum()
        print(f'  {DISPLAY_LABELS[sc_name]:<8} ({sc_name}): {n:,}')

    # Compute mean log2(CPM+1) and fraction expressing per subclass.
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

    # ------- plot -------
    n_rows = len(TARGET_SUBCLASSES)
    n_cols = len(GENES)
    fig, ax = plt.subplots(figsize=(0.9 * n_cols + 2.5, 0.6 * n_rows + 2))

    # Xiao-like colormap (white -> pink -> deep magenta/black).
    cmap = LinearSegmentedColormap.from_list(
        'xiao_pink',
        ['#ffffff', '#f6c6df', '#e23a8c', '#a0166a', '#3a0a2b'],
    )
    norm = Normalize(vmin=0, vmax=VMAX)

    # Dot size: % expressing clipped to [MIN_FRAC_DISPLAY, 1].
    # Use a fixed pixel scaling so dots are visually proportional to %.
    SIZE_MIN_PCT = MIN_FRAC_DISPLAY * 100
    SIZE_MAX_PCT = 100.0
    AREA_MIN = 30
    AREA_MAX = 600

    def size_for(p):
        # p is fraction in [0,1]
        pct = p * 100
        if pct < SIZE_MIN_PCT:
            return None  # don't draw
        t = (pct - SIZE_MIN_PCT) / (SIZE_MAX_PCT - SIZE_MIN_PCT)
        t = float(np.clip(t, 0, 1))
        return AREA_MIN + t * (AREA_MAX - AREA_MIN)

    for i, sc_name in enumerate(TARGET_SUBCLASSES):
        y = n_rows - 1 - i  # top-to-bottom
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

    # Highlight the L5 PT row (yellow) like Xiao.
    pt_y = n_rows - 1 - TARGET_SUBCLASSES.index('L5 ET CTX Glut')
    ax.add_patch(plt.Rectangle(
        (-1.2, pt_y - 0.4), 0.9, 0.8,
        facecolor='#f7d65c', edgecolor='none', zorder=-1, clip_on=False,
    ))

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(GENES, style='italic')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [DISPLAY_LABELS[s] for s in TARGET_SUBCLASSES][::-1]
    )
    ax.set_xlim(-1.3, n_cols - 0.3)
    ax.set_ylim(-0.7, n_rows - 0.3)
    ax.set_title('ACA, ORB and PL-ILA (cached; no ALM)', fontsize=11)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)

    # Color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.12, shrink=0.7)
    cbar.set_label('Gene expression\nlog$_2$(CPM + 1)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Size legend (separate axes to the right)
    legend_ax = fig.add_axes([0.78, 0.55, 0.18, 0.35])
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis('off')
    legend_ax.set_title('Percentage\nof cells', fontsize=9, loc='center')
    for k, pct in enumerate([20, 40, 60, 80, 100]):
        y = 0.85 - k * 0.18
        legend_ax.scatter(0.3, y, s=size_for(pct / 100), c='black')
        legend_ax.text(0.55, y, f'{pct}', va='center', fontsize=8)

    out = OUTPUT_DIR / 'dotplot_mPFC_xiao_comparison.png'
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved: {out}')


if __name__ == '__main__':
    main()
