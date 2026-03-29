"""Extract LSX expression from WMB-10Xv3-STR, save CSV, generate dotplot, commit+push."""
import pandas as pd, numpy as np, h5py, os, gc, json, subprocess, urllib.request, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import trim_mean
import anndata
from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import time as _time

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

lsx_cells = cell_extended[cell_extended['region_of_interest_acronym'] == 'LSX'].copy()
print(f"LSX cells: {len(lsx_cells):,}")

outdir = Path('scube1_lsx')
outdir.mkdir(exist_ok=True)
lsx_cells.to_csv(outdir / 'lsx_metadata.csv')
print("Saved lsx_metadata.csv")

gene_symbols = ['Scube1', 'Drd1', 'Drd2', 'Foxp2', 'Chat', 'Pvalb', 'Sst', 'Th', 'Aqp4', 'Mog']
gene_ids = gene[gene['gene_symbol'].isin(gene_symbols)].index.tolist()
gene_sym_map = dict(zip(gene.loc[gene_ids].index, gene.loc[gene_ids]['gene_symbol']))

manifest_path = download_base / abc_cache.current_manifest
with open(manifest_path) as f:
    manifest = json.load(f)

def git_commit_push(files, message):
    subprocess.run(['git', 'add'] + files, cwd='/home/user/abc_atlas_access', check=True)
    r = subprocess.run(['git', 'commit', '-m', message], cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
    if r.returncode != 0: print("  Nothing to commit"); return
    for attempt in range(4):
        r2 = subprocess.run(['git', 'push', '-u', 'origin', 'claude/regenerate-scube1-plots-8Bdh0'],
                           cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
        if r2.returncode == 0: print("  Pushed"); return
        _time.sleep(2 ** (attempt + 1))
    print(f"  Push failed: {r2.stderr}")

# Commit metadata
git_commit_push(['notebooks/scube1_lsx/lsx_metadata.csv'], 'Add LSX cell metadata CSV for Scube1 dotplot')

# Extract expression
h5ad_info = manifest['file_listing']['WMB-10Xv3']['expression_matrices']['WMB-10Xv3-STR']['log2']['files']['h5ad']
local_path = download_base / h5ad_info['relative_path']

if not local_path.exists():
    print(f"Downloading WMB-10Xv3-STR...")
    local_path.parent.mkdir(parents=True, exist_ok=True)
    for dl_attempt in range(5):
        try:
            urllib.request.urlretrieve(h5ad_info['url'], local_path)
            break
        except Exception as e:
            print(f"  Attempt {dl_attempt+1} failed: {e}")
            if local_path.exists(): os.remove(local_path)
            _time.sleep(2 ** (dl_attempt + 1))

print(f"Processing WMB-10Xv3-STR ({local_path.stat().st_size / 1e9:.1f} GB)...")
lsx_cell_labels = set(lsx_cells.index)

with h5py.File(local_path, 'r') as f:
    obs_names = f['obs']['cell_label'][:].astype(str)
    var_names = f['var']['gene_identifier'][:].astype(str)
    cell_mask = np.isin(obs_names, list(lsx_cell_labels))
    n_match = cell_mask.sum()
    print(f"  {n_match:,} LSX cells in file")

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
    df.to_csv(outdir / 'lsx_expression.csv')
    print(f"  Saved lsx_expression.csv ({df.shape})")

print("Deleting h5ad...")
os.remove(local_path)

git_commit_push(['notebooks/scube1_lsx/lsx_expression.csv'],
                f'Add LSX expression CSV ({n_match:,} cells x 10 genes)')

# === Generate dotplot ===
print("\nGenerating dotplot...")
expr = pd.read_csv(outdir / 'lsx_expression.csv', index_col=0)
meta = pd.read_csv(outdir / 'lsx_metadata.csv', index_col=0)
common = expr.index.intersection(meta.index)
expr, meta = expr.loc[common], meta.loc[common]

gene_list = list(expr.columns)
adata = anndata.AnnData(X=expr.values.astype(np.float32),
    obs=meta[['subclass','supertype','class','region_of_interest_acronym']].copy(),
    var=pd.DataFrame(index=gene_list))
adata.obs['subclass_short'] = adata.obs['subclass'].apply(lambda x: re.sub(r'^\d+\s+', '', str(x)))

# Exclude non-LSX noise
noise_patterns = [
    # Cortical
    'L2/3 IT', 'L4/5 IT', 'L5 IT', 'L5 ET', 'L5 NP', 'L6 ', 'L6b ',
    'IT EP-CLA', 'CLA-EPd', 'IT AON',
    'Pvalb Gaba', 'Pvalb chandelier', 'Sst Gaba', 'Vip Gaba',
    'Lamp5 Gaba', 'Lamp5 Lhx6', 'Sncg Gaba',
    # OB
    'OB ', 'OB-', 'Astro-OLF', 'OEC ',
    # Hippocampal
    'HPF ', 'DG ', 'DG-PIR',
    # Striatal (noise in LSX dissection)
    'STR D1 Gaba', 'STR D2 Gaba', 'STR D1 Sema5a', 'STR-PAL',
    'OT D3',
    # Hypothalamic
    'LHA ', 'LHA-AHN', 'AHN ', 'DMH ', 'VMH ', 'PVH', 'AVPV-MEPO',
    'HY Gnrh1', 'MPO-ADP', 'PH-',
    # Amygdala
    'CEA-', 'MEA-', 'MEA ', 'LA-BLA', 'COAa-PAA', 'RHP-COA',
    # Thalamic
    'TH Prkcd',
    # Midbrain
    'PAG', 'SC ', 'SNc-VTA', 'MRN',
    # Rare immune
    'DC NN', 'Lymphoid', 'Monocytes', 'ABC NN',
]

def is_noise(name): return any(pat in name for pat in noise_patterns)
noise_types = [s for s in adata.obs['subclass_short'].unique() if is_noise(s)]
print(f'Excluding {len(noise_types)} noise subclasses')
adata = adata[~adata.obs['subclass_short'].isin(noise_types)].copy()

# Also drop <10 cells
sub_counts = adata.obs.groupby('subclass_short', observed=True).size()
small = sub_counts[sub_counts < 10].index
if len(small) > 0:
    print(f'Excluding {len(small)} subclasses with <10 cells')
    adata = adata[~adata.obs['subclass_short'].isin(small)].copy()

adata.obs['subclass_short'] = pd.Categorical(adata.obs['subclass_short'])
print(f'Kept {adata.n_obs:,} cells, {adata.obs["subclass_short"].cat.categories.size} subclasses')

# Chat+ group
chat_mask = adata.X[:, list(adata.var_names).index('Chat')] > 0
ac = adata[chat_mask].copy()
ac.obs['subclass_short'] = 'Chat+ (manual select)'
adata_wc = anndata.concat([adata, ac])
cats = list(adata.obs['subclass_short'].cat.categories) + ['Chat+ (manual select)']
adata_wc.obs['subclass_short'] = pd.Categorical(adata_wc.obs['subclass_short'], categories=cats, ordered=True)

grp_counts = adata_wc.obs.groupby('subclass_short', observed=True).size().to_dict()
nc = len(cats)

# Dotplot
_CMAP = mcolors.LinearSegmentedColormap.from_list('wm', ['#FFFFFF','#FFCCEE','#FF66BB','#DD0088','#990066'])
gg = {'Primary': ['Scube1'], 'MSN markers': ['Drd1','Drd2','Foxp2'], 'Interneuron': ['Chat','Pvalb','Sst','Th'], 'Glia': ['Aqp4','Mog']}
fg, gs = [], []
for gn, gl in gg.items(): s = len(fg); fg.extend(gl); gs.append((s, len(fg)-1, gn))
ng = len(fg)
y_labels = [f'n={grp_counts.get(c,0):>6,} | {c}' for c in reversed(cats)]

edf = pd.DataFrame(adata_wc.X, index=adata_wc.obs.index, columns=adata_wc.var_names)
grp = adata_wc.obs['subclass_short'].values
tm, frac = np.zeros((nc,ng)), np.zeros((nc,ng))
for i, cat in enumerate(cats):
    m = grp == cat
    if m.sum() == 0: continue
    v = edf.loc[m, fg].values.astype(float)
    for j in range(ng): frac[i,j] = (v[:,j]>0).mean(); tm[i,j] = trim_mean(v[:,j], proportiontocut=0.25)

norm = mcolors.Normalize(vmin=0, vmax=11); ms, mi = 220, 12
fig, ax = plt.subplots(figsize=(16, max(6, nc*0.48)))
for i in range(nc):
    yi = nc-1-i
    for j in range(ng):
        f = frac[i,j]
        if f < 0.005: continue
        ax.scatter(j, yi, s=mi+(ms-mi)*min(f,1), c=[_CMAP(norm(min(tm[i,j],11)))], edgecolors='black', linewidths=0.5, zorder=3)
ax.set_yticks(range(nc)); ax.set_yticklabels(y_labels, fontsize=8, family='monospace')
ax.set_xticks(range(ng)); ax.set_xticklabels(fg, rotation=90, ha='center', fontsize=9)
ax.set_xlim(-0.5,ng-0.5); ax.set_ylim(-0.5,nc-0.5)
for s,e,gn in gs: ax.annotate(gn, xy=((s+e)/2,1.01), xycoords=('data','axes fraction'), ha='center', va='bottom', fontsize=10, fontweight='bold')
for j in range(ng+1): ax.axvline(j-0.5, color='gray', lw=0.3, alpha=0.3, zorder=1)
for i in range(nc+1): ax.axhline(i-0.5, color='gray', lw=0.3, alpha=0.3, zorder=1)
ax.set_title('Scube1 & markers — Lateral septum (LSX dissection)', fontsize=13, pad=30)
h = [Line2D([0],[0],marker='o',color='w',ls='None',ms=np.sqrt(mi+(ms-mi)*f)/1.7,mfc='black',mec='black',mew=0.5,label=f'{f}') for f in [0.2,0.4,0.6,0.8,1.0]]
leg = ax.legend(handles=h, title='Fraction of cell\nwith >0 read', loc='upper left', bbox_to_anchor=(1.02,1), frameon=True, fontsize=8, title_fontsize=9, handletextpad=0.5, labelspacing=1.2)
ax.add_artist(leg)
sm = plt.cm.ScalarMappable(cmap=_CMAP, norm=norm); sm.set_array([])
cax = inset_axes(ax, width='3%', height='30%', loc='lower left', bbox_to_anchor=(1.02,0,1,1), bbox_transform=ax.transAxes, borderpad=0)
cb = fig.colorbar(sm, cax=cax)
tl = list(range(0,12)); tlab = [str(t) for t in tl]; tlab[-1] = '>11'
cb.set_ticks(tl); cb.set_ticklabels(tlab); cb.set_label('Trimmed mean (25%-75%)\nLog2(CPM+1)', fontsize=9)
fig.savefig(outdir / 'dotplot_LSX_Scube1_by_subclass.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved dotplot')

git_commit_push(['notebooks/scube1_lsx/dotplot_LSX_Scube1_by_subclass.png',
                 'notebooks/extract_lsx_expression.py'],
                'Add Scube1 dotplot for lateral septum (LSX dissection)')

print("\nAll done!")
