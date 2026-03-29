"""Extract sAMY expression, generate dotplot, commit+push."""
import pandas as pd, numpy as np, h5py, os, gc, json, subprocess, urllib.request, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import trim_mean
import anndata, time as _time
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

ROI = 'sAMY'
roi_cells = cell_extended[cell_extended['region_of_interest_acronym'] == ROI].copy()
print(f"{ROI} cells: {len(roi_cells):,}")

outdir = Path('scube1_samy')
outdir.mkdir(exist_ok=True)
roi_cells.to_csv(outdir / 'samy_metadata.csv')

gene_symbols = ['Scube1', 'Drd1', 'Drd2', 'Foxp2', 'Chat', 'Pvalb', 'Sst', 'Th', 'Aqp4', 'Mog']
gene_ids = gene[gene['gene_symbol'].isin(gene_symbols)].index.tolist()
gene_sym_map = dict(zip(gene.loc[gene_ids].index, gene.loc[gene_ids]['gene_symbol']))

manifest_path = download_base / abc_cache.current_manifest
with open(manifest_path) as f:
    manifest = json.load(f)

def git_cp(files, msg):
    subprocess.run(['git', 'add'] + files, cwd='/home/user/abc_atlas_access', check=True)
    r = subprocess.run(['git', 'commit', '-m', msg], cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
    if r.returncode != 0: print("  Nothing to commit"); return
    for a in range(4):
        r2 = subprocess.run(['git', 'push', '-u', 'origin', 'claude/regenerate-scube1-plots-8Bdh0'],
                           cwd='/home/user/abc_atlas_access', capture_output=True, text=True)
        if r2.returncode == 0: print("  Pushed"); return
        _time.sleep(2 ** (a + 1))

git_cp(['notebooks/scube1_samy/samy_metadata.csv'], f'Add {ROI} metadata CSV')

# sAMY is in WMB-10Xv3-STR
h5info = manifest['file_listing']['WMB-10Xv3']['expression_matrices']['WMB-10Xv3-STR']['log2']['files']['h5ad']
local_path = download_base / h5info['relative_path']
if not local_path.exists():
    print("Downloading WMB-10Xv3-STR...")
    local_path.parent.mkdir(parents=True, exist_ok=True)
    for a in range(5):
        try: urllib.request.urlretrieve(h5info['url'], local_path); break
        except Exception as e:
            print(f"  Attempt {a+1} failed: {e}")
            if local_path.exists(): os.remove(local_path)
            _time.sleep(2 ** (a+1))

print(f"Processing ({local_path.stat().st_size/1e9:.1f} GB)...")
roi_labels = set(roi_cells.index)
with h5py.File(local_path, 'r') as f:
    obs = f['obs']['cell_label'][:].astype(str)
    var = f['var']['gene_identifier'][:].astype(str)
    mask = np.isin(obs, list(roi_labels))
    ci = np.where(mask)[0]; labels = obs[ci]
    print(f"  {len(ci):,} cells")
    gi, gn = [], []
    for gid in gene_ids:
        idx = np.where(var == gid)[0]
        if len(idx): gi.append(idx[0]); gn.append(gene_sym_map[gid])
    ed = np.zeros((len(ci), len(gi)), dtype=np.float32)
    X = f['X']; data, indices, indptr = X['data'], X['indices'], X['indptr']
    for ri, c in enumerate(ci):
        s, e = indptr[c], indptr[c+1]
        r_idx, r_dat = indices[s:e], data[s:e]
        for j, g in enumerate(gi):
            p = np.searchsorted(r_idx, g)
            if p < len(r_idx) and r_idx[p] == g: ed[ri, j] = r_dat[p]
    df = pd.DataFrame(ed, index=labels, columns=gn)
    df.to_csv(outdir / 'samy_expression.csv')
    print(f"  Saved ({df.shape})")

os.remove(local_path); gc.collect()
git_cp(['notebooks/scube1_samy/samy_expression.csv'], f'Add {ROI} expression CSV ({len(ci):,} cells)')

# === Dotplot ===
print("Generating dotplot...")
expr = pd.read_csv(outdir / 'samy_expression.csv', index_col=0)
meta = pd.read_csv(outdir / 'samy_metadata.csv', index_col=0)
common = expr.index.intersection(meta.index)
expr, meta = expr.loc[common], meta.loc[common]
adata = anndata.AnnData(X=expr.values.astype(np.float32),
    obs=meta[['subclass','supertype','class','region_of_interest_acronym']].copy(),
    var=pd.DataFrame(index=list(expr.columns)))
adata.obs['ss'] = adata.obs['subclass'].apply(lambda x: re.sub(r'^\d+\s+', '', str(x)))

# Noise: cortical, OB, hippocampal, thalamic, midbrain, hindbrain
noise_p = ['L2/3','L4/5','L5 ','L6 ','L6b','IT EP-CLA','CLA-EPd','IT AON',
    'Pvalb Gaba','Pvalb chan','Sst Gaba','Vip Gaba','Lamp5 Gaba','Lamp5 Lhx6','Sncg Gaba',
    'OB ','OB-','Astro-OLF','OEC','HPF ','DG ','DG-PIR',
    'TH ','RT-ZI','LGv','PVT','MH ','LH ','Pineal',
    'PAG','SC ','SNc-VTA','MRN','IC ','MB-',
    'PB ','NTS','MY ','P ','CB',
    'DC NN','Lymphoid','Monocytes','ABC NN',
    'HY Gnrh1','VMH ','PVR ','SCiw','SCig','SCsg','SCs ',
    'OB-STR-CTX','RHP-COA','LSX ','MS-SF',
]
noise = [s for s in adata.obs['ss'].unique() if any(p in s for p in noise_p)]
adata = adata[~adata.obs['ss'].isin(noise)].copy()
sc = adata.obs.groupby('ss', observed=True).size()
adata = adata[~adata.obs['ss'].isin(sc[sc < 10].index)].copy()
adata.obs['ss'] = pd.Categorical(adata.obs['ss'])
print(f"Kept {adata.n_obs:,} cells, {adata.obs['ss'].cat.categories.size} subclasses")

chat_mask = adata.X[:, list(adata.var_names).index('Chat')] > 0
ac = adata[chat_mask].copy(); ac.obs['ss'] = 'Chat+ (manual select)'
aw = anndata.concat([adata, ac])
cats = list(adata.obs['ss'].cat.categories) + ['Chat+ (manual select)']
aw.obs['ss'] = pd.Categorical(aw.obs['ss'], categories=cats, ordered=True)
gc_d = aw.obs.groupby('ss', observed=True).size().to_dict()

_CMAP = mcolors.LinearSegmentedColormap.from_list('wm', ['#FFFFFF','#FFCCEE','#FF66BB','#DD0088','#990066'])
gg = {'Primary':['Scube1'],'MSN markers':['Drd1','Drd2','Foxp2'],'Interneuron':['Chat','Pvalb','Sst','Th'],'Glia':['Aqp4','Mog']}
fg, gs = [], []
for gname, gl in gg.items(): s=len(fg); fg.extend(gl); gs.append((s,len(fg)-1,gname))
nc, ng = len(cats), len(fg)
yl = [f'n={gc_d.get(c,0):>6,} | {c}' for c in reversed(cats)]
edf = pd.DataFrame(aw.X, index=aw.obs.index, columns=aw.var_names)
grp = aw.obs['ss'].values
tm, frac = np.zeros((nc,ng)), np.zeros((nc,ng))
for i, cat in enumerate(cats):
    m = grp == cat
    if m.sum() == 0: continue
    v = edf.loc[m, fg].values.astype(float)
    for j in range(ng): frac[i,j]=(v[:,j]>0).mean(); tm[i,j]=trim_mean(v[:,j], proportiontocut=0.25)
norm = mcolors.Normalize(vmin=0, vmax=11); ms, mi = 220, 12
fig, ax = plt.subplots(figsize=(16, max(6, nc*0.48)))
for i in range(nc):
    yi = nc-1-i
    for j in range(ng):
        f = frac[i,j]
        if f < 0.005: continue
        ax.scatter(j, yi, s=mi+(ms-mi)*min(f,1), c=[_CMAP(norm(min(tm[i,j],11)))], edgecolors='black', linewidths=0.5, zorder=3)
ax.set_yticks(range(nc)); ax.set_yticklabels(yl, fontsize=8, family='monospace')
ax.set_xticks(range(ng)); ax.set_xticklabels(fg, rotation=90, ha='center', fontsize=9)
ax.set_xlim(-0.5,ng-0.5); ax.set_ylim(-0.5,nc-0.5)
for s,e,gname in gs: ax.annotate(gname, xy=((s+e)/2,1.01), xycoords=('data','axes fraction'), ha='center', va='bottom', fontsize=10, fontweight='bold')
for j in range(ng+1): ax.axvline(j-0.5, color='gray', lw=0.3, alpha=0.3, zorder=1)
for i in range(nc+1): ax.axhline(i-0.5, color='gray', lw=0.3, alpha=0.3, zorder=1)
ax.set_title('Scube1 & markers — Subcortical amygdala (sAMY dissection)\n(central, medial & cortical amygdala)', fontsize=13, pad=30)
h = [Line2D([0],[0],marker='o',color='w',ls='None',ms=np.sqrt(mi+(ms-mi)*f)/1.7,mfc='black',mec='black',mew=0.5,label=f'{f}') for f in [0.2,0.4,0.6,0.8,1.0]]
leg = ax.legend(handles=h, title='Fraction of cell\nwith >0 read', loc='upper left', bbox_to_anchor=(1.02,1), frameon=True, fontsize=8, title_fontsize=9, handletextpad=0.5, labelspacing=1.2)
ax.add_artist(leg)
sm = plt.cm.ScalarMappable(cmap=_CMAP, norm=norm); sm.set_array([])
cax = inset_axes(ax, width='3%', height='30%', loc='lower left', bbox_to_anchor=(1.02,0,1,1), bbox_transform=ax.transAxes, borderpad=0)
cb = fig.colorbar(sm, cax=cax)
tl=list(range(0,12)); tlab=[str(t) for t in tl]; tlab[-1]='>11'
cb.set_ticks(tl); cb.set_ticklabels(tlab); cb.set_label('Trimmed mean (25%-75%)\nLog2(CPM+1)', fontsize=9)
fig.savefig(outdir / 'dotplot_sAMY_Scube1_by_subclass.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved dotplot')
git_cp(['notebooks/scube1_samy/dotplot_sAMY_Scube1_by_subclass.png', 'notebooks/extract_samy_expression.py'],
       'Add Scube1 dotplot for subcortical amygdala (sAMY)')
print("Done!")
