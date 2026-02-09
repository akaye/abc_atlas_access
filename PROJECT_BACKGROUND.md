# Project Background: Neuromodulator Receptor Expression Atlas (BLA & mPFC)

This document describes the development history, technical approach, and complete file inventory for this project. It is intended to provide full context for continuing work in a new repository or Claude Code session.

---

## 1. Project Goal

The overarching goal is to create a cell-type-resolved atlas of neuromodulator receptor gene expression in two brain regions central to anxiety and fear circuitry: the **basolateral amygdala (BLA)** and the **medial prefrontal cortex (mPFC)**. The work uses the **Allen Brain Cell Atlas** (ABC Atlas) whole mouse brain dataset, analyzing both **10x single-nucleus RNA sequencing** and **MERFISH spatial transcriptomics** data. The ultimate aim is to understand how psychedelic drugs (targeting Htr2a, Htr1a, Htr2c, Htr7) and other psychiatric medications (targeting noradrenergic and dopaminergic receptors) differentially engage specific cell types within these circuits.

## 2. Data Source

- **Allen Brain Cell Atlas**: Release 20251031, accessed via `AbcProjectCache` API from the `abc_atlas_access` Python package
- **10x snRNA-seq**: 4,042,976 cells, 32,285 genes, organized by chemistry (10Xv2/v3) and dissection region
- **MERFISH**: 4,334,174 cells, 550-gene targeted panel, spatially registered to Allen CCF
- **S3 bucket**: `allen-brain-cell-atlas` (public)
- **28 receptor genes analyzed**: 14 serotonin (Htr1a/b/d/f, Htr2a/b/c, Htr3a/b, Htr4, Htr5a/b, Htr6, Htr7), 9 adrenergic (Adra1a/b/d, Adra2a/b/c, Adrb1/2/3), 5 dopamine (Drd1-5)
- **11 genes in MERFISH panel**: Htr1b, Htr1d, Htr2a, Htr3a, Htr7, Adra1a, Adra1b, Drd1, Drd2, Drd3, Drd5
- **Critically missing from MERFISH**: Htr1a, Htr2c (key psychedelic-relevant receptors)

## 3. Development History (Chronological)

### Phase 1: BLA 10x Dot Plot
1. Created initial notebook (`single_cell_dotplot_BLA_receptors.ipynb`) for BLA receptor expression
2. First attempt used all CTXsp cells (~122k across ~90 subclasses) — too broad
3. Refined to BLA-specific excitatory subclasses (LA-BLA-BMA-PA Glut, MEA-COA-BMA Ccdc42 Glut) plus shared GABAergic interneurons
4. Learned that 10x interneurons come from entire CTXsp dissection (cannot localize to BLA specifically)
5. Generated subclass-level and supertype-level dot plots
6. Added glia (astrocytes, microglia) for complete cellular profiling
7. **Key technical lesson**: Backed AnnData cannot chain views. Must use `adata[cell_idx, gene_idx]` with integer arrays from `np.where()`

### Phase 2: mPFC 10x Dot Plot
1. Created `single_cell_dotplot_mPFC_receptors.ipynb`
2. mPFC spans two dissection regions: PL-ILA-ORB (106k cells) + ACA (103k cells)
3. Expression data spread across 6 Isocortex h5ad files (v2: 1-4, v3: 1-2), each 8-12 GB
4. **Major challenge**: 30 GB disk cannot hold all files simultaneously
5. **Solution**: Process one h5ad at a time using helper script, extract subset, delete h5ad, repeat
6. Used h5py direct sparse matrix loading (much faster than backed anndata for this case)
7. Pre-extracted expression data to CSV for fast re-runs
8. 10 excitatory + 8 interneuron subclasses = 165,539 cells, 101 supertypes

### Phase 3: MERFISH Dot Plots (BLA + mPFC)
1. Created `single_cell_dotplot_BLA_MERFISH_receptors.ipynb`
   - Used CCF parcellation (`parcellation_structure == 'BLA'`) for spatial localization
   - 8,238 cells confirmed within BLA (vs. ~28k from dissection-based 10x)
   - Added substructure breakdown: BLAa, BLAp, BLAv — unique to MERFISH
   - Interneurons are now spatially confirmed within BLA (unlike 10x CTXsp)
2. Created `single_cell_dotplot_mPFC_MERFISH_receptors.ipynb`
   - Parcellation to PL, ILA, ACAd, ACAv areas
   - 64,474 cells with spatial localization
   - Area-specific breakdown showing PL vs ILA vs ACA expression differences

### Phase 4: CSV Caching
1. **Problem**: Each notebook run required re-downloading multi-GB h5ad files from S3
2. Multiple attempts at CSV caching timed out during execution
3. **Solution that worked**: Edit all 4 notebooks first (without executing), then execute one at a time in order from easiest to hardest
4. Execution order: (1) BLA MERFISH, (2) mPFC MERFISH, (3) BLA 10x, (4) mPFC 10x
5. Each notebook now checks for CSV existence before downloading h5ad files
6. Pattern: `if os.path.exists(csv_path): load CSV; else: extract from h5ad and save CSV`
7. Also saves metadata CSVs with subclass/supertype/class/neurotransmitter columns
8. Had to remove `*.csv` from `.gitignore` to track CSVs in git

### Phase 5: Cross-Modality Correlation (10x vs MERFISH)
1. Created `cross_modality_receptor_correlation.ipynb`
2. Compares 10x and MERFISH for same region using 11 shared genes
3. Analyses: scatter plots, per-gene correlations, per-cell-type correlations, rank comparisons, rank agreement heatmaps, consensus ranking dot plots
4. **Key finding**: Strong cross-platform agreement (BLA rho=0.85-0.87, mPFC rho=0.87-0.88)
5. Generates `consensus_receptor_rankings.csv` with confidence-weighted rankings
6. 55 high-confidence gene x cell type combinations identified

### Phase 6: Cross-Region Correlation (BLA vs mPFC)
1. Created `cross_region_receptor_correlation.ipynb`
2. Compares BLA vs mPFC for shared interneuron + glia types within each modality
3. MERFISH: 8 shared types, 11 genes; 10x: 10 shared types, 28 genes
4. Excitatory neurons excluded (region-specific taxonomy: BLA glut vs cortical layer types)
5. **Key findings**:
   - Interneuron profiles highly conserved (10x rho=0.94)
   - Htr2c most BLA-enriched receptor (log2 FC +1.6)
   - Adra1b most mPFC-enriched receptor (MERFISH log2 FC -1.2)
   - Adra2a strongly BLA-enriched (+1.1)
6. Generates enrichment heatmaps and `cross_region_enrichment.csv`

### Phase 7: White Paper
1. Used 3 parallel writing agents:
   - Agent 1: Introduction + Discussion (psychedelic emphasis, noradrenergic system)
   - Agent 2: Results (7 subsections with exact numbers from notebooks)
   - Agent 3: Methods (9 subsections) + Figure Legends (6 main + 11 supplementary)
2. First attempt used background agents whose transcripts were cleaned up before retrieval
3. Second attempt ran all 3 agents synchronously — all succeeded
4. Assembled into `white_paper.md` with Abstract
5. Generated `white_paper.pdf` (6.8 MB) using weasyprint with all 17 figures embedded

## 4. Technical Patterns & Lessons Learned

### Allen Brain Cell Atlas API
- `AbcProjectCache.from_s3_cache()` for S3 access
- Cell metadata: `abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata')`
- Gene metadata: `abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene')`
- Taxonomy: `abc_cache.get_metadata_dataframe(directory='WMB-taxonomy', file_name='cluster_to_cluster_annotation_membership_pivoted')`
- Join: `cell.join(cluster_details, on='cluster_alias')` gives neurotransmitter/class/subclass/supertype/cluster

### Expression Matrix Access
- 10x h5ad files organized by chemistry + region (e.g., `WMB-10Xv2-CTXsp`, `WMB-10Xv3-CTXsp`)
- MERFISH expression: `C57BL6J-638850/log2` (7.6 GB single file)
- **Backed AnnData gotcha**: Cannot chain views `adata[cells, :][:, genes]`. Must use single index `adata[cell_idx, gene_idx]` with integer arrays from `np.where()`
- For very large files (Isocortex): use h5py directly to load sparse CSR matrix, then subset
- h5ad sparse CSR stored as: `X/data`, `X/indices`, `X/indptr`; obs in `obs/cell_label`; var in `var/gene_identifier` (NOT `_index`)

### Region Definitions
- **10x BLA**: Filter by `region_of_interest_acronym == 'CTXsp'`, then select BLA-specific excitatory subclasses + shared interneurons
- **10x mPFC**: Filter by `region_of_interest_acronym` in `['PL-ILA-ORB', 'ACA']`
- **MERFISH BLA**: `parcellation_structure == 'BLA'`; substructures: BLAa, BLAp, BLAv
- **MERFISH mPFC**: `parcellation_substructure` in `['PL', 'ILA', 'ACAd', 'ACAv']`

### Mouse Gene Nomenclature
- First-letter capitalized: `Htr1a`, `Drd1`, `Adra2a` (not HTR1A or htr1a)

### Scanpy Dot Plots
- `sc.pl.dotplot(adata, var_names=gene_dict, groupby='subclass', standard_scale='var')`
- `var_names` accepts dict of gene groups for organized x-axis
- `standard_scale='var'` normalizes per gene for cross-gene comparison

### Disk Space Management
- 30 GB total disk; Isocortex h5ad files are 8-12 GB each (6 files)
- Cannot hold all simultaneously; process one at a time and delete
- CSV caching eliminates need to re-download after first extraction

## 5. Complete File Inventory

### Notebooks (in `notebooks/`)
| File | Description | Cells | Genes |
|------|-------------|-------|-------|
| `single_cell_dotplot_BLA_receptors.ipynb` | BLA 10x dot plots (subclass, supertype, with glia) | 54,113 | 28 |
| `single_cell_dotplot_mPFC_receptors.ipynb` | mPFC 10x dot plots (subclass, supertype, with glia) | 187,279 | 28 |
| `single_cell_dotplot_BLA_MERFISH_receptors.ipynb` | BLA MERFISH (subclass, supertype, substructure) | 8,238 | 11 |
| `single_cell_dotplot_mPFC_MERFISH_receptors.ipynb` | mPFC MERFISH (subclass, supertype, area) | 64,474 | 11 |
| `cross_modality_receptor_correlation.ipynb` | 10x vs MERFISH comparison | Both | 11 shared |
| `cross_region_receptor_correlation.ipynb` | BLA vs mPFC comparison | Shared types | 11/28 |

### CSV Data Files (in `notebooks/`)
| File | Rows | Description |
|------|------|-------------|
| `bla_10x_neuronal_expression.csv` | 28,810 | BLA neuronal expression (28 genes) |
| `bla_10x_glia_expression.csv` | 25,305 | BLA glia expression (28 genes) |
| `bla_10x_metadata.csv` | 54,114 | BLA cell metadata (subclass, supertype, etc.) |
| `mpfc_receptor_expression.csv` | 165,540 | mPFC neuronal expression (28 genes) |
| `mpfc_glia_receptor_expression.csv` | 21,741 | mPFC glia expression (28 genes) |
| `mpfc_10x_metadata.csv` | 187,280 | mPFC cell metadata |
| `bla_merfish_expression.csv` | 8,239 | BLA MERFISH expression (11 genes) |
| `bla_merfish_metadata.csv` | 8,239 | BLA MERFISH cell metadata |
| `mpfc_merfish_expression.csv` | 64,475 | mPFC MERFISH expression (11 genes) |
| `mpfc_merfish_metadata.csv` | 64,475 | mPFC MERFISH cell metadata |
| `consensus_receptor_rankings.csv` | 309 | Cross-modality consensus rankings |
| `cross_region_enrichment.csv` | 369 | Cross-region fold-change data |

### Output Figures (in `outputs/`)

**Main Dot Plots (9 files)**:
- `dotplot_BLA_receptors_by_subclass.png` / `_by_supertype.png` / `_with_glia.png`
- `dotplot_mPFC_receptors_by_subclass.png` / `_by_supertype.png` / `_with_glia.png`
- `dotplot_BLA_MERFISH_receptors_by_subclass.png` / `_by_supertype.png` / `_by_substructure.png`
- `dotplot_mPFC_MERFISH_receptors_by_subclass.png` / `_by_supertype.png` / `_by_area.png`

**Cross-Modality Figures (9 files)**:
- `correlation_scatter_10x_vs_merfish.png` — Overall scatter
- `correlation_per_gene.png` — Per-gene Spearman rho bars
- `correlation_per_celltype.png` — Per-cell-type correlation bars
- `rank_comparison_BLA.png` / `_mPFC.png` — Rank bump charts
- `rank_agreement_BLA.png` / `_mPFC.png` — Rank agreement heatmaps
- `consensus_rankings_BLA.png` / `_mPFC.png` — Consensus ranking dot plots

**Cross-Region Figures (8 files)**:
- `cross_region_scatter_merfish.png` / `_10x.png` — BLA vs mPFC scatters
- `cross_region_per_gene_merfish.png` / `_10x.png` — Per-gene regional correlations
- `cross_region_ranks_merfish.png` / `_10x.png` — Regional rank agreement
- `cross_region_enrichment_merfish.png` / `_10x.png` — Log2 fold-change heatmaps

### White Paper
- `white_paper.md` — Complete markdown source
- `white_paper.pdf` — PDF with all figures embedded (6.8 MB)
- `generate_pdf.py` — Script to regenerate PDF from markdown

### Other
- `RESULTS_SECTION.txt` — Raw results text (intermediate)
- `.gitignore` — Ignores `*.h5ad`, `notebooks/*.png`, `data/`; allows `*.csv`

## 6. Key Scientific Findings (Summary)

### Psychedelic-Relevant Receptors
- **Htr2a** (primary psychedelic target): Highest in Pvalb + Lamp5 interneurons (BLA), deep-layer IT neurons (mPFC). Cross-platform rho = 0.95-0.96.
- **Htr1a** (anxiolytic target, 10x only): Concentrated in Sst Chodl and Pvalb chandelier cells. May buffer against excessive Htr2a activation.
- **Htr2c** (opposes Htr2a, 10x only): Dramatically enriched in Sncg interneurons and BLA excitatory neurons. Most region-differentiated receptor (BLA >> mPFC, log2 FC +1.6).
- **Htr7** (fear extinction): Enriched in deep-layer corticothalamic neurons and chandelier cells. BLA-enriched (+0.7 in both modalities).

### Noradrenergic System
- **Adrb1** (propranolol target): Dominant in mPFC excitatory neurons. mPFC-enriched.
- **Adra1a/1b** (prazosin targets): Preferentially in interneurons (Lamp5, Pvalb, Vip). Adra1b strongly mPFC-enriched.
- **Adra2a/2c** (clonidine/guanfacine targets): BLA-enriched, especially in excitatory neurons.

### Cross-Validation
- 10x vs MERFISH agreement: rho = 0.79-0.88 overall; rho > 0.9 for key receptors
- BLA vs mPFC interneuron profiles: highly conserved (10x rho = 0.94)
- Cell type identity is primary determinant of expression; anatomy is secondary

## 7. Potential Next Steps

These were discussed but not yet implemented:
1. Extend to other anxiety-relevant regions (ventral hippocampus, BNST, PAG)
2. Integrate with Allen Connectivity Atlas for projection-specific analysis
3. Compare mouse vs human receptor profiles using ASAP/SeaAD datasets
4. Functional validation predictions for cell-type-specific manipulations
5. Co-expression analysis (which receptors are co-expressed in same cells)
6. Weighted receptor interaction networks per cell type
7. Comparison with bulk tissue pharmacological data

## 8. Environment & Dependencies

- Python 3.11
- `abc_atlas_access` 0.1.0 (Allen Institute package)
- `anndata` 0.10.3, `scanpy` 1.9.6
- `pandas` 2.1.4, `numpy` 1.26.2, `scipy` 1.11.4
- `matplotlib` 3.8.2, `weasyprint` 68.1
- Data from AWS S3: `allen-brain-cell-atlas` bucket (public, no credentials needed)
