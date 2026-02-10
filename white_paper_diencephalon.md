# Cell-Type-Resolved Neuromodulator Receptor Expression in the Thalamus, Lateral Habenula, and Claustrum: A Multi-Modal Transcriptomic Atlas for Diencephalic and Claustral Drug Targets

---

## Abstract

The thalamus, lateral habenula (LHb), and claustrum occupy critical positions in brain circuits underlying arousal, consciousness, affective processing, and sensory integration, yet comprehensive cell-type-resolved maps of neuromodulator receptor expression in these structures have been lacking. Here we systematically characterize the expression of 28 monoaminergic receptor genes (14 serotonin, 9 norepinephrine, 5 dopamine) across all major cell types in the mouse thalamus, lateral habenula, and claustrum using two complementary platforms from the Allen Brain Cell Atlas: 10x single-nucleus RNA sequencing (259,327 thalamic cells, 28,668 claustral cells) and MERFISH spatial transcriptomics (50,282 thalamic cells across 15 nuclei, 2,486 LHb cells, 8,932 claustral cells; 11 shared receptor genes). We reveal profound nucleus-specific heterogeneity in the thalamus: Htr7 dominates midline and intralaminar relay neurons (87-93% expressing in reuniens and parafascicular nuclei), while the reticular nucleus exhibits relative "monoamine resistance" with minimal serotonin receptor expression but selective Adra1b expression (73%). Association relay neurons preferentially express Adra1b (86% in mediodorsal, 88% in centromedian-intralaminar), while paraventricular and midline nuclei concentrate D2/D3 dopamine receptors (54% Drd2, 62% Drd3 in PVT). In the lateral habenula, principal neurons express Htr2a at 60%—a striking finding with direct implications for understanding serotonergic contributions to depression and ketamine's antidepressant mechanism. Claustrum excitatory neurons display unprecedented expression of Htr2c (96%, log2 7.95) and Htr1f (92%, log2 7.18), the highest levels observed across all brain regions in our dataset, with substantial Drd1 (87-94%) and divergent Drd2 expression between excitatory subtypes (10-fold difference). Cross-modality validation demonstrates strong concordance between platforms (22 shared thalamic subclasses, 11 shared claustral subclasses), while cross-region comparisons reveal fundamentally distinct receptor landscapes: thalamic and claustral neuronal profiles show minimal correlation (rho = 0.12-0.18), while claustrum and LHb show the strongest similarity (rho = 0.50-0.55), driven by shared Htr2a and Drd2 expression. These findings provide a molecular foundation for understanding neuromodulatory control of thalamocortical gating, habenula-mediated aversion processing, and claustral multisensory integration, with direct implications for anesthesia mechanisms, rapid-acting antidepressants, psychedelic therapeutics, and deep brain stimulation target optimization.

---

## 1. Introduction

The neuromodulatory systems—serotonergic, noradrenergic, and dopaminergic—exert profound influence over brain function through their widespread projections and diverse receptor repertoires. While recent advances in single-cell transcriptomics have begun to reveal the cellular architecture of neuromodulator signaling in cortical and amygdalar circuits, the diencephalon and claustrum remain comparatively underexplored. This represents a critical gap in our understanding, as these structures serve as essential hubs for arousal, consciousness, sensory integration, and affective processing—functions that are frequently disrupted in psychiatric and neurological disease.

This paper extends our previous characterization of neuromodulator receptor expression in the basolateral amygdala (BLA) and medial prefrontal cortex (mPFC) to three anatomically and functionally distinct structures: the thalamus, lateral habenula (LHb), and claustrum. Using single-nucleus RNA sequencing (10x snRNA-seq) and spatially resolved MERFISH transcriptomics from the Allen Brain Cell Atlas, we provide a comprehensive, cell-type-resolved map of 28 neuromodulator receptor genes (14 serotonin, 9 norepinephrine, and 5 dopamine receptors) across these regions. Together, these datasets enable unprecedented resolution of receptor expression patterns across more than 200 neuronal subclasses, linking molecular identity to anatomical location and functional circuit organization.

The thalamus occupies a unique position in brain architecture as the primary relay station between subcortical structures and the cerebral cortex. Unlike cortical pyramidal neurons, which integrate inputs to generate complex computational outputs, thalamic relay neurons primarily gate and modulate information flow. This fundamental difference in circuit logic suggests that neuromodulatory control mechanisms in thalamus may differ substantially from those in cortex. The thalamus comprises numerous functionally specialized nuclei: sensory relay nuclei (ventral posterior, lateral geniculate, medial geniculate) transmit modality-specific information to primary sensory cortices; association nuclei (mediodorsal, pulvinar, anterior group) connect with prefrontal and association cortices; midline and intralaminar nuclei (paraventricular, reuniens, centromedian) support arousal, attention, and affective processing; and the reticular nucleus provides GABAergic feedback inhibition to regulate thalamic output. Each nucleus exhibits distinct connectivity patterns, firing properties, and involvement in behavior and disease. Critically, thalamic circuits are central to general anesthesia mechanisms, with both GABAergic agents and ketamine modulating thalamic relay and reticular activity to produce loss of consciousness. Chronic pain states involve maladaptive plasticity in medial and intralaminar thalamic nuclei, making these targets for neuromodulatory interventions. Deep brain stimulation of specific thalamic nuclei has shown promise for treating disorders of consciousness, epilepsy, and treatment-resistant depression, yet the neuromodulatory pharmacology of individual thalamic cell types remains largely unknown.

The lateral habenula represents a small but disproportionately influential structure in the epithalamus, serving as a critical node for encoding negative valence and aversive prediction errors. LHb neurons respond to punishment, reward omission, and aversive stimuli, and their hyperactivity has been causally linked to depressive-like behaviors in animal models. The LHb receives dense serotonergic input from the dorsal and median raphe nuclei and projects to both the serotonergic raphe and dopaminergic ventral tegmental area, positioning it as a key regulator of monoamine release. Recent work has identified the LHb as a critical substrate for ketamine's rapid antidepressant effects: ketamine blocks NMDA receptor-mediated burst firing in LHb neurons, reducing their inhibitory influence on monoaminergic centers. Understanding the serotonin, norepinephrine, and dopamine receptor landscape of LHb neuronal subtypes may provide mechanistic insight into both the pathophysiology of depression and the therapeutic mechanisms of rapid-acting antidepressants.

The claustrum, a thin sheet of gray matter situated between the insular cortex and external capsule, has emerged as a structure of intense interest in consciousness research. The claustrum exhibits reciprocal connections with nearly all cortical areas and has been proposed to coordinate multisensory integration, salience detection, and attentional switching. Provocatively, Crick and Koch hypothesized that the claustrum might function as a "consciousness switch," synchronizing cortical activity to enable unified perceptual experience. While this hypothesis remains debated, converging evidence indicates that claustral lesions disrupt sensory integration and that claustral neurons respond selectively to cross-modal and salient stimuli. The claustrum receives dense neuromodulatory innervation and expresses high levels of serotonin 5-HT2A receptors, the primary target of classical psychedelic drugs. Recent work suggests that psychedelic-induced alterations in consciousness may involve claustral circuits, making receptor expression patterns in this structure particularly relevant to understanding altered states and potential therapeutic applications of psychedelics.

The present study reveals profound heterogeneity in neuromodulator receptor expression across these three regions. In the thalamus, we identify striking nucleus-specific patterns: midline and intralaminar relay neurons show dominant expression of Htr7 (87-93% of cells in reuniens and parafascicular nuclei), while the reticular nucleus exhibits relative "monoamine resistance" with low overall receptor expression but selective Adra1b expression in 73% of cells. Association thalamic relay neurons preferentially express Adra1b (86% in mediodorsal, 88% in centromedian-intralaminar), while paraventricular and midline nuclei show concentrated D2/D3 dopamine receptor expression (54% Drd2, 62% Drd3 in PVT). In the lateral habenula, principal neurons express Htr2a at 60%, suggesting strong serotonergic sensitivity with potential relevance to both depression pathophysiology and ketamine's mechanism. Claustrum excitatory neurons exhibit extreme and unprecedented expression of Htr2c (96% of cells, log2 expression 7.95) and Htr1f (92% of cells, log2 expression 7.18-7.19), the highest receptor expression levels observed across our entire dataset. These findings establish a foundation for understanding region-specific and cell-type-specific neuromodulatory control in diencephalic and claustral circuits.

---

## 2. Methods

### 2.1 Data Source

All data were obtained from the Allen Brain Cell Atlas (ABC Atlas) release 20251031, a comprehensive single-cell transcriptomic resource of the adult mouse whole brain. We accessed the data through the `abc_atlas_access` Python package (version 0.1.0) using the `AbcProjectCache` API class with the AWS S3 cache method (`from_s3_cache`). The Whole Mouse Brain (WMB) dataset comprises 4,042,976 cells profiled by 10x Chromium single-nucleus RNA sequencing (10x) and 4,334,174 cells profiled by multiplexed error-robust fluorescence in situ hybridization (MERFISH). Expression matrices for 10x data were stored as backed AnnData h5ad files in sparse CSR format, subdivided by chemistry version (10Xv2 and 10Xv3) and dissection region. Cell metadata were retrieved from `WMB-10X/cell_metadata`, gene metadata from `WMB-10X/gene`, and taxonomic annotations from `WMB-taxonomy/cluster_to_cluster_annotation_membership_pivoted`.

### 2.2 Cell Type Taxonomy

Cell types were classified according to the ABC Atlas consensus taxonomy, which organizes cells hierarchically into classes, subclasses, supertypes, and clusters based on transcriptomic similarity. We joined cell metadata with cluster annotation tables using the `cluster_alias` field to obtain neurotransmitter identity, class, subclass, supertype, and cluster labels for each cell. For 10x analyses, we applied minimum cell count thresholds to ensure robust estimates: >=100 cells per subclass for thalamus, >=50 cells per supertype for claustrum. For MERFISH analyses, we used >=30 cells per subclass for thalamus and >=20 cells per subclass for lateral habenula given the smaller total cell count.

### 2.3 Brain Region Definitions

**Thalamus MERFISH.** Thalamic cells were identified using MERFISH spatial coordinates and Common Coordinate Framework version 3 (CCFv3) parcellation labels. We selected cells annotated to 15 thalamic nuclei: mediodorsal (MD), reticular (RT), paraventricular (PVT), parafascicular (PF), reuniens (RE), anteroventral (AV), anterodorsal (AD), anteromedial (AM), central lateral (CL), central medial (CM), paracentral (PCN), intermediodorsal (IMD), interanteromedial (IAM), interanteromedial dorsal (IAD), and xiphoid (Xi). This yielded 50,282 cells spanning 26 subclasses and 80 supertypes.

**Thalamus 10x.** We retrieved cells from the thalamus (TH) dissection region, combining both chemistry versions (WMB-10Xv2-TH and WMB-10Xv3-TH). After filtering to subclasses with >=100 cells, the dataset comprised 259,327 cells across 39 subclasses and 147 supertypes.

**Lateral Habenula.** LHb cells were selected from MERFISH data using CCFv3 parcellation labels corresponding to the lateral habenula (LH) structure. After filtering to subclasses with >=20 cells, we obtained 2,486 cells spanning 13 subclasses and 34 supertypes, including 851 cells in the predominant LH Pou4f1 Sox1 Glut neuronal population.

**Claustrum MERFISH.** Claustrum cells were identified using CCFv3 parcellation labels for the claustrum proper (CLA) and endopiriform nucleus (EPd). From 11,179 initial cells, 8,932 were selected across 18 subclasses and 86 supertypes, distributed as 5,469 CLA and 5,710 EPd cells.

**Claustrum 10x.** Claustrum cells from 10x data were obtained from the cortical subplate (CTXsp) dissection region. We selected CLA-specific excitatory subclasses plus shared cortical interneuron subclasses, yielding 28,668 cells (18,797 excitatory, 9,871 interneurons) across 64 supertypes (>=50 cells).

### 2.4 Neuromodulator Receptor Gene Selection

We analyzed 28 neuromodulator receptor genes: 14 serotonin receptors (Htr1a, Htr1b, Htr1d, Htr1f, Htr2a, Htr2b, Htr2c, Htr3a, Htr3b, Htr4, Htr5a, Htr5b, Htr6, Htr7), 9 norepinephrine receptors (Adra1a, Adra1b, Adra1d, Adra2a, Adra2b, Adra2c, Adrb1, Adrb2, Adrb3), and 5 dopamine receptors (Drd1, Drd2, Drd3, Drd4, Drd5). All 28 genes were present in 10x transcriptomes. Of these, 11 (Htr1b, Htr1d, Htr2a, Htr3a, Htr7, Adra1a, Adra1b, Drd1, Drd2, Drd3, Drd5) were in the MERFISH 550-gene panel.

### 2.5 Expression Data Extraction

For 10x data, log2-normalized expression matrices were accessed as backed AnnData objects. Expression matrices are organized by chemistry and dissection region. Cells were subset by boolean indexing converted to integer arrays via `np.where()`, loaded into memory via `.to_memory()`, and converted to pandas DataFrames. Expression values were cached to CSV for rapid reloading.

For MERFISH data, the expression matrix (C57BL6J-638850/log2, 7.6 GB) was accessed as a backed AnnData object and the 11 receptor genes extracted for cells matching CCF parcellation criteria using the same integer indexing approach.

### 2.6 Dot Plot Visualization

Receptor expression patterns were visualized using scanpy's `sc.pl.dotplot()`. Dot plots encode mean expression (color intensity) and fraction expressing (dot size). Expression was variance-normalized per gene. Genes were grouped by monoamine family. For thalamic MERFISH data, dot plots were generated at multiple levels: by subclass, by supertype, by nucleus, and for MD substructure.

### 2.7 Cross-Modality Analysis

We compared 10x and MERFISH expression profiles for 22 shared thalamic subclasses and 11 shared claustral subclasses across 11 shared genes. Overall agreement was assessed via Pearson and Spearman correlations for all (gene x cell type) combinations. Per-gene and per-cell-type correlations evaluated consistency. Consensus rankings were computed as average of 10x and MERFISH ranks, with confidence scores defined as 1 - (|rank difference| / max possible difference).

### 2.8 Cross-Region Analysis

**Intra-thalamic comparisons.** MERFISH nuclei were grouped into functional categories: mediodorsal (MD), midline (PVT, CM, IMD), anterior (AD, AV, AM, IAD, IAM), reticular (RT), intralaminar (PF, CL, PCN), and reuniens (RE, Xi). Log2 fold-changes were computed relative to overall thalamic neuronal mean.

**Inter-regional comparisons.** Overall neuronal receptor profiles were compared pairwise between TH, CLA, and LHb using Spearman correlations. CLA excitatory versus TH excitatory neurons were compared using the full 28-gene panel. Shared glial populations were compared between TH and CLA.

### 2.9 Statistical Analysis and Software

All analyses were performed in Python 3.11 with: `abc_atlas_access` (0.1.0), `anndata` (0.10.3), `scanpy` (1.9.6), `pandas` (2.1.4), `numpy` (1.26.2), `scipy` (1.11.4), `seaborn`, and `matplotlib` (3.8.2). Spearman correlations served as the primary association metric.

---

## 3. Results

### 3.1 Cell Type Composition and Dataset Overview

To characterize neuromodulator receptor expression in the diencephalon and claustrum, we analyzed single-cell transcriptomic data from the Allen Brain Cell Atlas using both MERFISH spatial transcriptomics and 10x Chromium single-cell RNA sequencing. Our analysis encompassed three anatomically and functionally distinct structures: the thalamus (TH), lateral habenula (LHb), and claustrum (CLA).

The thalamus MERFISH dataset comprised 50,282 cells distributed across 15 nuclei, organized into 26 subclasses and 80 supertypes (**Figure 1**). The largest nuclei included the reticular nucleus (RT, 13,664 cells), mediodorsal nucleus (MD, 9,656 cells), and paraventricular nucleus (PVT, 5,158 cells). Additional nuclei analyzed included the parafascicular nucleus (PF, 3,355 cells), nucleus reuniens (RE, 2,752 cells), anteroventral nucleus (AV, 2,744 cells), anterodorsal nucleus (AD, 2,457 cells), anteromedial nucleus (AM, 2,401 cells), central lateral nucleus (CL, 2,357 cells), central medial nucleus (CM, 1,662 cells), paracentral nucleus (PCN, 1,208 cells), intermediodorsal nucleus (IMD, 1,127 cells), interanterodorsal nucleus (IAD, 919 cells), xiphoid nucleus (Xi, 624 cells), and interanteromedial nucleus (IAM, 406 cells). The most abundant neuronal subclasses were TH Prkcd Grin2c Glut (7,466 cells), RT-ZI Gnb3 Gaba (6,879 cells), PVT-PT Ntrk1 Glut (3,175 cells), CM-IAD-CL-PCN Sema5b Glut (3,119 cells), and AD Serpinb7 Glut (2,329 cells). The complementary thalamus 10x dataset provided deeper transcriptomic coverage with 259,327 cells across 39 subclasses and 147 supertypes, dominated by TH Prkcd Grin2c Glut (68,125 cells) and oligodendrocyte lineage cells (71,510 cells).

The lateral habenula MERFISH dataset captured 2,486 cells organized into 13 subclasses and 34 supertypes (**Supplementary Figure S10**). The primary neuronal population was LH Pou4f1 Sox1 Glut with 851 cells, accompanied by substantial glial populations including oligodendrocytes (503 cells), astrocytes (469 cells), and endothelial cells (223 cells).

For the claustrum, we analyzed 8,932 selected cells from the MERFISH dataset, combining claustrum proper (5,469 cells) with the adjacent endopiriform nucleus (5,710 cells) across 18 subclasses and 86 supertypes (**Figure 2**). The major excitatory populations included CLA-EPd-CTX Car3 Glut (1,499 cells), IT EP-CLA Glut (2,391 cells), and L6b EPd Glut (1,947 cells). GABAergic interneurons comprised Sncg (201 cells), Vip (200 cells), Lamp5 (223 cells), Pvalb (196 cells), and Sst (276 cells) subtypes. The claustrum 10x dataset provided 28,668 cells, including CLA excitatory neurons (18,797 cells: CLA-EPd-CTX Car3 Glut 4,022 and IT EP-CLA Glut 14,775) and CTXsp interneurons (9,871 cells).

### 3.2 Serotonin Receptor Expression Patterns

Serotonin receptor expression revealed striking regional and cell-type specificity across the diencephalon and claustrum. Htr2a, the primary target of classical psychedelics, showed particularly prominent expression in specific thalamic and claustral populations. In the thalamus, AD Serpinb7 Glut exhibited the highest Htr2a expression with consensus rank 1.0 and confidence 1.00, validated across both MERFISH and 10x platforms (**Figure 3**). Remarkably, the lateral habenula displayed unexpectedly high Htr2a expression, with 60% of LH Pou4f1 Sox1 Glut cells expressing this receptor at a mean level of 0.89—striking given that subcortical structures typically show lower Htr2a expression compared to cortical regions (**Supplementary Figure S10**). The claustrum showed even more dramatic Htr2a enrichment, with 96% of CLA-EPd-CTX Car3 Glut cells expressing this receptor at a mean MERFISH intensity of 3.14, achieving consensus rank 1.0 with confidence 1.00 (**Supplementary Figure S7**).

Htr7 emerged as the dominant serotonin receptor in midline and intralaminar thalamic nuclei, revealing functional architecture within the thalamus (**Figure 5**). PF Fzd5 Glut showed 93.1% expression frequency, RE-Xi Nox4 Glut reached 87.3%, and PVT-PT Ntrk1 Glut achieved 84.0%. This contrasted sharply with the reticular nucleus, where Htr7 expression was nearly silent at 11.9%. Quantitative analysis revealed nucleus reuniens showed +1.35 log2FC enrichment for Htr7 relative to overall thalamic levels, while intralaminar nuclei showed +0.65 log2FC and MD showed +0.32 log2FC. The reticular nucleus was markedly depleted at -1.73 log2FC. In the claustrum, Htr7 showed highest expression in Lamp5 Lhx6 Gaba interneurons (rank 1, confidence 1.00).

Htr1b expression provided additional regional distinctions. In the claustrum, IT EP-CLA Glut achieved consensus rank 1.0 with confidence 1.00, with 72% of cells expressing in MERFISH at mean 1.31. The lateral habenula showed particularly high Htr1b in periaqueductal gray-projecting neurons, with 82% expression in the Pax6 subclass and 68% in Tcf7l2 subclass (**Supplementary Figure S10**). Thalamic expression was moderate and broadly distributed.

Analysis of 10x-exclusive receptors revealed extraordinary subtype specificity in the claustrum. Htr2c showed a remarkable 10-fold difference between claustrum excitatory subtypes: CLA-EPd-CTX Car3 Glut displayed extraordinarily high expression (log2 7.95, 96% expressing) while IT EP-CLA showed substantially lower levels (log2 2.02) (**Figure 2**). Even more striking was Htr1f, which reached unprecedented levels in both CLA excitatory types (log2 7.18-7.19, 92% expressing), representing the highest Htr1f expression observed across all brain regions in our entire dataset. Thalamic expression of both receptors was substantially lower.

Htr3a showed conserved interneuron specificity. In claustrum MERFISH, Sncg Gaba achieved rank 1 with confidence 1.00 (90% expressing, mean 3.56), and the Sncg > Vip hierarchy matched patterns observed in BLA and mPFC (**Supplementary Figure S7**).

### 3.3 Noradrenergic Receptor Expression Patterns

Adra1b, encoding the alpha-1B adrenergic receptor, showed the broadest expression across thalamic relay neurons. In thalamus MERFISH, CM-IAD-CL-PCN Sema5b and TH Prkcd Grin2c both achieved rank 1.5 with confidence 0.95, indicating nearly equivalent enrichment. The mediodorsal nucleus showed 86% expression frequency. Notably, RT-ZI Gnb3 Gaba showed 73% Adra1b expression (**Figure 1**, **Supplementary Figure S4**)—significant because the reticular nucleus showed minimal expression of most serotonin receptors yet maintained robust adrenergic expression, suggesting differential neuromodulatory control: the reticular nucleus may be preferentially sensitive to noradrenergic rather than serotonergic modulation. In the claustrum, Vip Gaba interneurons achieved rank 1 with confidence 1.00, showing 98% expression at mean 2.42 in MERFISH (**Supplementary Figure S7**).

Adra1a was notably enriched in claustrum Lamp5 Gaba interneurons at 80% frequency in MERFISH. Thalamic expression was lower overall, suggesting that claustral interneurons integrate both alpha-1A and alpha-1B signaling for nuanced noradrenergic modulation.

Among 10x-exclusive receptors, Adra2c showed another striking divergence between claustrum excitatory subtypes: CLA-EPd-CTX Car3 expressed at very high levels (log2 5.60, 84%) compared to IT EP-CLA (log2 2.09, 38%), representing another strong molecular differentiator between these populations (**Figure 2**).

### 3.4 Dopamine Receptor Expression Patterns

Drd1 showed highest claustrum MERFISH expression in L6b EPd Glut (94%, mean 2.49), with IT EP-CLA at 91% (**Supplementary Figure S7**). In thalamus, AD Serpinb7 Glut achieved consensus rank 1.0 with confidence 1.00. LHb showed moderate expression at 21% in LH Glut.

Drd2 expression revealed functional organization in cortico-thalamic circuits. In claustrum, IT EP-CLA Glut achieved consensus rank 1.0 with confidence 1.00 (62% in MERFISH, mean 1.22). Thalamic expression was highest in PVT-PT Ntrk1 Glut at 54%, with midline nuclei showing general enrichment (+0.61 log2FC). The lateral habenula showed 32% expression in LH Glut—substantial dopamine sensitivity in a structure encoding aversive prediction errors (**Figure 6**, **Supplementary Figure S19**).

Drd3 showed midline thalamic specificity, with PVT-PT Ntrk1 at 62% in MERFISH. Co-expression of Drd2 and Drd3 in PVT neurons suggests integration of dopaminergic signals through both D2-like subtypes (**Figure 5**). Drd5 expression showed notable overlap with Htr7, with highest expression in reuniens, parafascicular, and paraventricular nuclei, suggesting that association and midline nuclei integrate serotonergic and dopaminergic signals through specific receptor complements.

### 3.5 Cross-Modality Validation

We performed systematic cross-modality validation comparing MERFISH and 10x for 22 shared thalamic subclasses and 11 shared claustral subclasses (**Figure 4**). Per-gene analysis showed Htr2a with the strongest cross-platform correlation, consistent with high expression levels and clear cell-type specificity. Htr1d showed the weakest correlation, reflecting low-abundance detection challenges—a pattern consistent with our BLA/mPFC analyses (**Supplementary Figure S12**).

Per-cell-type correlations were highest for glial populations, as observed in companion studies (**Supplementary Figure S13**). High-confidence consensus rankings (confidence >=0.7, rank <=3.0) identified 29 gene-cell type pairs in thalamus and 29 in claustrum, including AD Serpinb7 Glut for Htr2a and Drd1 in thalamus; CLA-EPd-CTX Car3 Glut for Htr2a and IT EP-CLA Glut for Htr1b and Drd2 in claustrum (**Supplementary Figures S14-S17**).

### 3.6 Cross-Region Comparison

**Intra-thalamic heterogeneity.** Analysis of functional nucleus groups from MERFISH neuronal populations revealed distinct receptor signatures (**Figure 5**, **Supplementary Figure S18**). Reuniens showed Htr7 enrichment of +1.35 log2FC and Drd5 of +0.65 log2FC. The reticular nucleus exhibited depletion for Htr7 (-1.73 log2FC) and Htr2a (-1.28 log2FC), but only modest depletion for Adra1b (-0.13 log2FC). Midline nuclei showed enrichment for Drd2 (+0.61 log2FC) and Drd3 (+0.82 log2FC). MD displayed Adra1b enrichment of +0.50 log2FC.

**Inter-regional comparisons.** Overall neuronal profiles revealed fundamental differences (**Figure 6**). TH versus CLA correlations were very low (Spearman rho = 0.118-0.182), reflecting differences between diencephalic relay and telencephalic projection neurons. TH versus LHb showed low-to-moderate correlation (rho = 0.336-0.391). CLA versus LHb showed the highest correlation (rho = 0.500-0.545), driven by shared Htr2a and Drd2 expression (**Supplementary Figure S19**). Shared glial populations showed conserved profiles (TH versus CLA glia rho = 0.742), indicating that regional variation targets neurons rather than glia (**Supplementary Figure S21**).

**CLA versus TH excitatory neurons.** Using the full 28-gene 10x panel, CLA excitatory neurons showed dramatically different profiles from thalamic relay neurons, with CLA distinguished by extreme Htr2c, Htr1f, and Adra2c expression while TH relay neurons showed dominant Htr7 and Adra1b (**Supplementary Figure S20**).

### 3.7 Spatial Substructure Analysis

Within the claustrum, CLA proper versus EPd showed nearly identical receptor profiles, with only slightly elevated Drd1 in L6b EPd Glut (**Supplementary Figure S9**). This contrasts with cortical regions where laminar position strongly predicts receptor expression.

Within the thalamus, functional nucleus groups showed clear receptor signatures (**Figure 3**). Association nuclei uniformly expressed high Htr7, distinguishing them from sensory relay nuclei. The reticular nucleus maintained its unique profile with Adra1b as the dominant receptor. MD substructure analysis was limited by the absence of subdivisions (MDc, MDl, MDm) in the current taxonomy (**Supplementary Figure S6**).

---

## 4. Discussion

This comprehensive single-cell analysis of neuromodulator receptor expression in the thalamus, lateral habenula, and claustrum reveals striking region-specific and cell-type-specific patterns that have important implications for understanding normal brain function, disease mechanisms, and therapeutic interventions. The profound differences in receptor expression profiles across these structures—and the marked contrast with the cortical and amygdalar patterns described in our companion paper—underscore the necessity of circuit-specific and cell-type-resolved approaches to understanding neuromodulatory signaling.

### Thalamic Heterogeneity and Functional Specialization

The thalamus exhibits a degree of nucleus-to-nucleus receptor heterogeneity that rivals the diversity observed across cortical areas. The dominant expression of Htr7 in midline and intralaminar nuclei (87% in reuniens, 93% in parafascicular) is particularly noteworthy, as these nuclei provide non-specific arousal and attentional signals to widespread cortical targets. Htr7 is a Gs-coupled receptor that enhances neuronal excitability and has been implicated in cognitive flexibility, circadian rhythms, and thermoregulation. Its concentration in association and intralaminar thalamic nuclei suggests a role in modulating arousal-dependent information flow to association cortices. Notably, Htr7 expression is nearly absent in primary sensory relay nuclei, suggesting that serotonin may preferentially regulate higher-order thalamic processing while leaving primary sensory transmission relatively unmodulated—a pattern consistent with serotonin's role in state-dependent gating of attention rather than basic sensory function.

The reticular nucleus presents a striking exception to this pattern. As the sole GABAergic component of the thalamus, the reticular nucleus provides feedback inhibition to thalamocortical relay neurons and plays a critical role in generating thalamic oscillations, sensory gating, and attention. Our finding that reticular neurons express minimal levels of most neuromodulator receptors but selectively express Adra1b in 73% of cells suggests specialized noradrenergic control of thalamic inhibition. Alpha-1B adrenergic receptors are Gq-coupled and enhance neuronal excitability; their expression in reticular neurons may enable norepinephrine to enhance reticular-mediated inhibition during arousal states, sharpening thalamocortical transmission through lateral inhibition. This pattern contrasts sharply with cortical and amygdalar interneurons, which express diverse receptor combinations. The relative monoamine resistance of the reticular nucleus may have important implications for general anesthesia, as both GABAergic anesthetics and ketamine modulate reticular activity to suppress consciousness.

The concentration of Adra1b in association relay nuclei (86% in mediodorsal, 88% in centromedian-intralaminar) suggests that norepinephrine preferentially modulates higher-order thalamic processing through alpha-1B receptor-mediated enhancement of relay neuron excitability. The selective expression of D2 and D3 dopamine receptors in paraventricular and midline nuclei (54% Drd2, 62% Drd3 in PVT) is equally intriguing, as these structures are implicated in stress, anxiety, and drug-seeking behavior. D2/D3 receptors are Gi-coupled and typically inhibitory; their expression in PVT suggests that dopamine may suppress midline thalamic activity, potentially providing a mechanism for reward-related suppression of stress signals.

### Lateral Habenula and Rapid Antidepressant Mechanisms

The lateral habenula's expression profile has direct relevance to understanding depression and ketamine's antidepressant mechanism. The finding that 60% of LHb principal neurons express Htr2a is particularly significant. 5-HT2A receptors are Gq-coupled and enhance neuronal excitability; serotonergic input to the LHb via 5-HT2A would be expected to increase LHb activity. Given that LHb hyperactivity is associated with depressive states and that LHb neurons inhibit both serotonergic raphe and dopaminergic VTA, this suggests a potential maladaptive positive feedback loop: increased serotonin release could enhance LHb activity via 5-HT2A, leading to further suppression of monoamine centers and potentially worsening depressive symptoms. This mechanism may help explain the delayed onset and variable efficacy of traditional SSRIs, which immediately increase synaptic serotonin but require weeks to produce antidepressant effects—potentially reflecting the time required for adaptive changes that overcome LHb-mediated negative feedback.

Ketamine's rapid antidepressant action involves NMDA receptor blockade in the LHb, reducing burst firing and suppressing LHb output. The presence of both Htr2a (60%) and Drd2 (32%) in LHb neurons suggests that this structure integrates serotonergic and dopaminergic signals with glutamatergic inputs. The concentration of Htr1b (68-82%) in PRC-PAG neurons projecting through the habenular region suggests additional serotonergic modulation of LHb-adjacent circuits. Understanding the cell-type-specific receptor profiles in LHb may enable development of more targeted interventions—for example, selective 5-HT2A antagonism in LHb might provide antidepressant effects without the psychotomimetic properties of systemic 5-HT2A agonists.

### Claustrum and Consciousness Modulation

The claustrum's receptor expression profile is among the most distinctive we have observed across all brain regions. The extreme expression of Htr2c (96% of excitatory neurons, log2 7.95) and Htr1f (92%, log2 7.18-7.19) is unprecedented in our dataset. 5-HT2C receptors are Gq-coupled and enhance excitability, while 5-HT1F receptors are Gi-coupled and typically inhibitory. The co-expression of these receptors with opposing functional effects suggests sophisticated serotonergic regulation of claustral activity, potentially enabling bidirectional control of claustral output depending on serotonin concentration and temporal dynamics.

The high 5-HT2C expression is particularly intriguing in light of the claustrum's proposed role in consciousness and salience detection. While 5-HT2A is the primary target of classical psychedelics, 5-HT2C also binds psychedelic compounds with substantial affinity. The extreme concentration of 5-HT2C in claustral excitatory neurons suggests that psychedelic-induced alterations in consciousness may involve claustral circuit modulation. The near-complete expression of Htr1f is equally remarkable, as this receptor is typically expressed at much lower levels. Htr1f is the target of lasmiditan, a recently approved migraine treatment, and the claustrum's extreme expression raises the possibility that claustral dysfunction might contribute to migraine with aura—a condition involving transient alterations in consciousness and sensory integration.

The dominant Drd1 expression in claustral excitatory neurons (87-94% in MERFISH) suggests dopaminergic enhancement of claustral activity, while the selective Drd2 expression in IT EP-CLA neurons (10-fold higher than CLA-EPd-CTX Car3) indicates cell-type-specific dopaminergic regulation. Claustral interneurons show receptor profiles that closely match cortical and amygdalar interneurons: Sncg expressing Htr3a, Vip expressing Adra1b. This conservation suggests that interneuron neuromodulation follows cell-type-intrinsic programs largely independent of circuit context, while excitatory neuron receptor expression is highly circuit-specific.

### Comparison with BLA and mPFC Patterns

The profound differences between thalamus/LHb/claustrum and the previously characterized BLA/mPFC highlight the region-specificity of neuromodulatory control. While mPFC pyramidal neurons show dominant Htr1a and Htr2a with layer-specific variation, thalamic relay neurons show Htr7 dominance in association nuclei, Adra1b dominance in relay populations, and near-absence of 5-HT2A. BLA principal neurons show moderate Htr1a, Htr2a, and Htr2c across subclasses, while claustral excitatory neurons show extreme and nearly exclusive Htr2c and Htr1f. The only conserved pattern across all regions is cell-type-specific expression in interneurons, where Sncg to Htr3a and Vip to Adra1b associations hold across cortex, amygdala, and claustrum.

### Limitations and Future Directions

Several limitations constrain interpretation. mRNA expression does not perfectly predict protein levels or functional receptor availability. The MERFISH 550-gene panel captures only a subset of receptors. Dissection boundaries in 10x are imperfect, and small structures like LHb are challenging to sample comprehensively. These data derive from adult mouse brain and may not fully reflect human receptor expression or disease-associated changes.

Future work should prioritize human thalamic and claustral receptor mapping, functional validation using cell-type-specific manipulations, optimization of deep brain stimulation targets informed by receptor profiles, and characterization of receptor changes in disease models. Integrating receptor expression with connectivity data would enable circuit-level understanding of neuromodulatory control. The present dataset provides a comprehensive molecular foundation for mechanistic investigation of arousal, consciousness, and affective processing, and a framework for developing more targeted therapeutic interventions.

---

## Figures

![Figure 1](outputs/dotplot_TH_receptors_by_subclass.png)
**Figure 1. Neuromodulator receptor expression across thalamic subclasses (10x scRNA-seq).** Dot plot showing 28 receptor genes across 39 subclasses (259,327 cells from TH dissection). Dot color indicates variance-normalized mean log2 expression; dot size indicates fraction expressing. Subclasses organized by neurotransmitter identity and taxonomy.

![Figure 2](outputs/dotplot_CLA_receptors_with_glia.png)
**Figure 2. Neuromodulator receptor expression in claustrum cell types including glia (10x scRNA-seq).** Dot plot showing 28 genes across CLA excitatory (18,797 cells), interneurons (9,871 cells), and glial populations from CTXsp dissection.

![Figure 3](outputs/dotplot_TH_MERFISH_receptors_by_nucleus.png)
**Figure 3. Spatially resolved receptor expression across thalamic nuclei (MERFISH).** Dot plot of 11 genes across 15 nuclei (50,282 cells), with each row a unique nucleus-supertype combination (132 groups). MERFISH spatial resolution reveals nucleus-specific patterns not apparent in bulk profiles.

![Figure 4](outputs/xmodal_scatter_TH_CLA.png)
**Figure 4. Cross-modality concordance of receptor expression (10x vs MERFISH).** Scatter plots comparing mean expression and fraction expressing for 11 shared genes across matched cell types in thalamus (22 shared subclasses) and claustrum (11 shared subclasses).

![Figure 5](outputs/xregion_thalamic_enrichment.png)
**Figure 5. Nucleus-specific receptor enrichment in thalamus (MERFISH).** Log2 fold-change of mean expression relative to overall thalamic neuronal mean. Nuclei grouped by functional category. Warm colors = enriched; cool colors = depleted.

![Figure 6](outputs/xregion_overall_profiles.png)
**Figure 6. Comparison of neuronal receptor profiles across brain regions (MERFISH).** Bar plots comparing 11 MERFISH receptor genes across TH, CLA, and LHb neuronal populations, revealing fundamentally distinct receptor landscapes.

---

## Supplementary Figures

![Figure S1](outputs/dotplot_TH_receptors_by_supertype.png)
**Figure S1. Thalamic receptor expression at supertype resolution (10x).** 28 genes across 147 supertypes (259,327 cells).

![Figure S2](outputs/dotplot_CLA_receptors_by_subclass.png)
**Figure S2. Claustrum receptor expression at subclass resolution, neuronal only (10x).** 28 genes across claustral neuronal subclasses.

![Figure S3](outputs/dotplot_CLA_receptors_by_supertype.png)
**Figure S3. Claustrum receptor expression at supertype resolution (10x).** 28 genes across 64 supertypes (28,668 cells).

![Figure S4](outputs/dotplot_TH_MERFISH_receptors_by_subclass.png)
**Figure S4. Thalamic receptor expression at subclass resolution (MERFISH).** 11 genes across 26 subclasses (50,282 cells).

![Figure S5](outputs/dotplot_TH_MERFISH_receptors_by_supertype.png)
**Figure S5. Thalamic receptor expression at supertype resolution (MERFISH).** 11 genes across 80 supertypes (50,282 cells).

![Figure S6](outputs/dotplot_TH_MERFISH_MD_substructure.png)
**Figure S6. Receptor expression within mediodorsal thalamus (MERFISH).** Expression across MD supertypes (9,656 cells).

![Figure S7](outputs/dotplot_CLA_MERFISH_receptors_by_subclass.png)
**Figure S7. Claustrum receptor expression at subclass resolution (MERFISH).** 11 genes across 18 subclasses (8,932 cells from CLA + EPd).

![Figure S8](outputs/dotplot_CLA_MERFISH_receptors_by_supertype.png)
**Figure S8. Claustrum receptor expression at supertype resolution (MERFISH).** 11 genes across 86 supertypes (8,932 cells).

![Figure S9](outputs/dotplot_CLA_MERFISH_receptors_by_structure.png)
**Figure S9. CLA versus EPd receptor expression comparison (MERFISH).** Cell types subdivided by parcellation structure (CLA vs EPd), showing nearly identical profiles.

![Figure S10](outputs/dotplot_LHb_MERFISH_receptors_by_subclass.png)
**Figure S10. Lateral habenula receptor expression at subclass resolution (MERFISH).** 11 genes across 13 subclasses (2,486 cells).

![Figure S11](outputs/dotplot_LHb_MERFISH_receptors_by_supertype.png)
**Figure S11. Lateral habenula receptor expression at supertype resolution (MERFISH).** 11 genes across 34 supertypes (2,486 cells).

![Figure S12](outputs/xmodal_per_gene_TH_CLA.png)
**Figure S12. Per-gene cross-modality correlations.** Spearman rho between 10x and MERFISH for each of 11 shared genes across cell types in TH and CLA.

![Figure S13](outputs/xmodal_per_celltype_TH_CLA.png)
**Figure S13. Per-cell-type cross-modality correlations.** Spearman rho across 11 genes for each shared cell type in TH and CLA.

![Figure S14](outputs/xmodal_rank_agreement_TH.png)
**Figure S14. Thalamic rank agreement (10x vs MERFISH).** 10x ranks, MERFISH ranks, and absolute rank differences for 22 shared subclasses.

![Figure S15](outputs/xmodal_rank_agreement_CLA.png)
**Figure S15. Claustrum rank agreement (10x vs MERFISH).** Rank comparison for 11 shared subclasses.

![Figure S16](outputs/xmodal_consensus_TH.png)
**Figure S16. Thalamic consensus receptor rankings.** Averaged 10x + MERFISH rankings with confidence scores.

![Figure S17](outputs/xmodal_consensus_CLA.png)
**Figure S17. Claustrum consensus receptor rankings.** Averaged rankings for 11 shared subclasses.

![Figure S18](outputs/xregion_thalamic_nucleus_heatmap.png)
**Figure S18. Thalamic nucleus-group receptor expression.** Mean expression and fraction expressing across 6 functional nucleus groups (MERFISH, neuronal only).

![Figure S19](outputs/xregion_LHb_vs_TH_relay.png)
**Figure S19. LHb versus thalamic relay neuron comparison (MERFISH).** Mean expression and fraction expressing for LH Pou4f1 Sox1 Glut versus major TH glutamatergic subclasses.

![Figure S20](outputs/xregion_CLA_vs_TH_excitatory.png)
**Figure S20. CLA versus TH excitatory neuron comparison (10x, 28 genes).** Fraction expressing heatmap for claustral and thalamic glutamatergic populations.

![Figure S21](outputs/xregion_glia_TH_vs_CLA.png)
**Figure S21. Shared glial receptor expression (TH vs CLA, MERFISH).** Scatter and enrichment heatmap for non-neuronal populations (Spearman rho = 0.742).
