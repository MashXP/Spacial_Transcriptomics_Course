# Spatial Transcriptomics Sequencing Service Extraction

**Source:** `Spatial Transcriptomics Sequencing Service_20250917 (1).pdf`
**Date:** 2026-05-19

---

## Batch 1: Pages 1-5

### Page 1: Title Slide
- **Company:** Tri-I Biotech, Inc. 源資國際生物科技股份有限公司
- **Title:** Spatial Transcriptomics Sequencing Service
- **Presenter:** 源資國際生物科技股份有限公司 陳良坤
- **Contact:** TEL: +886-2-27053433, FAX: +886-2-27053431, URL: http://www.tri-ibiotech.com.tw

### Page 2: Company Overview & Services
- **Company Tagline:** Total Solution Provider of Systems Biology and Chemoinformatics.
- **Services Offered by Tri-I Biotech:**
    1. Sanger Sequencing
    2. NGS Sequencing
    3. Oligo synthesis
    4. Genotyping & SNP analysis
    5. Gene Expression analysis (qPCR & RNA-seq)
    6. Bioinformatics analysis
- **Partner Products & Tools:**
    - **Systems Biology:** GeneGO (USA), Softgenetics (USA), Funakoshi (Japan).
    - **Chemoinformatics:** Certara (USA), ACD Labs (Canada), Gaussian (USA).
    - **Bioinformatics:** Genious (USA), Bioinformatics solution (Canada), Gene Codes (USA), Prism (USA).
- **Accreditation:** TAF (Taiwan Accreditation Foundation).

### Page 3: 10x Genomics Visium System
- **Topic:** Spatial transcriptomics – Visium system
- **Partner:** 10x Genomics

### Page 4: Workflow Overview
- **Title:** Transform histological sections into whole transcriptome data
- **Subtitle:** End-to-end solution for single cell-scale spatial transcriptomics
- **Workflow Steps:**
    1. **Sample prep & imaging:** Tissue section on glass slide.
    2. **CytAssist capture & spatial barcoding:** Transfer of RNA to Visium slide.
    3. **Reverse transcription or probe extension:** Synthesis of cDNA.
    4. **Library construction:** Preparation for sequencing.
    5. **Sequencing:** High-throughput sequencing.
    6. **Data analysis & visualization:** Mapping data to tissue image.

### Page 5: Academic Reference & Applications
- **Reference:** Cao et al., "Spatial Transcriptomics: A Powerful Tool in Disease Understanding and Drug Discovery", *Theranostics* 2024, Vol. 14, Issue 7. DOI: 10.7150/thno.95908.
- **Figure 1 Summary (Application in Drug R&D):**
    - **Diseases Study:** Understanding disease mechanisms.
    - **Precision Medicine:** Personalized medicine and therapeutic response.
    - **Drug Development:** Target discovery, drug discovery, pharmacological models, pathway prediction, and drug action networks.

## Batch 2: Pages 6-10

### Page 6: Comparison of Visium Types
- **Title:** Three types of Visium
- **Comparison Table Summary:**
    - **Resolution:**
        - **Visium WT Panel:** 55 µm per spot (1-10 cell level).
        - **Visium HD WT Panel:** 2 µm per bin (single cell level).
        - **Visium HD 3':** 2 µm per bin (single cell level).
    - **Applicable Species:**
        - **Visium WT Panel:** Human & Mouse (Whole Transcriptome Probe).
        - **Visium HD WT Panel:** Human & Mouse (Whole Transcriptome Probe).
        - **Visium HD 3':** All Species agnostic 3' whole transcriptome (Poly-A mRNA capture).
    - **Applicable Samples:**
        - **Visium WT Panel:** FFPE, fresh-frozen, fixed-frozen tissues.
        - **Visium HD WT Panel:** FFPE, fresh-frozen, fixed-frozen tissues.
        - **Visium HD 3':** Only fresh frozen tissues.
    - **Acceptable Specimen Types:**
        - **Visium WT Panel:** FFPE wax block, FFPE H&E stained tissue slide, FFPE IF stained tissue slide.
        - **Visium HD WT Panel:** FFPE wax block, FFPE H&E stained tissue slide, FFPE IF stained tissue slide.
        - **Visium HD 3':** OCT H&E stained tissue slide.

### Page 7: Workflow & QC Requirements
- **Title:** Workflow
- **Process Flow:** Sample Record -> Sample QC -> 10X Visium Experiment -> NGS Sequencing -> Data Analysis.
- **QC Requirements:**
    - If QC fails, resend sample.
    - **RNA Extraction** required.
    - **RNA Quality Requirements:**
        - Visium WT Panel: DV200 > 30%
        - Visium HD WT Panel: DV200 > 30%
        - Visium HD 3' Panel: RIN >= 7
- **NGS Platform:** Illumina NovaSeq X Plus (PE150).

### Page 8: Demo Report Transition
- **Title:** Demo report for Single cell RNA-seq service (Transition slide).

### Page 9: Space Ranger Summary Report - Metrics
- **Title:** Space Ranger analysis_ Summary report
- **Sample:** Visium_HD_Human_Colon_Cancer (Human Whole Transcriptome Probe Set).
- **Key Metrics:**
    - Number of 8 µm binned squares under tissue: **516,880**
    - Mean reads per 8 µm bin: **1,395.3**
    - Mean UMIs per 8 µm bin: **523.0**
    - Total genes detected: **18,118**
- **Mapping Metrics:**
    - Reads Mapped to Probe Set: **98.8%**
    - Reads Mapped Confidently to Probe Set: **98.2%**
    - Reads Mapped Confidently to the Filtered Probe Set: **96.3%**
    - Reads Half-Mapped to Probe Set: **0.5%**
    - Reads Split-Mapped to Probe Set: **0.0%**
- **Sequencing Metrics:**
    - Number of Reads: **721,219,716**
    - Valid Barcodes: **90.0%**
    - Valid UMI Sequences: **99.9%**
    - Sequencing Saturation: **56.5%**
    - Q30 Bases in Barcode: **94.8%**
    - Q30 Bases in Probe Read: **94.1%**
    - Q30 Bases in UMI: **95.2%**
    - Fraction of Bins Under Tissue 8 µm: **73.6%**
    - UMIs per sq mm of Tissue: **8,172,085**
    - Reads per sq mm of Tissue: **19,457,827**
    - Fraction Reads in Squares Under Tissue: **98.9%**
- **Cell Segmentation Metrics:**
    - Number of Cells: **220,703**
    - Reads in Cells: **80.0%**
    - UMIs in Cells: **87.9%**
    - Mean Reads per Cell: **2,586.2**
    - Median Genes per Cell: **580.0**
    - Median UMIs per Cell: **759.0**
    - Median Cell Area (µm²): **96.0**
    - Median Nucleus Area (µm²): **28.0**
    - Maximum Nucleus Diameter (pixels): **256**

### Page 10: Space Ranger Summary Report - Cell Segmentation
- **Title:** Space Ranger analysis_ Summary report
- **Content:** Screenshot of Cell Segmentation.
- **Details:** Shows a tissue section with 13 spatial clusters mapped.

## Batch 3: Pages 11-15

### Page 11: Loupe Browser - Aggregation
- **Title:** Aggregation (Loupe Browser)
- **Comparison Table Summary:**
    - **Colorectal Cancer:**
        - Total Reads (single end): **721,219,716**
        - Mapping Rate (Reads mapped to probe set): **98.8%**
        - Uniquely Mapped Rate: **98.2%**
        - Cell Number (Segmentation): **220,703**
        - Median Genes per cell: **581**
        - Median UMI counts per cell: **759**
        - Total genes detected: **18,132**
    - **Lymph Node:**
        - Total Reads (single end): **443,996,591**
        - Mapping Rate (Reads mapped to probe set): **98.9%**
        - Uniquely Mapped Rate: **96.3%**
        - Cell Number (Segmentation): **295,776**
        - Median Genes per cell: **68**
        - Median UMI counts per cell: **72**
        - Total genes detected: **18,132**

### Page 12: Cluster Cell Count Distribution
- **Content:** Comparison table and chart of cell percentage across 29 clusters.
- **Key Findings:**
    - **Total Cells:** Colorectal Cancer = 220,694 (Note: table total shows 220,694, Page 11 table shows 220,703); Lymph Node = 295,776.
    - **Colorectal Cancer Largest Clusters:** Cluster 5 (15.787%, 34,841 cells), Cluster 6 (13.523%, 29,845 cells), Cluster 9 (9.655%, 21,307 cells), Cluster 10 (8.101%, 17,879 cells), Cluster 11 (7.128%, 15,731 cells).
    - **Lymph Node Largest Clusters:** Cluster 1 (20.337%, 60,153 cells), Cluster 2 (19.143%, 56,619 cells), Cluster 3 (18.430%, 54,511 cells), Cluster 4 (12.755%, 37,725 cells).
    - Other clusters (e.g., Clusters 5, 9, 11, 13, 14, 16, 19, 21, 24, 25, 26, 27, 29) have 0% or near-zero cell representation in Lymph Node, showing high sample specificity.

### Page 13: UMAP Cluster Visualization
- **Title:** Cluster Analysis - UMAP圖 (UMAP Plot)
- **Visuals:** Shows UMAP projection of all cells colored by sample (Colorectal Cancer vs Lymph Node) and colored by cluster, as well as split UMAPs per sample. Shows distinct separation between Colorectal Cancer cells and Lymph Node cells.

### Page 14: Spatial Mapping of Clusters
- **Title:** Aggregation Cluster
- **Visuals:** Overlays UMAP-identified clusters onto the physical H&E tissue slides for Lymph Node and Colorectal Cancer.
    - Lymph Node shows distinct follicular/germinal center structures mapped to specific clusters.
    - Colorectal Cancer shows complex tumor microenvironment zoning mapped to specific clusters.

### Page 15: Significant Genes Per Cluster
- **Title:** Top50 _ all significant genes per cluster
- **Content:** Large data table containing gene expression metrics (Average, Log2 Fold Change, P-value) for each cluster.
- **Top Genes Listed:**
    - **CHGB** (Chromogranin B): Feature ENSG00000089199, high average in Cluster 1 (79.84, Log2FC: 0.60) and Cluster 2 (70.94, Log2FC: 0.35).
    - **TFF1** (Trefoil Factor 1): Feature ENSG00000160182, high average in Cluster 1 (6.01, Log2FC: 0.72), Cluster 2 (5.03, Log2FC: 0.72), Cluster 3 (5.71, Log2FC: 0.72).
    - **GC** (Vitamin D Binding Protein): Feature ENSG00000145321, high average in Cluster 1 (15.84, Log2FC: 0.57).
    - Other key significant features: **CLEC3A**, **SLC18A2**, **MUC13**, **AGR3**, **CXCL13**, **AGR2**, **CHGA**, mitochondrial genes (**MT-ND4**, **MT-CYB**, **MT-CO3**, **MT-CO2**, **MT-ND2**), **SLPI**, **CGA**, **TFF3**, **HMGCS2**, **LEPR**, **TMEM176A**, **SCG2**, **PCSK1**, **RNASE4**, **TTR**.

## Batch 4: Pages 16-20

### Page 16: Feature Visualization - OSM
- **Title:** Selected features (UMAP)
- **Gene:** OSM (Oncostatin M) - Feature ID: feature_17127.
- **Visuals:**
    - UMAP plot shows localized expression in a small cluster.
    - Violin plot shows high expression specifically in Cluster 29.
- **Other Selected Features Listed:**
    - **TMSB4X**: feature_17428
    - **IGHM**: feature_12369
    - **VIM**: feature_8729
- **Violin Plots for Cluster 1:** Shows expression distribution of marker genes for Cluster 1 across all clusters.

### Page 17: Cell Annotation Resources
- **Title:** Cell annotation
- **Content:** Screenshot of **CellMarker 2.0** database interface (http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.html).
- **Database Statistics:**
    - **Human:** 429 Tissue types, 1,715 Cell types, 278 Cancer types, 4,334 Cell marker sets, 15,737 Cell markers.
    - **Mouse:** 399 Tissue types, 1,434 Cell types, 94 Cancer types, 2,185 Cell marker sets, 11,385 Cell markers.

### Page 18: Cell Type Assignment & Distribution
- **Content:** Detailed tables mapping clusters to cell types and summarizing cell type abundance.
- **Cluster to Cell Type Mapping:**
    - Clusters 1, 15 -> **Naive B cell**
    - Clusters 2, 3, 7 -> **Naive T cell**
    - Clusters 4, 28 -> **Memory B cell**
    - Cluster 5 -> **Transit-amplifying cell**
    - Cluster 6 -> **Macrophage**
    - Clusters 9, 25, 26 -> **Coloncyte**
    - Clusters 10, 23 -> **Plasma cell**
    - Clusters 11, 16 -> **Colonic stem cell**
    - Clusters 12, 13, 14 -> **Enteroendocrine cell**
    - Clusters 17, 21 -> **Cancer-associated fibroblast**
    - Cluster 18 -> **Low-density neutrophil**
    - Clusters 19, 27, 29 -> **Monocyte**
    - Cluster 20 -> **B cell**
    - Cluster 22 -> **Fibroblast**
    - Cluster 24 -> **Paneth cell**
- **Cell Type Distribution Summary (Key Samples):**
    - **Naive T cell:** 0.329% in Colorectal Cancer vs **46.585%** in Lymph Node.
    - **Naive B cell:** 0.036% in Colorectal Cancer vs **24.093%** in Lymph Node.
    - **Memory B cell:** 0.134% in Colorectal Cancer vs **13.019%** in Lymph Node.
    - **Enteroendocrine cell:** **18.306%** in Colorectal Cancer vs 0.041% in Lymph Node.
    - **Transit-amplifying cell:** **15.787%** in Colorectal Cancer vs 0.000% in Lymph Node.
    - **Macrophage:** **13.523%** in Colorectal Cancer vs 0.095% in Lymph Node.
    - **Coloncyte:** **12.612%** in Colorectal Cancer vs 0.000% in Lymph Node.

### Page 19: Example of Cell Type Changes & Subtype Analysis
- **Title:** Changes in cell types | Subtype clustering analysis
- **Content:** Academic reference figures (Zakharov PN et al., J Exp Med. 2020) demonstrating:
    - Dynamic changes in cell type proportions over time (4w, 8w, 15w).
    - Subclustering of T cells and Innate Lymphoid Cells (ILCs) into CD4, CD8, gamma delta T, ILC2, ILC3, NK, and NKT cells.
    - Violin plots for specific marker genes (Cd3e, Cd4, Cd8a, Gata3, Eomes, etc.).

### Page 20: Gene Expression Comparison Table
- **Title:** Comparison of Gene Expression
- **Content:** Screenshot of a large data table comparing expression metrics (Log2FoldChange, Means) between "Treat vs Control" across different clusters.
- **Example Genes:** SCNN1D, ACAP3, PUSL1, INTS11, CPTP, TAS1R3, DVL1, MXRA8, AURKAIP1, CCNL2.

## Batch 5: Pages 21-25

### Page 21: Pathway Enrichment Analysis (KEGG)
- **Title:** Pathway enrichment Analysis
- **Content:** Screenshot of Kyoto Encyclopedia of Genes and Genomes (**KEGG**) database interface (https://www.genome.jp/kegg/) and a table of enriched pathways.
- **Top Enriched Pathways (Sorted by P-value):**
    - **Olfactory transduction** (hsa04740): 11 DEGs (1.36%), 7.00E-10 P-value, 2.13E-07 FDR.
    - **Calcium signaling pathway** (hsa04020): 56 DEGs (6.91%), 6.62E-07 P-value, 0.010 FDR.
    - **Neuroactive ligand-receptor interaction** (hsa04080): 70 DEGs (8.63%), 8.51E-06 P-value, 0.072 FDR.
    - **Protein processing in endoplasmic reticulum** (hsa04141): 3 DEGs (0.37%), 9.48E-06 P-value, 0.476 FDR.
    - **Cell cycle** (hsa04110): 2 DEGs (0.25%), 8.55E-05 P-value, 0.476 FDR.
    - **Spliceosome** (hsa03040): 3 DEGs (0.37%), 9.37E-05 P-value, 0.476 FDR.
    - **Cholinergic synapse** (hsa04725): 28 DEGs (3.45%), 0.0003 P-value, 1.389 FDR.
    - **Insulin secretion** (hsa04911): 23 DEGs (2.84%), 0.0004 P-value, 1.389 FDR.
    - **Ribosome biogenesis in eukaryotes** (hsa03008): 2 DEGs (0.25%), 0.0004 P-value, 1.389 FDR.
    - **Autophagy - animal** (hsa04140): 4 DEGs (0.49%), 0.0007 P-value, 2.167 FDR.

### Page 22: Enriched Pathway Scatter & Volcano Plots
- **Title:** Enriched Pathway - Scatter Plot | Volcano plot
- **Visuals:**
    - **Scatter Plot (Statistics of Pathway Enrichment):** Displays top 20 enriched pathways by P-value. X-axis shows Rich Factor; dots are sized by gene number and colored by P-value.
    - **Volcano Plot:** Plots log2(Fold Change) vs -log10(P-value) for all genes. Red dots represent significantly up-regulated DEGs, green dots represent down-regulated DEGs, and grey dots represent non-significant genes (Not DEGs).

### Page 23: Expression Heatmap & Dotplot
- **Title:** Heatmap | Dotplot of features Analysis
- **Software:** Seurat
- **Visuals:**
    - **Heatmap:** Displays single-cell expression profiles of top marker genes across major cell types (Naive CD4 T, Memory CD4 T, CD14+ Mono, B, CD8 T, FCGR3A+ Mono, NK, DC, Platelet).
    - **Dotplot of Features Analysis:** Compares average expression levels and percentage of expressing cells for key features (Bmp4, Cdkn1b, Edn1, Fgfr3, Flna, Hspa1a, Hspa1b, Igf1r, Il1a, Insr, Klhdc8b, Lrp5, Pdgfb, Rab11a, Rhoa, Ube2s, Ywhah) between Treated vs Control states across distinct cell groups.

### Page 24: Trajectory & Pseudotime Analysis
- **Title:** Trajectory Analysis (Pseudotime)
- **Software:** Monocle2 package (https://cole-trapnell-lab.github.io/monocle-release/docs/)
- **Visuals:** Trajectory tree showing component 1 vs component 2:
    - Colored by time point hours (0, 24, 48, 72).
    - Colored by cellular differentiation state (1 to 5).
    - Colored by pseudotime timeline (0 to 25).

### Page 25: Cell-Cell Communication (CellChat)
- **Title:** CellChat
- **Concept:** Systematic analysis of cell-cell communication. CellChat quantifies the signaling communication probability between cell groups using a mass action-based model. It accounts for ligand-receptor complexes with multi-subunit structures and co-factors (agonists, antagonists, co-receptors).
- **CellChatDB Database Details:**
    - Contains ~3,300 Ligand-Receptor (L-R) pairs.
    - **Interaction Categories:** Secreted Signaling (12%), ECM-Receptor (38%), Cell-Cell Contact (20%), Non-protein Signaling (30%).
    - **Structure:** Heterodimers (41%), Others (59%).
- **Source Citation:** Nat Protoc 20, 180-219 (2025). URL: https://github.com/jinworks/CellChat.

## Batch 6: Pages 26-27

### Page 26: CellChat Analysis Examples
- **Title:** CellChat
- **Content:** Academic reference figures showing CellChat output analysis (Heliyon 10(15):e35263 (2024)).
- **Visuals:**
    - **Network Diagrams (Number of interactions):** Compares Control vs OA (Osteoarthritis) states, showing interaction density between cell groups (RegC, FC, preFC, HomC, RepC, EC, HTC).
    - **Pathway Networks:** Compares TGFb and FN1 signaling pathway networks between Control and OA.
    - **Role Heatmaps:** Identifies senders, receivers, mediators, and influencers for TGFb and FN1 signaling across cell groups.

### Page 27: Company Website & Contact Info
- **Title:** https://www.tri-ibiotech.com.tw/
- **Content:** Screenshot of the Tri-I Biotech website detailing NGS & Third Generation Sequencing technologies.
- **Key Mentions:**
    - Human Genome Project (HGP) using Sanger Sequencing.
    - Next Generation Sequencing (NGS) introduced in 2007 (Illumina system).
    - Third Generation Sequencing (Long-Read Sequencing) including PacBio (SMRT) in 2010 and Oxford Nanopore Technologies (ONT) in 2012.
    - Mention of Nanopore platforms: Flongle, MinION, GridION, and PromethION 2 Solo.
- **Contact:**
    - **Website QR Code** and **Official LINE QR Code**.
    - **Phone:** (02)-2695-4311, 0800-884311.

---





