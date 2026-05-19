# Spatial Transcriptomics Analysis: From Wet Lab to Advanced Geocomputational Modeling
## An End-to-End Curriculum (Theory & Practice)

This curriculum provides a structured, step-by-step pathway through spatial transcriptomics (STx), starting from the wet-lab bench and proceeding through upstream processing, downstream bioinformatics, and advanced spatial statistics. Each module is divided into **Theory** (conceptual, physical, and mathematical principles) and **Practice** (protocol steps, command-line arguments, and R/Bioconductor implementations).

---

## Program Workflow Overview

![[workflow_diagram.jpg]]

---

![[phase0.jpg]]
## Module 0: Foundational Prerequisites (Genomics, S4 Objects, and Spatial Math)

### I. Theory

#### A. Molecular Biology & Sequencing Chemistry
*   **Transcript Structure & Capture**: Poly-adenylated mRNA transcripts carry a 3' poly(A) tail. We exploit this tail by using complementary poly(dT) oligonucleotide sequences anchored to a solid substrate (slide or bead) to hybridize and capture mRNA during reverse transcription (RT).
*   **Reverse Transcription Kinetics**: Turning RNA into single-stranded complementary DNA (cDNA) using a reverse transcriptase enzyme. For spatial transcriptomics, this reaction occurs *in situ* (directly inside the permeabilized tissue section).
*   **Illumina Library Architecture**:
    *   **P5 / P7 Adapters**: Terminal sequences that hybridize to the flow cell oligos for bridge amplification.
    *   **Sample Indexes (i5 / i7)**: Barcodes that label each pooled library, allowing multiple samples to run on a single flow cell lane (multiplexing).
    *   **Read 1 (R1)**: Sequences the spatial barcode (giving the coordinates of the transcript on the slide) and the Unique Molecular Identifier (UMI).
    *   **Read 2 (R2)**: Sequences the cDNA insert (the biological transcript itself).
*   **UMI Math vs. Raw Counts**: PCR amplification introduces duplicate sequences. To count *molecules* rather than PCR duplicates, each captured transcript receives a random 10–12 base pair sequence (UMI). After sequencing, reads with identical spatial barcodes, gene mappings, and UMIs are collapsed into a single count.
    
$$\text{True Expression Count} = \text{Count of unique [Spatial Barcode + UMI + Gene ID]} \text{ combinations}$$

#### B. R/Bioconductor S4 Object Architecture
*   **The S4 Class System**: Unlike basic R objects (lists, vectors), Bioconductor relies on S4 objects, which enforce strict class structures, validation rules, and inheritance chains.
*   **SummarizedExperiment Inheritance**:
    *   `SummarizedExperiment`: The base class containing a matrix of counts (`assays`), row metadata (`rowData` for genes), and column metadata (`colData` for samples/spots).
    *   `SingleCellExperiment`: Extends the base class by adding slots for dimensionality reductions (`reducedDims`) and alternative experiments (`altExps`).
    *   `SpatialExperiment`: Extends `SingleCellExperiment` by adding coordinate-specific slots (`spatialCoords`) and tissue image tracking (`imgData`).
*   **Memory Footprint of Sparse Matrices**: Spatial datasets can contain millions of zero values (sparsity >90%). Storing this as a standard dense matrix wastes system memory. We use Compressed Sparse Column (`dgCMatrix`) formats, storing only non-zero values and their indices.

#### C. 2D Coordinate Geometry & Tessellation
*   **Cartesian Grids**: Every spatial transcriptomics spot is registered onto a 2D plane with ($x,y$) pixel coordinates.
*   **Euclidean Distance**: Distance between two spots $i$ and $j$:

$$d_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}$$

*   **Tessellation Patterns**: Visium utilizes a hexagonal grid pattern (where every spot has 6 equidistant neighbors), whereas other platforms use square grid arrangements (4 orthogonal, 4 diagonal neighbors). Hexagonal packing maximizes space coverage efficiency.

---

### II. Practice

#### A. Working with Sparse Matrices in R
In this section, we create a small sparse matrix to demonstrate how they save computer memory compared to standard dense matrices.

```R
library(Matrix)

# 1. Create a dense matrix with lots of zeros
dense_mat <- matrix(0, nrow = 5, ncol = 5)
dense_mat[1, 2] <- 4
dense_mat[3, 5] <- 10
dense_mat[5, 1] <- 1

print(dense_mat)

# 2. Convert to dgCMatrix (Compressed Sparse Column)
sparse_mat <- as(dense_mat, "dgCMatrix")
print(sparse_mat)

# 3. Check memory size comparison (in bytes)
object.size(dense_mat)
object.size(sparse_mat)
```

#### B. Manually Constructing a SummarizedExperiment Object
Understanding S4 structure by assembling an object from its component matrices:

```R
library(SummarizedExperiment)

# 1. Create dummy assay data (counts)
counts_matrix <- matrix(rpois(50, lambda = 2), nrow = 10, ncol = 5)
rownames(counts_matrix) <- paste0("GENE_", 1:10)
colnames(counts_matrix) <- paste0("SPOT_", 1:5)

# 2. Create row metadata (Gene characteristics)
gene_metadata <- DataFrame(
  gene_id = paste0("ENSG00000", 1:10),
  biotype = rep(c("protein_coding", "lncRNA"), each = 5)
)
rownames(gene_metadata) <- rownames(counts_matrix)

# 3. Create column metadata (Spot characteristics)
spot_metadata <- DataFrame(
  barcode = c("AAAC", "AAAG", "AAAT", "AAAC", "AAAG"),
  condition = c("Control", "Control", "Treated", "Treated", "Treated")
)
rownames(spot_metadata) <- colnames(counts_matrix)

# 4. Construct the S4 object
se <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  rowData = gene_metadata,
  colData = spot_metadata
)

# 5. Extract components using accessor functions
head(assay(se, "counts"))
head(colData(se))
head(rowData(se))
```

---

## Phase I: Wet Lab (Sample Prep to Sequencing)
![[phase1.jpg]]

### Module 1: Tissue Handling, Quality Control, and Staining

#### I. Theory
*   **Tissue Preservation Dynamics**: Freezing rate impact on cellular morphology. OCT (Optimal Cutting Temperature) compound polymerization. Why direct liquid nitrogen immersion causes tissue cracking (due to Leidenfrost effect) and the use of dry-ice cooled isopentane.
*   **Substrate Chemistry**: Charge interactions on Matsunami Platinum Pro slides vs. Schott Nexterion Slide H (functionalized for difficult-to-adhere tissue types).
*   **Pre-analytical RNA Integrity**: Theoretical principles of RNA degradation in tissue blocks. RIN (RNA Integrity Number) calculation based on 28S/18S rRNA ratios vs. DV200 (percentage of RNA fragments >200 nt) as a robust metric for FFPE samples.

#### II. Practice
*   **Cryosectioning & Mounting**:
    *   Operating the cryostat at $-20^\circ\text{C}$ to $-15^\circ\text{C}$.
    *   Cutting sections at $10\ \mu\text{m}$ (Fresh Frozen) or $5\ \mu\text{m}$ (FFPE).
    *   Adhering sections to Visium slides within the designated $6.5 \times 6.5\ \text{mm}$ capture area.
*   **H&E Staining Protocol**:
    *   Fixing sections in cold methanol (Fresh Frozen) or baking and deparaffinizing (FFPE).
    *   Applying Gill II Hematoxylin (nuclear stain), bluing buffer, and alcoholic Eosin (cytoplasmic stain).
*   **Imaging & Decoverslipping**:
    *   High-resolution tile-scan imaging using $20\times$ or $40\times$ brightfield magnification.
    *   Decoverslipping in $1\times\text{PBS}$ at room temperature to avoid peeling tissue.
    *   Destaining using a light hydrochloric acid / ethanol wash sequence.

#### Related Resources
*   [[Visium_HD_Workflow#Page 2: Product Introduction|Visium_HD_Workflow.md: Sectioning & Slide Preparation]]
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 7: Workflow & QC Requirements|Spatial_Transcriptomics_Sequencing_Service.md: Tissue and Specimen QC Requirements]]

---

### Module 2: CytAssist Biochemistry & Capture Chemistry

#### I. Theory
*   **CytAssist Capture Mechanism**: Using a transition instrument to bridge paraffin-embedded or fresh-frozen sections from standard glass slides to spatial capture slides.
*   **Oligo Capture Chemistry**: Design of the Visium capture probes. Spatial barcodes (18 nt) providing geographical coordinates ($x,y$), Unique Molecular Identifiers (UMIs, 12 nt) to resolve PCR duplication, and poly(dT)VN tails to capture poly-adenylated transcripts.
*   **Permeabilization & Diffusion Kinetics**: Balancing enzyme concentration and incubation time to release transcripts without causing lateral diffusion (which blurs spatial resolution).

```
Spatial Capture Barcode Primer Design:
5'- [Slide Surface] - Spacer - Read 1 Linker - Spatial Barcode (18bp) - UMI (12bp) - Poly(dT)VN - 3'
```

#### II. Practice
*   **CytAssist Run Execution**:
    *   Aligning the tissue slide with the Visium slide inside the CytAssist alignment cassette.
    *   Adjusting the alignment guide parameters (Cassette size: $6.5\ \text{mm}$ or $11\ \text{mm}$).
    *   Executing the CytAssist thermal cycle program.
*   **In Situ cDNA Synthesis**:
    *   Running on-slide reverse transcription using template-switching oligos (TSO).
    *   Applying Second Strand Synthesis reagents to generate double-stranded cDNA.
    *   Performing on-slide chemical elution using potassium hydroxide ($\text{KOH}$) and neutralisation with $1\ \text{M}$ Tris-HCl (pH 7.0).

#### Related Resources
*   [[Visium_HD_Workflow#Page 8 (Slide 9): Visium HD Cassette Kit|Visium_HD_Workflow.md: CytAssist Cassette Setup & Run Details]]
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 4: Workflow Overview|Spatial_Transcriptomics_Sequencing_Service.md: CytAssist Transfer Overview]]

---

### Module 3: Library Prep & Sequencing Configuration

#### I. Theory
*   **cDNA Fragmentation Mechanics**: Physical shear vs. enzymatic fragmentation. Standardizing insert sizes to fit the Illumina sequencing window.
*   **Library Amplification Principles**: Polymerase bias during index incorporation. Calculating PCR cycle counts as a function of cDNA yield to minimize PCR duplicates while obtaining sufficient library concentration.
*   **Sequencing Metric Math**: Calculating required read depth based on tissue coverage.

$$\text{Required Reads} = \text{Capture Area Coverage \%} \times \text{Minimum Platform Read Depth}$$

*   For **Visium HD 3'**, the minimum depth is $550\text{M}$ read pairs per fully covered array.

#### II. Practice
*   **cDNA Clean-up and QC**:
    *   Performing SPRIselect bead clean-ups at $0.6\times$ (right-side selection) and $0.8\times$ (left-side selection).
    *   Assessing cDNA fragment size distribution via TapeStation D5000 (expected peak: $500\text{--}1000\ \text{bp}$).
*   **Adapter Ligation & Indexing**:
    *   Enzymatic fragmentation, end repair, and A-tailing.
    *   Ligation of TruSeq Read 2 adapters.
    *   Executing Sample Index PCR using unique dual index (UDI) plates.
*   **Final Library QC & Sequencing Setup**:
    *   Validating index PCR output (expected peak: $\sim 450\ \text{bp}$).
    *   Quantifying library concentrations via KAPA Library Quantification qPCR.
    *   Configuring Illumina sequencing parameters:
        *   **Read 1**: $43$ cycles (captures spatial barcode and UMI).
        *   **Index 1 (i7)**: $10$ cycles.
        *   **Index 2 (i5)**: $10$ cycles.
        *   **Read 2**: $75$ cycles (captures cDNA insert).

#### Related Resources
*   [[Visium_HD_Workflow#Page 10 (Slide 11): Reagent Kits|Visium_HD_Workflow.md: Library Construction Reagents & Thermal Cycling Parameters]]
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 7: Workflow & QC Requirements|Spatial_Transcriptomics_Sequencing_Service.md: NovaSeq Sequencing Setup and Platforms]]

---

## Phase II: Upstream Processing
![[phase2.jpg]]
### Module 4: Space Ranger Pipeline Processing

#### I. Theory
*   **Alignment Algorithms**: Mapping reads to the reference transcriptome (e.g., STAR aligner).
*   **Barcode Demultiplexing**: Error-correction tolerances for spatial barcodes (allowing up to 1 mismatch).
*   **Image Alignment & Registration**: Cross-correlation algorithms mapping histological coordinates to the spatial grid of the array using fiducial markers.
*   **Binning Paradigms**: Summarizing sub-single cell transcripts into spatial squares. Visium HD maps raw coordinates into $2\ \mu\text{m}$ spatial bins, and pools them into $8\ \mu\text{m}$ and $16\ \mu\text{m}$ matrices.

#### II. Practice
*   **Running the Pipeline**:
    *   Executing `spaceranger count` via command line:
        ```bash
        spaceranger count \
          --id=sample_JBO019 \
          --transcriptome=/path/to/ref/GRCh38 \
          --probe-set=/path/to/probe_set.csv \
          --fastqs=/path/to/fastqs \
          --image=/path/to/histology_image.tiff \
          --slide=V12D34-567 \
          --area=A1 \
          --loupe-alignment=alignment.json \
          --create-bam=false
        ```
*   **QC Metrics Assessment**:
    *   Opening `web_summary.html`.
    *   Checking "Q30 Bases in Barcode" and "Q30 Bases in RNA Read" ($\ge 90\%$).
    *   Inspecting alignment logs to verify tissue detection and spot registration overlay.

#### Related Resources
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 9: Space Ranger Summary Report - Metrics|Spatial_Transcriptomics_Sequencing_Service.md: Space Ranger Web Summary Metric Definitions]]

---

### Module 5: Interactive Visual Exploration

#### I. Theory
*   **Unsupervised Spatial Clustering Principles**: K-means and Graph-based clustering models. How visual layouts assist pathobiology review.
*   **Marker Gene Selection Concepts**: Differential expression calculation on spatial clusters to distinguish anatomical layers vs. cell-type-specific boundaries.

#### II. Practice
*   **Loupe Browser Interface**:
    *   Opening the `.cloupe` output file.
    *   Adjusting tissue image opacity, zoom, and spatial coordinate grids.
*   **Interactive Spot Partitioning**:
    *   Manually annotating histological structures using the lasso tool.
    *   Exporting custom spot annotations as `.csv` metadata files.
    *   Running differential expression inside Loupe to identify local tissue layer marker genes.

#### Related Resources
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 11: Loupe Browser - Aggregation|Spatial_Transcriptomics_Sequencing_Service.md: Loupe Cluster Visualizations and Differential Genes]]

---

## Phase III: Downstream Bioinformatics (R/Bioconductor)
![[phase3.jpg]]

### Module 6: Data Classes, Quality Control, and Normalisation

#### I. Theory
*   **Computational Classes**: Structure of the Bioconductor `SpatialExperiment` (SPE) container. Assays (counts, logcounts), column data (spot metadata, `in_tissue` status), row data (gene metadata), and spatial coordinates.
*   **QC Biology**: Biological significance of UMI sum dropouts, low expressed gene counts, and high mitochondrial percentages (indicative of cellular membrane rupture and leakage of cytoplasmic mRNA).
*   **Normalisation Mechanics**: Why simple library-size scaling is preferred in STx over complex scRNA-seq methods (spots containing multiple heterogeneous cell types). Log-transformation variance stabilization.

```
SpatialExperiment Structure:
+-------------------------------------------------------------+
|                     SpatialExperiment                       |
+-------------------------------------------------------------+
|  assays()       --> counts, logcounts                       |
|  colData()      --> spot coordinates, metadata, sum, mito%  |
|  rowData()      --> gene name, ID, biological variance      |
|  spatialCoords()--> x, y pixel coordinates                  |
|  imgData()      --> lowres/hires histology images            |
+-------------------------------------------------------------+
```

#### II. Practice
*   **Data Ingestion & QC Filtering**:
    ```R
    library(SpatialExperiment)
    library(scater)
    library(scran)

    # 1. Loading the SpatialExperiment dataset
    spe <- read10xVisium(samples = "outs_dir", type = "sparse", data = "filtered")

    # 2. Calculating mitochondrial proportion
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

    # 3. Setting threshold filters based on distribution plots
    qc_lib_size <- colData(spe)$sum < 700
    qc_detected <- colData(spe)$detected < 500
    qc_mito <- colData(spe)$subsets_mito_percent > 28
    
    # 4. Filter spots
    spe$discard <- qc_lib_size | qc_detected | qc_mito
    spe <- spe[, !spe$discard]
    ```
*   **Normalisation**:
    ```R
    # Calculate size factors and compute log-normalized counts
    spe <- computeLibraryFactors(spe)
    spe <- logNormCounts(spe)
    ```

#### Related Resources
*   [[practical-session-1#1.4.1 SpatialExperiment class|practical-session-1.md: Creating and Navigating the SpatialExperiment Object]]
*   [[practical-session-2#2.1 Spot-level Quality Control|practical-session-2.md: Diagnostic Plots and Spot-level QC Implementation]]

---

### Module 7: Feature Selection & Spatially Variable Genes (SVGs)

#### I. Theory
*   **Highly Variable Genes (HVGs)**: Decomposition of variance into biological and technical components using fitted mean-variance trends across spots.
*   **Spatially Variable Genes (SVGs)**: Concept of spatial autocorrelation (Tobler's First Law of Geography: near things are more related than distant things).
*   **Autocorrelation Metrics**:
    *   *Moran's $I$*: Measure of global spatial clustering (ranging from $-1$ to $+1$).
    *   *Geary's $C$*: Measure of local spatial dissimilarity (ranging from $0$ to $2$).
    *   *SPARK-X*: Non-parametric test for spatial expression patterns robust to cellular resolution scale.

#### II. Practice
*   **HVG Selection**:
    ```R
    # Fit mean-variance trend
    dec <- modelGeneVar(spe)
    
    # Select top 10% of biologically variable genes
    top_hvgs <- getTopHVGs(dec, prop = 0.1)
    ```
*   **SVG Autocorrelation Modeling**:
    ```R
    # Run SPARK-X or Moran's I on the logcounts matrix using spatial coordinates
    # Identify genes with significant spatial patterns (FDR < 0.05)
    # Intersect top HVGs and SVGs to compile a robust feature set for clustering
    ```

#### Related Resources
*   [[practical-session-2#2.3 Selecting genes|practical-session-2.md: Gene Feature Selection & HVG Trend Plotting]]

---

### Module 8: Dimensionality Reduction & Unsupervised Clustering

#### I. Theory
*   **Dimensionality Reduction**: The curse of dimensionality in high-attribute spaces. PCA (linear, preserves global distances) vs. UMAP (non-linear, preserves local manifolds for visualization).
*   **Graph-based Clustering**: Louvain/Walktrap algorithm execution on Shared Nearest Neighbor (SNN) graphs built from the top Principal Components (PCs).
*   **Spatial Clustering vs. Non-Spatial**: Incorporating coordinates into clustering algorithms to identify spatially contiguous tissues.

#### II. Practice
*   **PCA and UMAP Execution**:
    ```R
    # Set seed for reproducibility
    set.seed(987)
    
    # Run PCA on Selected HVGs
    spe <- runPCA(spe, subset_row = top_hvgs)
    
    # Run UMAP on PCA components
    spe <- runUMAP(spe, dimred = "PCA")
    ```
*   **Graph-Based Clustering**:
    ```R
    # Build SNN graph
    g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
    
    # Cluster using Walktrap
    g_walk <- igraph::cluster_walktrap(g)
    colLabels(spe) <- factor(g_walk$membership)
    ```
*   **Differential Marker Detection**:
    ```R
    # Identify cluster-specific marker genes
    markers <- findMarkers(spe, test = "binom", direction = "up")
    
    # Plot top marker gene log-fold changes as heatmap
    library(pheatmap)
    logFCs <- getMarkerEffects(markers[[1]])
    pheatmap(logFCs[1:20, ])
    ```

#### Related Resources
*   [[practical-session-2#2.4 Dimensionality reduction|practical-session-2.md: Dimensionality Reduction & Graph-Based Clustering Implementation]]

---

## Phase IV: Advanced Integration
![[phase4.jpg]]
### Module 9: Single-Cell Reference Integration & Deconvolution

#### I. Theory
*   **Deconvolution Math**: Reconstructing bulk-like spot expression matrices as linear combinations of single-cell reference signatures.
*   **Algorithms**:
    *   *RCTD (Robust Cell Type Decomposition)*: Uses a Poisson log-linear model to fit cell type profiles while accounting for platform-specific effects.
    *   *Cell2location*: A Bayesian hierarchical model estimating cell-type densities across spatial coordinate arrays.

#### II. Practice
*   **Running RCTD in R**:
    *   Loading the single-cell reference `SingleCellExperiment` object.
    *   Specifying cell-type annotation labels.
    *   Formatting spatial data to match RCTD constraints:
        ```R
        library(spacexr)
        
        # Create Reference & Spatial object
        query <- SpatialObject(spe@assays@data$counts, spatialCoords(spe))
        reference <- Reference(sc_counts, sc_cell_types)
        
        # Run RCTD in doublet or multi-cell mode
        myRCTD <- create.RCTD(query, reference, max_cores = 4)
        myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
        results <- myRCTD@results
        ```

#### Related Resources
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 17: Cell Annotation Resources|Spatial_Transcriptomics_Sequencing_Service.md: Cell Annotation Tools and Deconvolution References]]

---

### Module 10: Trajectories & Cell-Cell Communication

#### I. Theory
*   **Spatial Pseudotime**: Tracking developmental trajectories along spatial gradients rather than temporal coordinates.
*   **Spatial Cell-Cell Communication**: Modeling signaling pathways across tissue domains. Integrating distance decay limits to represent physical limits on paracrine and juxtacrine signaling.

```
Spatial Signaling Constraint:
P_ij = L_i * R_j * exp(-d_ij / d_max)
Where P_ij is interaction probability, L_i is Ligand expression at spot i, 
R_j is Receptor expression at spot j, d_ij is spatial distance, d_max is decay limit.
```

#### II. Practice
*   **Spatial Trajectories with Monocle3**:
    *   Importing `SpatialExperiment` assays into a Monocle3 cell dataset.
    *   Running path learning algorithms to construct trajectory graphs across coordinates.
*   **CellChat Spatial Modeling**:
    *   Creating a CellChat object from spatial coordinates and normalized counts.
    *   Applying the `CellChatDB` ligand-receptor database:
        ```R
        library(CellChat)
        cellchat <- createCellChat(object = spe, assay = "logcounts", datatype = "spatial")
        cellchat <- identifyOverexpressedGenes(cellchat)
        cellchat <- computeCommunProb(cellchat, type = "truncatedMean", distance.use = TRUE)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        ```

#### Related Resources
*   [[Spatial_Transcriptomics_Sequencing_Service#Page 25: Cell-Cell Communication (CellChat)|Spatial_Transcriptomics_Sequencing_Service.md: CellChat Communication Probabilities & Role Analysis]]

---

## Phase V: Geocomputational Statistics & Geographically Weighted Models
![[phase5.jpg]]
### Module 11: Geocomputational Structures, Neighbor Graphs, and Spatial Weights

#### I. Theory
*   **Simple Features (`sf`) Specifications**: Representing spatial boundaries. `sf` integrates standard geometries (`sfg`) and coordinate systems (`sfc`) as columns inside attribute tables.
*   **Neighbour Connectivity Models**:
    *   *Adjacency-based*: Contiguity (Rook vs. Queen configurations).
    *   *Distance-based*: K-nearest neighbors ($k$-NN) or threshold distance bands.
*   **Spatial Weight Standardization**:
    *   Row-standardised weights (**W**): Normalizing neighbor counts so weights sum to $1$ for each location, preventing spots with high neighbor density from dominating models.

#### II. Practice
*   **Generating sf and Graphs in SFE**:
    ```R
    library(SpatialFeatureExperiment)
    library(spdep)

    # 1. Loading SFE container (which incorporates sf geometries)
    sfe <- read10xVisiumSFE(samples = "JBO019_Results", type = "sparse")

    # 2. Visualizing tissue annotation polygons
    ggplot() + 
      geom_sf(data = colGeometry(sfe, "spotHex"), aes(fill = colData(sfe)$annotation)) + 
      theme_void()

    # 3. Constructing k-nearest neighbor graphs with W standardization
    sfe <- addSpatialNeighGraphs(sfe, sample_id = "JBO019", type = "knearneigh", style = "W", k = 6)
    ```

#### Related Resources
*   [[practical-session-3#3.2 Background|practical-session-3.md: Geocomputational simple features & spatial neighbor configuration]]

---

### Module 12: Geographically Weighted Principal Components Analysis (GWPCA)

#### I. Theory
*   **Global vs. Local Models**: Global methods (e.g., standard PCA) assume gene-to-PC loading ratios are stationary across a tissue. Local methods (e.g., GWPCA) allow these loadings to vary geographically.
*   **Spatial Kernel Function**: Mathematical weights applied based on distance from location $i$.
*   **Percentage of Total Variation (PTV)**: Local measures of how much variance is explained by the local PCs, revealing coordinates of high molecular heterogeneity.

#### II. Practice
*   **Local GWPCA Execution**:
    ```R
    library(GWmodel)

    # 1. Extract coordinates and local expression matrix from SFE
    coords <- spatialCoords(sfe)
    expr_matrix <- t(assays(sfe)$logcounts[top_hvgs[1:100], ]) # Using top 100 features

    # 2. Run geographically weighted PCA using a bi-square kernel
    # (Adjust bandwidth and run local matrix decompositions)
    # gwpca_results <- gwpca(data = expr_matrix, coords = coords, bw = 50, kernel = "bisquare")
    ```
*   **Evaluating Local Loadings**:
    *   Identifying leading genes (highest local loading value) for each coordinate point.
    *   Mapping PTV scores across tissue structures to highlight spatial discrepancies.

#### Related Resources
*   [[practical-session-4#4.1 Geographically Weighted Principal Components Analysis (GWPCA)|practical-session-4.md: Local PC Calibrations & Leading Gene Interpretations]]

---

## Summary of Resources

For comprehensive reference material, check the underlying documentation source files directly:

*   **Wet-Lab Protocols & Platforms**:
    *   [[Visium_HD_Workflow|Visium_HD_Workflow.md]]
        *   [[Visium_HD_Workflow#Page 2: Product Introduction|Product & Workflow Overview]]
        *   [[Visium_HD_Workflow#Page 6 (Slide 6): Components, Shipping, and Storage|Kit Storage & Shipping Temperatures]]
        *   [[Visium_HD_Workflow#Page 17 (Slide 20): Tested Glass Slides|Glass Slide & Tissue QC Allowances]]
        *   [[Visium_HD_Workflow#Page 36 (Slide 41): Step 3: Poly(A) RNA Capture, RT, & Denaturation (Capture)|On-Slide RT, Denaturation & Elution Protocols]]
        *   [[Visium_HD_Workflow#Page 41 (Slide 46): Step 6: Library Construction (Fragmentation, End Repair, & A-tailing)|Library Fragmentation & Sample Indexing Cycles]]
    *   [[Spatial_Transcriptomics_Sequencing_Service|Spatial_Transcriptomics_Sequencing_Service.md]]
        *   [[Spatial_Transcriptomics_Sequencing_Service#Page 6: Comparison of Visium Types|Visium Platform Comparisons & Resolutions]]
        *   [[Spatial_Transcriptomics_Sequencing_Service#Page 7: Workflow & QC Requirements|Sequencing QC Requirements & Platforms]]
        *   [[Spatial_Transcriptomics_Sequencing_Service#Page 9: Space Ranger Summary Report - Metrics|Space Ranger Summary Metrics]]
        *   [[Spatial_Transcriptomics_Sequencing_Service#Page 11: Loupe Browser - Aggregation|Loupe Browser Cluster Interpretations]]
        *   [[Spatial_Transcriptomics_Sequencing_Service#Page 25: Cell-Cell Communication (CellChat)|CellChat Ligand-Receptor Communications]]
*   **R/Bioconductor Computational Code**:
    *   [[practical-session-1|practical-session-1.md]]
        *   [[practical-session-1#1.4.1 SpatialExperiment class|SpatialExperiment Container Architecture]]
        *   [[practical-session-1#1.3 Import 10X Visium data|Loading 10x Space Ranger Outputs]]
    *   [[practical-session-2|practical-session-2.md]]
        *   [[practical-session-2#2.1 Spot-level Quality Control|Spot-Level Quality Control & Filtering]]
        *   [[practical-session-2#2.2 Normalisation of counts|Log-Transformation Normalisation]]
        *   [[practical-session-2#2.3 Selecting genes|Highly Variable Genes (HVG) Selection]]
        *   [[practical-session-2#2.4 Dimensionality reduction|PCA, UMAP & Walktrap Clustering]]
        *   [[practical-session-2#2.5 Clustering|Binomial Marker Gene Selection]]
    *   [[practical-session-3|practical-session-3.md]]
        *   [[practical-session-3#3.3 Data structures preparation|Simple Features (sf) and Hexagonal Geometries]]
        *   [[practical-session-3#3.8 Neighbour graph and distance matrix|Neighbor Contiguities & Graph Construction]]
    *   [[practical-session-4|practical-session-4.md]]
        *   [[practical-session-4#4.1 Geographically Weighted Principal Components Analysis (GWPCA)|Geographically Weighted PCA (GWPCA) Concepts]]
        *   [[practical-session-4#4.5 Run GWPCA|Local Model PCA Execution]]
        *   [[practical-session-4#4.7 Identify the leading genes in each location|Leading Genes Extraction & PTV Mapping]]
