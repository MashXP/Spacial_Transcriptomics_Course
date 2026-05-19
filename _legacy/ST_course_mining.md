
# Course 1: [Spatial Transcriptomics Course | Learn Spatial Gene Expression Analysis](https://www.arraygen.com/spatial_transcritome.php)

|                                    |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Module-NGS Spatial Transcriptomics | 📘 **Introduction to Spatial Transcriptomics**  <br>- What is spatial transcriptomics?  <br>- Comparison with bulk and single-cell RNA-seq  <br>- Key applications (cancer, neuroscience, immunology, etc.)  <br>- Overview of current platforms  <br>- Tissue sectioning and imaging  <br>- RNA capture and library preparation  <br>- Understanding barcoded arrays or in situ hybridization  <br>- Best practices to avoid batch effects and RNA degradation  <br>- Understanding NGS and different file formats  <br>  <br>📘 **Linux Basics & Environment Setup**  <br>- Linux Command Line Basics  <br>- Installing Spatial Transcriptomics Analysis Tools  <br>- Using Conda and Shell Scripting  <br>  <br>📘 **Introduction to R/Bioconductor**  <br>- Installing packages with CRAN and Bioconductor  <br>- Data types and standardized data container  <br>- Data manipulation  <br>  <br>📘 **Data Preprocessing & Quality Control**  <br>- Image alignment and registration  <br>- Running pipeline  <br>- Filtering low-quality spots  <br>- Visualization of QC metrics (UMIs, gene counts, mitochondrial genes)  <br>  <br>📘 **Downstream Analysis – Spatial Expression Mapping**  <br>- Loading data into Seurat or Scanpy  <br>- Normalization and scaling  <br>- PCA, UMAP/t-SNE for dimensionality reduction  <br>- Identifying spatial patterns of gene expression  <br>  <br>📘 **Clustering & Spatial Domain Detection**  <br>- Unsupervised clustering of spots  <br>- Spatially informed clustering  <br>- Identifying marker genes per region or domain  <br>- Integrating histological images with expression data  <br>  <br>📘 **Spatial Differential Expression**  <br>- Detecting spatially varying genes  <br>- Statistical models (SpatialDE, SPARK, etc.)  <br>- Visualizing expression gradients  <br>  <br>📘 **Cell–Cell Interaction & Ligand-Receptor Analysis**  <br>- Inferring spatial interactions between cells  <br>- Interpreting interaction networks  <br>  <br>📘 **Additional Post Analysis**  <br>- Different plots (Heatmap/volcano plot)  <br>- Functional annotation and pathway enrichment (clusterProfiler, Enrichr, KEGG)  <br>- Network analysis using STRING-db and Cytoscape for spatial DEGs |
# Course 2: [Spatial transcriptomics with R/Bioconductor - physalia-courses](https://www.physalia-courses.org/courses-workshops/spatial-transcriptomics/)
## Session content

**Day 1 The WHYs and the HOWs of spatial transcriptomics - 2-6 PM Berlin time**  
  
    Intro  
    Lecture 1 The spatial dimension  
    Lecture 2 Mapping transcripts in space  
    Lab 1 Computational set up & raw spatial data  
  
**Day 2 From transcripts to single cells****- 2-6 PM Berlin time**  
  
    Lecture 3 Imaging-based spatial data analysis  
    Lecture 4 Sequencing-based spatial data analysis  
    Lab 2 Cell segmentation & interactive data exploration  
  
**Day 3 Spatial mapping of cell types and cell states****- 2-6 PM Berlin time**  
  
    Lecture 5 Shared downstream data processing and QC  
    Lecture 6 Leveraging paired single-cell and spatial omics  
    Lab 3 Single-cell resolved spatial analyses  
  
**Day 4 The cellular niche****- 2-6 PM Berlin time**  
  
    Lecture 7 Cells live and interact in multicellular niches  
    Lab 4 Neighbourhood and niche analysis  
    Lab 5 Spatial cell-cell communication analysis  
  
**Day 5 A special  future****- 2-6 PM Berlin time**  
  
    Lecture 8 Towards 4D spatial multiomics  
    Workshop: Design your own spatial project  
    Hackathon: Start digging in your spatial data  
    Conclusion

# Course 3: [Introduction to Spatial Transcriptomics Data Analysis – Functional Genomics Center Zurich | ETH Zurich](https://fgcz.ch/education/bioinformaticsTraining/spatial_course.html)
## **Course Programme**

**General Overview of Spatial Transcriptomics**  
• Brief history of spatial transcriptomics  
• Overview of different categories of current spatial transcriptomics methods  
• Detailed introduction to Spatial Transcriptomics at FGCZ  
  -Strengths and limitations of different technologies  
• Lab tour for spatial transcriptomics instruments

  
**Introduction to Spatial Transcriptomics Analysis**  
• Data preprocessing and quality control for Visium HD and Xenium data  
• Spatial visualization of gene expression patterns  
• Cell type annotation using scRNA-seq reference mapping  
• Spatial clustering and region identification  
• Differential expression analysis in spatial context  
• Analysis of spatial relation between cell types  
• Hands-on analysis using Seurat V5 spatial workflows

  
**Optional: Own Data Analysis**  
• Consultation on participant-provided datasets  
• Discussion of experimental design considerations  
• Guided analysis of individual research questions  
• Troubleshooting common challenges in spatial data  
• Implementation of customized analysis workflows  
• Integration of spatial data with other omics data types

**Learning outcomes**  
At the end of the course, participants will be able to:  
• Summarize the principles of current spatial transcriptomics methods and data analysis concepts  
• Critically evaluate spatial transcriptomics data and analysis results  
• Apply fundamental visualization and interpretation techniques to spatially resolved transcriptomics data  
• Discuss how spatial transcriptomics can be integrated into a research experiment  
• Make an informed decision, which Spatial Transcriptomics technology to apply to a specific research question

# Course 4: [Spatial transcriptomics data analysis: theory and practice](https://bookdown.org/sjcockell/ismb-tutorial-2023/index.html)

[[mining_toc]]
## Learning objectives
This is very good resource


Participants in this tutorial will gain understanding of the core technologies for undertaking a spatial transcriptomics experiment, and the common tools used for the analysis of this data. In particular, participants will appreciate the strengths of geospatial data analysis methods in relation to this type of data. Specific learning objectives will include:

1. Describe and discuss core technologies for spatial transcriptomics
2. Make use of key computational technologies to process and analyse STx data
3. Apply an analysis strategy to obtain derived results and data visualisations
4. Appreciate the principles underlying spatial data analysis
5. Understand some of the methods available for spatial data analysis
6. Apply said methods to an example STx data set