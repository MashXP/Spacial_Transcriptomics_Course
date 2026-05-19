

<!-- PAGE: index.html -->

# Spatial transcriptomics data analysis: theory and practice

_Eleftherios Zormpas, Dr Simon J. Cockell_

_2023-07-20_

# Welcome[](index.html#welcome)

This book will guide you through the practical steps of the in-person tutorial IP2 for the ISMB/ECCB 2023 conference in Lyon named: _âSpatial transcriptomics data analysis: theory and practiceâ_.

## Abstract[](index.html#abstract)

Recent technological advances have led to the application of RNA Sequencing _in situ_. This allows for whole-transcriptome characterisation, at approaching single-cell resolution, while retaining the spatial information inherent in the intact tissue. Since tissues are congregations of intercommunicating cells, identifying local and global patterns of spatial association is imperative to elucidate the processes which underlie tissue function. Performing spatial data analysis requires particular considerations of the distinct properties of data with a spatial dimension, which gives rise to an association with a different set of _statistical_ and _inferential_ considerations.

In this comprehensive tutorial, we will introduce users to spatial transcriptomics (STx) technologies and current pipelines of STx data analysis inside the **Bioconductor** framework. Furthermore, we will introduce attendees to the underlying features of spatial data analysis and how they can effectively utilise space to extract in-depth information from STx datasets.

## Learning objectives[](index.html#learning-objectives)

Participants in this tutorial will gain understanding of the core technologies for undertaking a spatial transcriptomics experiment, and the common tools used for the analysis of this data. In particular, participants will appreciate the strengths of geospatial data analysis methods in relation to this type of data. Specific learning objectives will include:

  1. Describe and discuss core technologies for spatial transcriptomics
  2. Make use of key computational technologies to process and analyse STx data
  3. Apply an analysis strategy to obtain derived results and data visualisations
  4. Appreciate the principles underlying spatial data analysis
  5. Understand some of the methods available for spatial data analysis
  6. Apply said methods to an example STx data set




<!-- PAGE: practical-session-1.html -->

# Chapter 1 Practical session 1[](practical-session-1.html#practical-session-1)

In this practical session you will familiarise yourself with some example spatial transcriptomics (STx) data and the common features of such data.

## 1.1 About this documentation[](practical-session-1.html#about-this-documentation)

This handbook is designed to walk you through the practical elements of todayâs tutorial. All of the code you need to accomplish the basic tasks thoughout the day is presented in full, there are some âstretch goalsâ in some of the tutorials where only a template is provided. This is not a typing tutorial, so feel free to copy and paste where necessary. The tutorials are written in `rmarkdown` and presented on [bookdown.org](https://bookdown.org/).

## 1.2 Posit Cloud[](practical-session-1.html#posit-cloud)

You should have received an invite to the [Posit Cloud](https://posit.cloud/) Space for todayâs tutorial. Accepting this invite will give you access to the 4 RStudio projects for the 4 sessions we will run today. Each project has the required packages pre-installed, and the data files uploaded. These projects are set up as âAssignmentsâ so that you get your own copy of the workspace.

## 1.3 Import 10X Visium data[](practical-session-1.html#import-10x-visium-data)

In this tutorial we will be using data from the [STexampleData](https://bioconductor.org/packages/STexampleData) package that contains a small collection of STx datasets from different technologies, including SlideSeq V2, seqFISH and 10x Genomics Visium. These datasets are provided in the `SpatialExperiment` format - described below.

The specific dataset used for this tutorial is a single sample from the dorsolateral prefrontal cortex (DLPFC) - a 10x Genomics Visium dataset that was published by Maynard et al. (2021).

Here, we show how to load the data from the `STexampleData` package.
    
    
    [](practical-session-1.html#cb1-1)library(SpatialExperiment)
    [](practical-session-1.html#cb1-2)library(STexampleData)
    [](practical-session-1.html#cb1-3)library(ggplot2)
    [](practical-session-1.html#cb1-4)library(ggspavis)
    [](practical-session-1.html#cb1-5)
    [](practical-session-1.html#cb1-6)# Load the object
    [](practical-session-1.html#cb1-7)spe <- Visium_humanDLPFC()

## 1.4 Explore data types[](practical-session-1.html#explore-data-types)

There is a long history of encapsulating expression data in S3 and S4 objects in R, going back to the `ExpressionSet` class in Biobase which was designed to store a matrix of microarray data alongside associated experimental metadata. This concept of storing all the relevant data and metadata in a single object has persisted through the development of RNA-Seq analysis (e.g.Â `SummarizedExperiment`) and into the age of single-cell transcriptomics (e.g.Â `SingleCellExperiment` \- see below).

### 1.4.1 SpatialExperiment class[](practical-session-1.html#spatialexperiment-class)

For the first part of this tutorial (practical sessions 1 and 2), we will be using the [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) S4 class from Bioconductor as the main data structure for storing and manipulating datasets.

`SpatialExperiment` is a specialized object class that supports the storing of spatially-resolved transcriptomics datasets within the Bioconductor framework. It builds on the [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) class (Amezquita et al. 2020) for single-cell RNA sequencing data, which itself extends the `RangedSummarizedExperiment` class. More specifically, `SpatialExperiment` has extra customisations to store spatial information (i.e., spatial coordinates and images).

An overview of the `SpatialExperiment` object structure is is presented in [1.1](practical-session-1.html#fig:SpExp-overview). In brief, the `SpatialExperiment` object consists of the below five parts:

  1. `assays`: gene expression counts  

  2. `rowData`: information about features, usually genes  

  3. `colData`: information on spots (non-spatial and spatial metadata)  

  4. `spatialCoords`: spatial coordinates  

  5. `imgData`: image data



**NOTE:** For spot-based STx data (i.e., 10x Genomics Visium), a single `assay` named `counts` is used.

![Overview of the `SpatialExperiment` object class structure.](images/SpatialExperiment.png)

Figure 1.1: Overview of the `SpatialExperiment` object class structure. 

For more details, see the related publication from Righelli et al., 2021 describing the `SpatialExperiment` (Righelli et al. 2022).

### 1.4.2 Inspect the object[](practical-session-1.html#inspect-the-object)
    
    
    [](practical-session-1.html#cb2-1)## Check the object's structure
    [](practical-session-1.html#cb2-2)spe
    
    
    ## class: SpatialExperiment 
    ## dim: 33538 4992 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
    ##   ENSG00000268674
    ## rowData names(3): gene_id gene_name feature_type
    ## colnames(4992): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
    ##   TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
    ## colData names(7): barcode_id sample_id ... ground_truth cell_count
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):
    ## spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
    ## imgData names(4): sample_id image_id data scaleFactor
    
    
    [](practical-session-1.html#cb4-1)## Check number of features/genes (rows) and spots (columns)
    [](practical-session-1.html#cb4-2)dim(spe)
    
    
    ## [1] 33538  4992
    
    
    [](practical-session-1.html#cb6-1)## Check names of 'assay' tables
    [](practical-session-1.html#cb6-2)assayNames(spe)
    
    
    ## [1] "counts"

### 1.4.3 Counts table and gene metadata[](practical-session-1.html#counts-table-and-gene-metadata)
    
    
    [](practical-session-1.html#cb8-1)## Have a look at the counts table
    [](practical-session-1.html#cb8-2)assay(spe)[1:6,1:4]
    
    
    ## 6 x 4 sparse Matrix of class "dgTMatrix"
    ##                 AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1
    ## ENSG00000243485                  .                  .                  .
    ## ENSG00000237613                  .                  .                  .
    ## ENSG00000186092                  .                  .                  .
    ## ENSG00000238009                  .                  .                  .
    ## ENSG00000239945                  .                  .                  .
    ## ENSG00000239906                  .                  .                  .
    ##                 AAACACCAATAACTGC-1
    ## ENSG00000243485                  .
    ## ENSG00000237613                  .
    ## ENSG00000186092                  .
    ## ENSG00000238009                  .
    ## ENSG00000239945                  .
    ## ENSG00000239906                  .

As we can see here the counts table is an object of class `dgTMatrix` which is a sparse matrix. This is because much like scRNA-seq data, STx data include many zeros. As a result, to make the counts table as light as possible we resort to using sparse matrices. This next code chunk examines a part of the matrix that includes genes with some level of expression:
    
    
    [](practical-session-1.html#cb10-1)assay(spe)[20:40, 2000:2010]
    
    
    ## 21 x 11 sparse Matrix of class "dgTMatrix"
    ##                                      
    ## ENSG00000223764 . . . . . . . . . . .
    ## ENSG00000187634 . . . . . . . . . . .
    ## ENSG00000188976 . . 2 . . . . . . 1 1
    ## ENSG00000187961 . . . . . . . . . . .
    ## ENSG00000187583 . . . . . . . . . . .
    ## ENSG00000187642 . . . . . . . . . . .
    ## ENSG00000272512 . . . . . . . . . . .
    ## ENSG00000188290 1 . . . . . . . . 2 .
    ## ENSG00000187608 . 1 . . . . 2 . . 1 .
    ## ENSG00000224969 . . . . . . . . . . .
    ## ENSG00000188157 . 1 . . 2 . . . . 1 .
    ## ENSG00000273443 . . . . . . . . . . .
    ## ENSG00000237330 . . . . . . . . . . .
    ## ENSG00000131591 . . . . . . . . . 1 .
    ## ENSG00000223823 . . . . . . . . . . .
    ## ENSG00000272141 . . . . . . . . . . .
    ## ENSG00000205231 . . . . . . . . . . .
    ## ENSG00000162571 . . . . . . . . . . .
    ## ENSG00000186891 . . . 1 . . . . . . .
    ## ENSG00000186827 . . . . . . . . . . .
    ## ENSG00000078808 . 1 2 . 1 . . . . 1 .
    
    
    [](practical-session-1.html#cb12-1)assay(spe)[33488:33508, 2000:2010]
    
    
    ## 21 x 11 sparse Matrix of class "dgTMatrix"
    ##                                                     
    ## ENSG00000160294  .  .   .  .   .   .  .  .  .   .  .
    ## ENSG00000228137  .  .   .  .   .   .  .  .  .   .  .
    ## ENSG00000239415  .  .   .  .   .   .  .  .  .   .  .
    ## ENSG00000182362  .  .   .  .   .   .  .  .  1   .  .
    ## ENSG00000160298  .  .   .  .   .   .  .  .  .   .  .
    ## ENSG00000160299  .  .   1  .   1   .  .  .  .   .  .
    ## ENSG00000160305  .  .   .  .   .   2  .  .  .   .  .
    ## ENSG00000160307  1  3   1  1   4   5  1  1  .   2  1
    ## ENSG00000160310  .  .   .  .   1   .  .  .  .   2  .
    ## ENSG00000198888 17 44  71 16 154  97 12 14 32 167  6
    ## ENSG00000198763 16 59  64 11 116  63 11 12 18 123  6
    ## ENSG00000198804 37 85 155 25 252 176 24 27 38 335 12
    ## ENSG00000198712 23 79 120 23 214 170 22 25 48 242 10
    ## ENSG00000228253  2  .   3  .   1   .  .  1  1   6  .
    ## ENSG00000198899 20 39  93  9 136 108 20 18 25 165  7
    ## ENSG00000198938 27 59 133 20 216 120 22 26 43 232  9
    ## ENSG00000198840  5 27  33  5  71  39  8 11 12  78  .
    ## ENSG00000212907  2  .   4  2   7   5  .  1  1   9  .
    ## ENSG00000198886 15 65  95  9 183  98 18 19 33 178  7
    ## ENSG00000198786  2 10  10  3  20  14  1  2  2  25  4
    ## ENSG00000198695  1  1   3  .   2   2  .  .  .   1  .

The levels of expression of different genes in the same spots differ significantly with many low values being present. We have to remember here that this data is not as yet normalized, and is therefore affected by systematic factors such as library size. Nonetheless, what is demonstrated here is typical for STx data (as it is for scRNA-seq data) - many genes will show low expression in individual spots.

To continue our exploration of the information stored in the `SpatialExperiment` object:
    
    
    [](practical-session-1.html#cb14-1)## Have a look at the genes metadata
    [](practical-session-1.html#cb14-2)head(rowData(spe))
    
    
    ## DataFrame with 6 rows and 3 columns
    ##                         gene_id   gene_name    feature_type
    ##                     <character> <character>     <character>
    ## ENSG00000243485 ENSG00000243485 MIR1302-2HG Gene Expression
    ## ENSG00000237613 ENSG00000237613     FAM138A Gene Expression
    ## ENSG00000186092 ENSG00000186092       OR4F5 Gene Expression
    ## ENSG00000238009 ENSG00000238009  AL627309.1 Gene Expression
    ## ENSG00000239945 ENSG00000239945  AL627309.3 Gene Expression
    ## ENSG00000239906 ENSG00000239906  AL627309.2 Gene Expression

### 1.4.4 Coordinates table and spot metadata[](practical-session-1.html#coordinates-table-and-spot-metadata)

The data that distinguished a `SpatialExperiment` object is the coordinate data which describes the spatial location of each spot.
    
    
    [](practical-session-1.html#cb16-1)## Check the spatial coordinates
    [](practical-session-1.html#cb16-2)head(spatialCoords(spe))
    
    
    ##                    pxl_col_in_fullres pxl_row_in_fullres
    ## AAACAACGAATAGTTC-1               3913               2435
    ## AAACAAGTATCTCCCA-1               9791               8468
    ## AAACAATCTACTAGCA-1               5769               2807
    ## AAACACCAATAACTGC-1               4068               9505
    ## AAACAGAGCGACTCCT-1               9271               4151
    ## AAACAGCTTTCAGAAG-1               3393               7583
    
    
    [](practical-session-1.html#cb18-1)## spot-level metadata
    [](practical-session-1.html#cb18-2)head(colData(spe))
    
    
    ## DataFrame with 6 rows and 7 columns
    ##                            barcode_id     sample_id in_tissue array_row
    ##                           <character>   <character> <integer> <integer>
    ## AAACAACGAATAGTTC-1 AAACAACGAATAGTTC-1 sample_151673         0         0
    ## AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1 sample_151673         1        50
    ## AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1 sample_151673         1         3
    ## AAACACCAATAACTGC-1 AAACACCAATAACTGC-1 sample_151673         1        59
    ## AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1 sample_151673         1        14
    ## AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1 sample_151673         1        43
    ##                    array_col ground_truth cell_count
    ##                    <integer>  <character>  <integer>
    ## AAACAACGAATAGTTC-1        16           NA         NA
    ## AAACAAGTATCTCCCA-1       102       Layer3          6
    ## AAACAATCTACTAGCA-1        43       Layer1         16
    ## AAACACCAATAACTGC-1        19           WM          5
    ## AAACAGAGCGACTCCT-1        94       Layer3          2
    ## AAACAGCTTTCAGAAG-1         9       Layer5          4

### 1.4.5 Image metadata[](practical-session-1.html#image-metadata)

Finally, the `SpatialExperiment` object also contains the image data from the STx experiment, giving the coordinates we looked at in the previous section some context in terms of the tissue of origin.
    
    
    [](practical-session-1.html#cb20-1)## Have a look at the image metadata
    [](practical-session-1.html#cb20-2)imgData(spe)
    
    
    ## DataFrame with 2 rows and 4 columns
    ##       sample_id    image_id   data scaleFactor
    ##     <character> <character> <list>   <numeric>
    ## 1 sample_151673      lowres   ####   0.0450045
    ## 2 sample_151673       hires   ####   0.1500150

As well as this (fairly basic) metadata, the `spe` object also contains the image itself, which the `SpatialExperiment` class allows us to access, like so:
    
    
    [](practical-session-1.html#cb22-1)## retrieve the image
    [](practical-session-1.html#cb22-2)spi <- getImg(spe)
    [](practical-session-1.html#cb22-3)## "plot" the image
    [](practical-session-1.html#cb22-4)plot(imgRaster(spi))

![](_main_files/figure-html/01_plot-image-1.png)

We can also use the scaling factors in the `imgData` to plot the locations of the Visium spots over the image. The position of a point in an image does not map directly to the spot location in cartesian coordinates, as it is the top-left of an image that is (0,0), not the bottom-left. In order to manage this, we need to transform the y-axis coordinates.
    
    
    [](practical-session-1.html#cb23-1)## "Plot" the image
    [](practical-session-1.html#cb23-2)plot(imgRaster(spi))
    [](practical-session-1.html#cb23-3)## Extract the spot locations
    [](practical-session-1.html#cb23-4)spot_coords <- spatialCoords(spe) %>% as.data.frame
    [](practical-session-1.html#cb23-5)## Scale by low-res factor
    [](practical-session-1.html#cb23-6)lowres_scale <- imgData(spe)[imgData(spe)$image_id == 'lowres', 'scaleFactor']
    [](practical-session-1.html#cb23-7)spot_coords$x_axis <- spot_coords$pxl_col_in_fullres * lowres_scale
    [](practical-session-1.html#cb23-8)spot_coords$y_axis <- spot_coords$pxl_row_in_fullres * lowres_scale
    [](practical-session-1.html#cb23-9)## lowres image is 600x600 pixels
    [](practical-session-1.html#cb23-10)dim(imgRaster(spi))
    
    
    ## [1] 600 600
    
    
    [](practical-session-1.html#cb25-1)## flip the Y axis
    [](practical-session-1.html#cb25-2)spot_coords$y_axis <- abs(spot_coords$y_axis - (ncol(imgRaster(spi)) + 1))
    [](practical-session-1.html#cb25-3)points(x=spot_coords$x_axis, y=spot_coords$y_axis)

![](_main_files/figure-html/01_plot-spots-1.png)

An equivalent plot, using `ggplot2` as the plotting library:
    
    
    [](practical-session-1.html#cb26-1)ggplot(mapping = aes(1:600, 1:600)) +
    [](practical-session-1.html#cb26-2)  annotation_raster(imgRaster(spi), xmin = 1, xmax = 600, ymin = 1, ymax = 600) +
    [](practical-session-1.html#cb26-3)  geom_point(data=spot_coords, aes(x=x_axis, y=y_axis), alpha=0.2) + xlim(1, 600) + ylim(1, 600) +
    [](practical-session-1.html#cb26-4)  coord_fixed() + 
    [](practical-session-1.html#cb26-5)  theme_void()

![](_main_files/figure-html/01_ggplot-spots-1.png) We can also extract additional metadata to make these plots more informative - for instance, the annotation from `colData` that flags whether a spot is âon tissueâ or not can be used to colour the spots like so:
    
    
    [](practical-session-1.html#cb27-1)## Add the annotation to the coordinate data frame
    [](practical-session-1.html#cb27-2)spot_coords$on_tissue <- as.logical(colData(spe)$in_tissue)
    [](practical-session-1.html#cb27-3)
    [](practical-session-1.html#cb27-4)ggplot(mapping = aes(1:600, 1:600)) +
    [](practical-session-1.html#cb27-5)  annotation_raster(imgRaster(spi), xmin = 1, xmax = 600, ymin = 1, ymax = 600) +
    [](practical-session-1.html#cb27-6)  geom_point(data=spot_coords, aes(x=x_axis, y=y_axis, colour=on_tissue), alpha=0.2) + xlim(1, 600) + ylim(1, 600) +
    [](practical-session-1.html#cb27-7)  coord_fixed() + 
    [](practical-session-1.html#cb27-8)  theme_void()

![](_main_files/figure-html/01_ggplot_ontissue-1.png)

Having to manually extract the relevant information from the `SpatialExperiment` object to generate plots like this does not generally make sense, and defies the point of using a data class that can encapsulate this information. We can instead use a package like `ggspavis`, which is explicitly built for generating visualisations of STx data directly from the `SpatialExperiment` object. We will make extensive use of this package during the next tutorial as we work through quality control processes for STx data. The pre-built nature of these plots is convenient, though it prevents users from achieving tasks like adding the tissue image to the plot. For many users the convenience will outweigh any issues this presents, though it is worth being aware of how to build visualisations from the ground up.
    
    
    [](practical-session-1.html#cb28-1)plotSpots(spe, in_tissue = NULL, annotate='in_tissue', size=0.5)

![](_main_files/figure-html/01_ggspavis-ontissue-1.png)

## 1.5 Conclusion[](practical-session-1.html#conclusion)

This first practical session has been a pretty straight-forward examination of an example Visium dataset. Weâve demonstrated where in this object the data and metadata are stored, how to extract it and make simple use of it.

### References[](references.html#references)

Amezquita, Robert A., Aaron T. L. Lun, Etienne Becht, Vince J. Carey, Lindsay N. Carpp, Ludwig Geistlinger, Federico Marini, et al. 2020. âOrchestrating single-cell analysis with Bioconductor.â _Nat Methods_ 17 (February): 137â45. <https://doi.org/10.1038/s41592-019-0654-x>. 

Maynard, Kristen R., Leonardo Collado-Torres, Lukas M. Weber, Cedric Uytingco, Brianna K. Barry, Stephen R. Williams, Joseph L. Catallini, et al. 2021. âTranscriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex.â _Nat Neurosci_ 24 (March): 425â36. <https://doi.org/10.1038/s41593-020-00787-0>. 

Righelli, Dario, Lukas M. Weber, Helena L. Crowell, Brenda Pardo, Leonardo Collado-Torres, Shila Ghazanfar, Aaron T. L. Lun, Stephanie C. Hicks, and Davide Risso. 2022. âSpatialExperiment: infrastructure for spatially-resolved transcriptomics data in R using Bioconductor.â _Bioinformatics_ 38 (11): 3128â31. <https://doi.org/10.1093/bioinformatics/btac299>. 


<!-- PAGE: practical-session-2.html -->

# Chapter 2 Practical session 2[](practical-session-2.html#practical-session-2)

Having previously introduced some of the Bioconductor ecosystem for storing and manipulating STx data, in this second session we will focus on some of the most common STx analysis tasks - particularly quality control assessment and associated spot- and gene- level filtering. We will also consider some global methods of STx analysis, including dimensionality reduction and clustering. All of the methods demonstrated here continue to focus on interoperable packages available via Bioconductor.
    
    
    [](practical-session-2.html#cb29-1)## Load packages {-}
    [](practical-session-2.html#cb29-2)library(SpatialExperiment)
    [](practical-session-2.html#cb29-3)library(STexampleData)
    [](practical-session-2.html#cb29-4)library(ggspavis)
    [](practical-session-2.html#cb29-5)library(ggplot2)
    [](practical-session-2.html#cb29-6)library(scater)
    [](practical-session-2.html#cb29-7)library(scran)
    [](practical-session-2.html#cb29-8)library(igraph)
    [](practical-session-2.html#cb29-9)library(pheatmap)
    [](practical-session-2.html#cb29-10)library(ggExtra)

  * [`ggspavis`](https://bioconductor.org/packages/release/bioc/html/ggspavis.html) is a Bioconductor package that includes visualization functions for spatially resolved transcriptomics datasets stored in `SpatialExperiment` format from spot-based (e.g., 10x Genomics Visium) platforms (Weber and Crowell (2022)).

  * [`scater`](https://bioconductor.org/packages/release/bioc/html/scater.html) is also a Bioconductor package that is a selection of tools for doing various analyses of scRNA-seq gene expression data, with a focus on quality control and visualization which has extended applications to STx data too. It is based on the `SingleCellExperiment` and `SpatialExperiment` classes and thus is interoperable with many other Bioconductor packages such as [`scran`](Spot-level%20quality%20control%20\(sQC\)%20procedures%20are%20employed%20to%20eliminate%20low-quality%20spots%20before%20conducting%20further%20analyses), [`scuttle`](https://bioconductor.org/packages/release/scuttle) and [`iSEE`](https://bioconductor.org/packages/release/iSEE).



    
    
    [](practical-session-2.html#cb30-1)## Reload the example dataset
    [](practical-session-2.html#cb30-2)spe <- Visium_humanDLPFC()
    
    
    ## see ?STexampleData and browseVignettes('STexampleData') for documentation
    
    
    ## loading from cache

## 2.1 Spot-level Quality Control[](practical-session-2.html#spot-level-quality-control)

Considered quality control (QC) procedures are essential for analysing any high-throughput data in molecular biology. The removal of noise and low quality data from complex datasets can improve the reliability of downsrtream analyses. STx is no different in this regard, and QC can be undertaken in 2 main places - spot-level and gene-level. Here, we focus on spot-level QC.

Spot-level quality control (sQC) procedures are employed to eliminate low-quality spots before conducting further analyses. Low-quality spots may result from issues during library preparation or other experimental procedures, such as a high percentage of dead cells due to cell damage during library preparation, or low mRNA capture efficiency caused by ineffective reverse transcription or PCR amplification. Keeping these spots usually leads to creating problems during downstream analyses.

We can identify low-quality spots using several characteristics that are also used in cell-level QC for scRNA-sq data, including:

  1. **library size** (total of UMI counts per spot will vary due to sequencing _-like different samples in a bulk RNA-seq-_ , or due to number of cells in the spot)
  2. **number of expressed genes** (i.e.Â number of genes with non-zero UMI counts per spot)
  3. **proportion of reads mapping to mitochondrial genes** (a high proportion indicates putative cell damage)



Low library size or low number of expressed features can indicate poor mRNA capture rates, e.g.Â due to cell damage and missing mRNAs, or low reaction efficiency. A high proportion of mitochondrial reads indicates cell damage, e.g.Â partial cell lysis leading to leakage and missing cytoplasmic mRNAs, with the resulting reads therefore concentrated on the remaining mitochondrial mRNAs that are relatively protected inside the mitochondrial membrane. Unusually high numbers of cells per spot can indicate problems during cell segmentation.

The idea of using scRNA-seq QC metrics in STx data comes from the fact that if we remove space and effectively treat each spot as a single cell, the two datasets share common features. We need to bear in mind, however, that the expected distributions for high-quality _spots_ are different (compared to high-quality _cells_ in scRNA-seq), since spots may contain zero, one, or multiple cells.

A few publications for further reading that can help you understand the quality controls: McCarthy et al. (2017) and Amezquita et al. (2020).

### 2.1.1 Plot tissue map[](practical-session-2.html#plot-tissue-map)

The dorso-lateral prefrontal cortex (DLPFC) is a functional brain region in primates involved in executive function. It consists of six layers of neurons that differ in their cell types, density and connections. The DLPFC dataset we looked at in session one, and will be here using comes with manual annotation of these layers (and the adjacent white matter - WM) by the authors Maynard et al. (2021). We can plot the tissue map with and without the annotations to get a complete view.
    
    
    [](practical-session-2.html#cb33-1)## Plot spatial coordinates without annotations
    [](practical-session-2.html#cb33-2)plotSpots(spe)
    [](practical-session-2.html#cb33-3)
    [](practical-session-2.html#cb33-4)## Plot spatial coordinates with annotations
    [](practical-session-2.html#cb33-5)plotSpots(spe,
    [](practical-session-2.html#cb33-6)          annotate = "ground_truth")

![](_main_files/figure-html/02_plot-maps-gTruth-1.png)![](_main_files/figure-html/02_plot-maps-gTruth-2.png)

### 2.1.2 Calculating QC metrics[](practical-session-2.html#calculating-qc-metrics)

We will calculate the three main QC metrics described above using methods from the `scater` (McCarthy et al. 2017) package, and investigate their influence on the DLPFC dataset with some plots from `ggspavis`, along with some additional plots of our own.

At present, the dataset contains both on- and off-tissue spots - we plotted these in the previous practical. For any future analysis though we are only interested in the on-tissue spots. Therefore, before we run any calculations we want to remove the off-tissue spots.

**_NOTE_** : the on- or off-tissue information for each spot can be found in the `colData` of the `spe` object and in the `in_tissue` column where _0 = off-tissue_ and _1 = on-tissue_.
    
    
    [](practical-session-2.html#cb34-1)## Dataset dimensions before the filtering
    [](practical-session-2.html#cb34-2)dim(spe)
    
    
    ## [1] 33538  4992
    
    
    [](practical-session-2.html#cb36-1)## Subset to keep only on-tissue spots
    [](practical-session-2.html#cb36-2)spe <- spe[, colData(spe)$in_tissue == 1]
    [](practical-session-2.html#cb36-3)dim(spe)
    
    
    ## [1] 33538  3639

The next thing we need to do before we make decisions on how to quality _âtrimâ_ the dataset is to calculate the percentage per spot of mitochodrial gene expression and store this information inside the `colData`. First of all, find the mitochrondrial genes - their gene names start with âMT-â or âmt-â.
    
    
    [](practical-session-2.html#cb38-1)## Classify genes as "mitochondrial" (is_mito == TRUE) 
    [](practical-session-2.html#cb38-2)## or not (is_mito == FALSE)
    [](practical-session-2.html#cb38-3)is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    [](practical-session-2.html#cb38-4)rowData(spe)$gene_name[is_mito]
    
    
    ##  [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3" 
    ##  [8] "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"

Then find what proportion of reads in a spotâs library are attributable to the expression of these genes. This uses a function, `addPerCellQC()` from `scater` (which in this instance is actually a wrapper around `scuttle`).
    
    
    [](practical-session-2.html#cb40-1)## Calculate per-spot QC metrics and store in colData
    [](practical-session-2.html#cb40-2)spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
    [](practical-session-2.html#cb40-3)head(colData(spe))
    
    
    ## DataFrame with 6 rows and 13 columns
    ##                            barcode_id     sample_id in_tissue array_row
    ##                           <character>   <character> <integer> <integer>
    ## AAACAAGTATCTCCCA-1 AAACAAGTATCTCCCA-1 sample_151673         1        50
    ## AAACAATCTACTAGCA-1 AAACAATCTACTAGCA-1 sample_151673         1         3
    ## AAACACCAATAACTGC-1 AAACACCAATAACTGC-1 sample_151673         1        59
    ## AAACAGAGCGACTCCT-1 AAACAGAGCGACTCCT-1 sample_151673         1        14
    ## AAACAGCTTTCAGAAG-1 AAACAGCTTTCAGAAG-1 sample_151673         1        43
    ## AAACAGGGTCTATATT-1 AAACAGGGTCTATATT-1 sample_151673         1        47
    ##                    array_col ground_truth cell_count       sum  detected
    ##                    <integer>  <character>  <integer> <numeric> <numeric>
    ## AAACAAGTATCTCCCA-1       102       Layer3          6      8458      3586
    ## AAACAATCTACTAGCA-1        43       Layer1         16      1667      1150
    ## AAACACCAATAACTGC-1        19           WM          5      3769      1960
    ## AAACAGAGCGACTCCT-1        94       Layer3          2      5433      2424
    ## AAACAGCTTTCAGAAG-1         9       Layer5          4      4278      2264
    ## AAACAGGGTCTATATT-1        13       Layer6          6      4004      2178
    ##                    subsets_mito_sum subsets_mito_detected subsets_mito_percent
    ##                           <numeric>             <numeric>            <numeric>
    ## AAACAAGTATCTCCCA-1             1407                    13              16.6351
    ## AAACAATCTACTAGCA-1              204                    11              12.2376
    ## AAACACCAATAACTGC-1              430                    13              11.4089
    ## AAACAGAGCGACTCCT-1             1316                    13              24.2223
    ## AAACAGCTTTCAGAAG-1              651                    12              15.2174
    ## AAACAGGGTCTATATT-1              621                    13              15.5095
    ##                        total
    ##                    <numeric>
    ## AAACAAGTATCTCCCA-1      8458
    ## AAACAATCTACTAGCA-1      1667
    ## AAACACCAATAACTGC-1      3769
    ## AAACAGAGCGACTCCT-1      5433
    ## AAACAGCTTTCAGAAG-1      4278
    ## AAACAGGGTCTATATT-1      4004

After calculating a required metric, we need to apply a cut-off threshold for the metric to decide whether or not to keep each spot. It is important to consider an individual dataset on its own merits, as it might need slightly different cut-off values to be applied. As a result we cannot rely on identifying a single value to use every time and we need to rely on plotting these metrics and making a decision on a dataset-by-dataset basis.

### 2.1.3 Library size threshold plot[](practical-session-2.html#library-size-threshold-plot)

We can plot a histogram of the library sizes across spots. The library size is the number of UMI counts in each spot. We can find this information in the `sum` column in the `colData`.
    
    
    [](practical-session-2.html#cb42-1)## Density and histogram of library sizes
    [](practical-session-2.html#cb42-2)ggplot(data = as.data.frame(colData(spe)),
    [](practical-session-2.html#cb42-3)       aes(x = sum)) +
    [](practical-session-2.html#cb42-4)  geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-2.html#cb42-5)                 colour = "black", 
    [](practical-session-2.html#cb42-6)                 fill = "grey") +
    [](practical-session-2.html#cb42-7)  geom_density(alpha = 0.5,
    [](practical-session-2.html#cb42-8)               adjust = 1.0,
    [](practical-session-2.html#cb42-9)               fill = "#A0CBE8",
    [](practical-session-2.html#cb42-10)               colour = "#4E79A7") +
    [](practical-session-2.html#cb42-11)  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb42-12)  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb42-13)  xlab("Library size") + 
    [](practical-session-2.html#cb42-14)  ylab("Density") + 
    [](practical-session-2.html#cb42-15)  theme_classic()

![](_main_files/figure-html/02_plot-libSize-histo-1.png)

As we can see there are no obvious issues with the library sizes. An example of an issue could be a high frequency of small libraries which would indicate poor experimental output. Generally we do not want to keep spots with too small libraries.

If the dataset we are analysing contains the number of cells that are present in each spot (this one does), then it makes sense to also plot the library sizes against the number of cells per spot. In that way we are making sure that we donât remove any spots that may have biological meaning. In many cases though the datasets do not have such information unless we can generate it using a nuclei segmentation tool to extract this information from the H&E images.

The horizontal red line (argument `threshold` in the `plotQC` function) shows a first guess at a possible filtering threshold for library size based on the above histogram.
    
    
    [](practical-session-2.html#cb43-1)## Scatter plot, library size against number of cells per spot
    [](practical-session-2.html#cb43-2)plotQC(spe, type = "scatter", 
    [](practical-session-2.html#cb43-3)       metric_x = "cell_count", metric_y = "sum", 
    [](practical-session-2.html#cb43-4)       threshold_y = 700)

![](_main_files/figure-html/02_plot-libSizeVScelNo-1.png)

**NOTE** : The `ggspavis` plots for QC are convenient, but not very configurable. As can be seen from the âmissingâ bin in the top histogram here, the default configuration provided is not always the best. A ggplot2 alternative (using `ggExtra` to provide the marginal histograms) is also provided here.
    
    
    [](practical-session-2.html#cb44-1)p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=sum)) +
    [](practical-session-2.html#cb44-2)  geom_point(size=0.5) + 
    [](practical-session-2.html#cb44-3)  geom_smooth(se=FALSE) +
    [](practical-session-2.html#cb44-4)  geom_hline(yintercept = 700, colour='red') + 
    [](practical-session-2.html#cb44-5)  theme_minimal()
    [](practical-session-2.html#cb44-6)ggMarginal(p, type='histogram', margins = 'both')

![](_main_files/figure-html/02_ggplot-libSizeVScelNo-1.png)

We need to keep in mind here that the threshold is, to an extent, arbitrary. It is therefore important to look at the number of spots that are left out of the dataset by this choice of cut-off value, and also have a look at their putative spatial patterns. If we filtered out spots with biological relevance, then we should observe some patterns on the tissue map that correlate with some of the known biological structures of the tissue. If we do observe such a phenomenon, we have probably set our threshold too high (i.e.Â not permissive enough).
    
    
    [](practical-session-2.html#cb45-1)## Select library size threshold
    [](practical-session-2.html#cb45-2)qc_lib_size <- colData(spe)$sum < 700
    [](practical-session-2.html#cb45-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb45-4)table(qc_lib_size)
    
    
    ## qc_lib_size
    ## FALSE  TRUE 
    ##  3628    11
    
    
    [](practical-session-2.html#cb47-1)## Add threshold in colData
    [](practical-session-2.html#cb47-2)colData(spe)$qc_lib_size <- qc_lib_size
    [](practical-session-2.html#cb47-3)
    [](practical-session-2.html#cb47-4)## Check putative spatial patterns of removed spots
    [](practical-session-2.html#cb47-5)plotQC(spe, type = "spots", 
    [](practical-session-2.html#cb47-6)       discard = "qc_lib_size")

![](_main_files/figure-html/02_libSize-thresh-1.png)

As an optional exercise, try to illustrate what happens if we set the threshold too high (i.e., 2000 UMI counts).

**NOTE:** For reference, remember the ground truth layers in this dataset [that we plotted](practical-session-2.html#plot-tissue-map) at the beginning of this session.
    
    
    [](practical-session-2.html#cb48-1)## Select library size threshold
    [](practical-session-2.html#cb48-2)code...
    [](practical-session-2.html#cb48-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb48-4)code...
    [](practical-session-2.html#cb48-5)## Add threshold in colData
    [](practical-session-2.html#cb48-6)code...
    [](practical-session-2.html#cb48-7)
    [](practical-session-2.html#cb48-8)## Check putative spatial patterns of removed spots
    [](practical-session-2.html#cb48-9)plotQC(...)

### 2.1.4 Number of expressed genes[](practical-session-2.html#number-of-expressed-genes)

As we did with the library sizes, we can plot a histogram of the number of expressed genes across spots. A gene is âexpressedâ in a spot if it has at least one count in it. We can find this information in the `detected` column in the `colData`.

We will follow the same logic for the plots as we did for the library size earlier.
    
    
    [](practical-session-2.html#cb49-1)## Density and histogram of expressed genes
    [](practical-session-2.html#cb49-2)ggplot(data = as.data.frame(colData(spe)),
    [](practical-session-2.html#cb49-3)       aes(x = detected)) +
    [](practical-session-2.html#cb49-4)  geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-2.html#cb49-5)                 colour = "black", 
    [](practical-session-2.html#cb49-6)                 fill = "grey") +
    [](practical-session-2.html#cb49-7)  geom_density(alpha = 0.5,
    [](practical-session-2.html#cb49-8)               adjust = 1.0,
    [](practical-session-2.html#cb49-9)               fill = "#A0CBE8",
    [](practical-session-2.html#cb49-10)               colour = "#4E79A7") +
    [](practical-session-2.html#cb49-11)  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb49-12)  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb49-13)  xlab("Genes expressed in each spot") + 
    [](practical-session-2.html#cb49-14)  ylab("Density") + 
    [](practical-session-2.html#cb49-15)  theme_classic()

![](_main_files/figure-html/02_plot-genesInSpot-histo-1.png)
    
    
    [](practical-session-2.html#cb50-1)# plot number of expressed genes vs. number of cells per spot
    [](practical-session-2.html#cb50-2)p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=detected)) +
    [](practical-session-2.html#cb50-3)  geom_point(size=0.5) + 
    [](practical-session-2.html#cb50-4)  geom_smooth(se=FALSE) +
    [](practical-session-2.html#cb50-5)  geom_hline(yintercept = 500, colour='red') + 
    [](practical-session-2.html#cb50-6)  theme_minimal()
    [](practical-session-2.html#cb50-7)ggMarginal(p, type='histogram', margins = 'both')

![](_main_files/figure-html/02_genesInSpot-scatter-1.png)

Finally, again as before, we apply the chosen threshold to flag spots with (in this case) fewer than 500 expressed genes.
    
    
    [](practical-session-2.html#cb51-1)## Select expressed genes threshold
    [](practical-session-2.html#cb51-2)qc_detected <- colData(spe)$detected < 500
    [](practical-session-2.html#cb51-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb51-4)table(qc_detected)
    
    
    ## qc_detected
    ## FALSE  TRUE 
    ##  3628    11
    
    
    [](practical-session-2.html#cb53-1)## Add threshold in colData
    [](practical-session-2.html#cb53-2)colData(spe)$qc_detected <- qc_detected
    [](practical-session-2.html#cb53-3)
    [](practical-session-2.html#cb53-4)## Check for putative spatial pattern of removed spots
    [](practical-session-2.html#cb53-5)plotQC(spe, type = "spots", 
    [](practical-session-2.html#cb53-6)       discard = "qc_detected")

![](_main_files/figure-html/02_genesInSpot-thresh-1.png)

Again, an optional exercise is provided to see the effects of an over-enthusiastic filter - to illustrate what happens if we set the threshold too high (i.e., 1000 expressed genes).

**NOTE:** For reference, remember the ground truth layers in this dataset [that we plotted](practical-session-2.html#plot-tissue-map) at the beginning of this session.
    
    
    [](practical-session-2.html#cb54-1)## Select library size threshold
    [](practical-session-2.html#cb54-2)code...
    [](practical-session-2.html#cb54-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb54-4)code...
    [](practical-session-2.html#cb54-5)## Add threshold in colData
    [](practical-session-2.html#cb54-6)code...
    [](practical-session-2.html#cb54-7)
    [](practical-session-2.html#cb54-8)## Check putative spatial patterns of removed spots
    [](practical-session-2.html#cb54-9)plotQC(...)

### 2.1.5 Percentage of mitochondrial expression[](practical-session-2.html#percentage-of-mitochondrial-expression)

As we briefly touched on at the beginning, a high proportion of mitochondrial reads indicates low cell quality, probably due to cell damage.

We calculated this data earlier on in this session, and can now investigate the percentage of mitochondrial expression across spots by looking at the column `subsets_mito_percent` in the `colData`.
    
    
    [](practical-session-2.html#cb55-1)## Density and histogram of percentage of mitochondrial expression
    [](practical-session-2.html#cb55-2)ggplot(data = as.data.frame(colData(spe)),
    [](practical-session-2.html#cb55-3)       aes(x = subsets_mito_percent)) +
    [](practical-session-2.html#cb55-4)  geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-2.html#cb55-5)                 colour = "black", 
    [](practical-session-2.html#cb55-6)                 fill = "grey") +
    [](practical-session-2.html#cb55-7)  geom_density(alpha = 0.5,
    [](practical-session-2.html#cb55-8)               adjust = 1.0,
    [](practical-session-2.html#cb55-9)               fill = "#A0CBE8",
    [](practical-session-2.html#cb55-10)               colour = "#4E79A7") +
    [](practical-session-2.html#cb55-11)  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb55-12)  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb55-13)  xlab("Percentage of mitochondrial expression") + 
    [](practical-session-2.html#cb55-14)  ylab("Density") + 
    [](practical-session-2.html#cb55-15)  theme_classic()

![](_main_files/figure-html/02_plot-mitoPercent-histo-1.png)

In this instance, a higher percentage of mitochondrial expression is the thing to avoid, so the threshold is an upper bound, rather than the lower bounds we have observed so far. Our suggestion this time is to cut-off at 28%.
    
    
    [](practical-session-2.html#cb56-1)# plot mitochondrial read proportion vs. number of cells per spot
    [](practical-session-2.html#cb56-2)p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=subsets_mito_percent)) +
    [](practical-session-2.html#cb56-3)  geom_point(size=0.5) + 
    [](practical-session-2.html#cb56-4)  geom_smooth(se=FALSE) +
    [](practical-session-2.html#cb56-5)  geom_hline(yintercept = 28, colour='red') + 
    [](practical-session-2.html#cb56-6)  theme_minimal()
    [](practical-session-2.html#cb56-7)ggMarginal(p, type='histogram')

![](_main_files/figure-html/02_mitoPercent-scatter-1.png)
    
    
    [](practical-session-2.html#cb57-1)## Select expressed genes threshold
    [](practical-session-2.html#cb57-2)qc_mito <- colData(spe)$subsets_mito_percent > 28
    [](practical-session-2.html#cb57-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb57-4)table(qc_mito)
    
    
    ## qc_mito
    ## FALSE  TRUE 
    ##  3622    17
    
    
    [](practical-session-2.html#cb59-1)## Add threshold in colData
    [](practical-session-2.html#cb59-2)colData(spe)$qc_mito <- qc_mito
    [](practical-session-2.html#cb59-3)
    [](practical-session-2.html#cb59-4)## Check for putative spatial pattern of removed spots
    [](practical-session-2.html#cb59-5)plotQC(spe, type = "spots", 
    [](practical-session-2.html#cb59-6)       discard = "qc_mito")

![](_main_files/figure-html/02_mitoPercent-thresh-1.png)

Again, try to illustrate what happens if we set the threshold too low (i.e., 20 0r 25%).

**NOTE:** For reference, remember the ground truth layers in this dataset [that we plotted](practical-session-2.html#plot-tissue-map) at the beginning of this session.
    
    
    [](practical-session-2.html#cb60-1)## Select library size threshold
    [](practical-session-2.html#cb60-2)code...
    [](practical-session-2.html#cb60-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb60-4)code...
    [](practical-session-2.html#cb60-5)## Add threshold in colData
    [](practical-session-2.html#cb60-6)code...
    [](practical-session-2.html#cb60-7)
    [](practical-session-2.html#cb60-8)## Check putative spatial patterns of removed spots
    [](practical-session-2.html#cb60-9)plotQC(...)

### 2.1.6 Number of cells per spot[](practical-session-2.html#number-of-cells-per-spot)

As previously mentioned, number of cells per spot is an attribute that not all datasets include. Nonetheless, it can be useful to further control the quality of the dataset prior to any downstream analysis. Of course, the number of cells per spot depends on the tissue type and organism and according to [10X Genomics](https://kb.10xgenomics.com/hc/en-us/articles/360035487952-How-many-cells-are-captured-in-a-single-spot-), each spot typically contains between 0 and 10 cells.

The DPFLC dataset does contain information on the number of cells per spot (acquired by processing and cell segmentation of high-resolution histology images obtained prior on-slide cDNA synthesis, see Maynard et al. (2021) for details). To investigate the number of cells in each spot looking for any outlier values that could indicate problems we need to take a look in the column `cell_count` in `colData`.
    
    
    [](practical-session-2.html#cb61-1)## Density and histogram of the number of cells in each spot
    [](practical-session-2.html#cb61-2)ggplot(data = as.data.frame(colData(spe)),
    [](practical-session-2.html#cb61-3)       aes(x = cell_count)) +
    [](practical-session-2.html#cb61-4)  geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-2.html#cb61-5)                 binwidth = 1,
    [](practical-session-2.html#cb61-6)                 colour = "black", 
    [](practical-session-2.html#cb61-7)                 fill = "grey") +
    [](practical-session-2.html#cb61-8)  geom_density(alpha = 0.5,
    [](practical-session-2.html#cb61-9)               adjust = 1.5,
    [](practical-session-2.html#cb61-10)               fill = "#A0CBE8",
    [](practical-session-2.html#cb61-11)               colour = "#4E79A7") +
    [](practical-session-2.html#cb61-12)  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb61-13)  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb61-14)  xlab("Number of cells per spot") + 
    [](practical-session-2.html#cb61-15)  ylab("Density") + 
    [](practical-session-2.html#cb61-16)  theme_classic()

![](_main_files/figure-html/02_plot-cellsPerSpot-histo-1.png)
    
    
    [](practical-session-2.html#cb62-1)## Have a look at the values
    [](practical-session-2.html#cb62-2)table(colData(spe)$cell_count)
    
    
    ## 
    ##   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
    ##  84 211 483 623 617 541 421 287 140  92  50  25  18  10   9   3   8   2   1   2 
    ##  20  21  22  23  25  26  27 
    ##   3   2   1   1   2   2   1
    
    
    [](practical-session-2.html#cb64-1)# plot number of expressed genes vs. number of cells per spot
    [](practical-session-2.html#cb64-2)p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=detected)) +
    [](practical-session-2.html#cb64-3)  geom_point() + 
    [](practical-session-2.html#cb64-4)  geom_smooth(se=FALSE) +
    [](practical-session-2.html#cb64-5)  geom_vline(xintercept = 10, colour='red') + 
    [](practical-session-2.html#cb64-6)  theme_minimal()
    [](practical-session-2.html#cb64-7)ggMarginal(p, type='histogram')

![](_main_files/figure-html/02_cellsPerSpot-scatter-1.png)

As we can see from both the histogram and the scatter plot there is a tail of very high values, which could indicate problems for these spots. More specifically, we can see from the scatter plot that most of the spots with very high cell counts also tend to have lower numbers of expressed genes. This indicates problems with the experiment on these spots, and they should be removed.
    
    
    [](practical-session-2.html#cb65-1)## Select expressed genes threshold
    [](practical-session-2.html#cb65-2)qc_cell_count <- colData(spe)$cell_count > 10
    [](practical-session-2.html#cb65-3)## Check how many spots are filtered out
    [](practical-session-2.html#cb65-4)table(qc_cell_count)
    
    
    ## qc_cell_count
    ## FALSE  TRUE 
    ##  3549    90
    
    
    [](practical-session-2.html#cb67-1)## Add threshold in colData
    [](practical-session-2.html#cb67-2)colData(spe)$qc_cell_count <- qc_cell_count
    [](practical-session-2.html#cb67-3)
    [](practical-session-2.html#cb67-4)## Check for putative spatial pattern of removed spots
    [](practical-session-2.html#cb67-5)plotQC(spe, type = "spots", 
    [](practical-session-2.html#cb67-6)       discard = "qc_cell_count")

![](_main_files/figure-html/02_cellsPerSpot-thresh-1.png)

While there is a spatial pattern to the discarded spots, it does not appear to be correlated with the known biological features (cortical layers). The discarded spots are typically at the edges of the tissue. It seems plausible that something has gone wrong with the cell segmentation on the edges of the images, so it makes sense to remove these spots.

### 2.1.7 Remove low-quality spots[](practical-session-2.html#remove-low-quality-spots)

All the steps so far have flagged spots with potential issues - before proceeding with analysis, we want to remove these spots from our SpatialExperiment object. Since we have calculated different spot-level QC metrics and selected thresholds for each one, we can combine them to identify a set of low-quality spots, and remove them from our `spe` object in a single step.

We can also check once more that the combined set of discarded spots does not correspond to any obvious biologically relevant group of spots.
    
    
    [](practical-session-2.html#cb68-1)## Check the number of discarded spots for each metric
    [](practical-session-2.html#cb68-2)apply(cbind(qc_lib_size, qc_detected, qc_mito, qc_cell_count), 2, sum)
    
    
    ##   qc_lib_size   qc_detected       qc_mito qc_cell_count 
    ##            11            11            17            90
    
    
    [](practical-session-2.html#cb70-1)## Combine together the set of discarded spots
    [](practical-session-2.html#cb70-2)discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
    [](practical-session-2.html#cb70-3)## Store the set in the object
    [](practical-session-2.html#cb70-4)colData(spe)$discard <- discard
    [](practical-session-2.html#cb70-5)
    [](practical-session-2.html#cb70-6)## Check the spatial pattern of combined set of discarded spots
    [](practical-session-2.html#cb70-7)plotQC(spe, type = "spots", 
    [](practical-session-2.html#cb70-8)       discard = "discard")

![](_main_files/figure-html/02_checkQC-thresh-1.png)

Since this dataset has also manual annotation ([remember](practical-session-2.html#plot-tissue-map))) we see that there are locations that are not annotated (marked with `NA`). We could further remove those locations to reduce potential noise and further increase the quality of the dataset.
    
    
    [](practical-session-2.html#cb71-1)## Select locations without annotation
    [](practical-session-2.html#cb71-2)qc_NA_spots <- is.na(colData(spe)$ground_truth)
    [](practical-session-2.html#cb71-3)## Combine together the set of discarded spots
    [](practical-session-2.html#cb71-4)discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count | qc_NA_spots
    [](practical-session-2.html#cb71-5)## Store the set in the object
    [](practical-session-2.html#cb71-6)colData(spe)$discard <- discard
    [](practical-session-2.html#cb71-7)
    [](practical-session-2.html#cb71-8)## Check the spatial pattern of combined set of discarded spots
    [](practical-session-2.html#cb71-9)plotQC(spe, type = "spots", 
    [](practical-session-2.html#cb71-10)       discard = "discard")

![](_main_files/figure-html/02_notAnnotSpots-1.png)
    
    
    [](practical-session-2.html#cb72-1)## remove combined set of low-quality spots
    [](practical-session-2.html#cb72-2)spe <- spe[, !colData(spe)$discard]

## 2.2 Normalisation of counts[](practical-session-2.html#normalisation-of-counts)

### 2.2.1 Background[](practical-session-2.html#background)

Normalisation is applied in STx data for the same reason as any other RNA-Seq technique - the differences observed in the count data can arise from a range of systematic factors, not just a physiologically-relevant change in expression. The primary systematic effect is that of library size (or in the case of STx, counts/UMIs per spot). `scater` corrects for library size by scaling the sizes across all spots such that the mean library size is 1. Normalized counts are then calculated as a ratio of observed count to library size factor.

Secondly, a log-transformation is applied to the scaled counts - this transformation is commonly applied as it stabilises the variance across the range of transcriptomics data (otherwise the variance is dominated by highly expressed genes) and it facilitates comparisons of expression by rendering positive and negative changes symmetrical and found by subtraction rather than division. Since \\(log2(0)\\) is undefined, a _pseudocount_ is added to each observed count to avoid this error - a pseudocount of 1 is typically applied, as \\(log2(0+1) = 0\\).

Here we will be using methods from the `scater` (McCarthy et al. 2017) and `scran` (Lun, McCarthy, and Marioni 2016) packages that calculate logcounts using library size factors. The library size factors approach is arguably the simplest approach for STx data. Other approaches used in scRNA-seq are more difficult to justify their use in STx because of two main reasons:

  1. Spots can contain multiple cells of different cell-types.
  2. Datasets can include multiple tissue samples which will lead to different clusterings.



### 2.2.2 Log-tranformation of counts[](practical-session-2.html#log-tranformation-of-counts)
    
    
    [](practical-session-2.html#cb73-1)## Calculate library size factors
    [](practical-session-2.html#cb73-2)spe <- computeLibraryFactors(spe)
    [](practical-session-2.html#cb73-3)## Have a look at the size factors
    [](practical-session-2.html#cb73-4)summary(sizeFactors(spe))
    
    
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1514  0.6326  0.9011  1.0000  1.2849  3.7500

As described above, the mean size factor is 1.0.
    
    
    [](practical-session-2.html#cb75-1)## Density and histogram of library sizes
    [](practical-session-2.html#cb75-2)ggplot(data = data.frame(sFact = sizeFactors(spe)), 
    [](practical-session-2.html#cb75-3)       aes(x = sFact)) +
    [](practical-session-2.html#cb75-4)  geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-2.html#cb75-5)                 colour = "black", 
    [](practical-session-2.html#cb75-6)                 fill = "grey") +
    [](practical-session-2.html#cb75-7)  geom_density(alpha = 0.5,
    [](practical-session-2.html#cb75-8)               adjust = 1.0,
    [](practical-session-2.html#cb75-9)               fill = "#A0CBE8",
    [](practical-session-2.html#cb75-10)               colour = "#4E79A7") +
    [](practical-session-2.html#cb75-11)  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb75-12)  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-2.html#cb75-13)  xlab("Library size") + 
    [](practical-session-2.html#cb75-14)  ylab("Density") + 
    [](practical-session-2.html#cb75-15)  theme_classic()

![](_main_files/figure-html/02_plot-labfact-histo-1.png)

The log-transformation that takes place is a log2-transformation and in order to avoid _\- Infinity_ values we add a pseudo value of 1. Both the log2- transformation and the pseudocount of 1 are defaults in this method.
    
    
    [](practical-session-2.html#cb76-1)## Calculate logcounts and store in the spe object
    [](practical-session-2.html#cb76-2)spe <- logNormCounts(spe)
    [](practical-session-2.html#cb76-3)
    [](practical-session-2.html#cb76-4)## Check that a new assay has been added
    [](practical-session-2.html#cb76-5)assayNames(spe)
    
    
    ## [1] "counts"    "logcounts"

## 2.3 Selecting genes[](practical-session-2.html#selecting-genes)

### 2.3.1 Background[](practical-session-2.html#background-1)

Gene selection - or alternatively âfeature selectionâ - is applied to identify genes that are likely to be informative for downstream analyses. The most common feature selection method is the definition of highly variable genes (HVGs). The assumption is that since we quality-controlled and normalised our dataset, the genes with high variability are the ones that contain high levels of biological variability too. Since here we have a spatial dataset we can also try to identify spatially variable genes too (SVGs).

It is important to note that HVGs are identified solely from the gene expression data. Spatial information does not play a role in finding HVGs. STx data pose a dilemma; does the meaningful spatial information reflect only spatial distribution of major cell types or does it reflect additional important spatial features? If we believe the former, relying on HVGs can be enough. If the second also holds true though, it is important to identify SVGs as well.

### 2.3.2 Highly Variable Genes (HVGs)[](practical-session-2.html#highly-variable-genes-hvgs)

Here we will be using methods from the `scran` package (Lun, McCarthy, and Marioni 2016) to identify a set of HVGs. Again, here we need to remember that `scran` methods were developed for scRNA-seq and we are performing the analysis under the assumption that the spots of an STx experiment can be treated as single cells.

In this dataset, the mitochondrial genes are too highly expressed and are not of major biological interest. As a result, if we are to identify true HVGs, we first need to remove the mitochondrial genes.
    
    
    [](practical-session-2.html#cb78-1)## Remove mitochondrial genes
    [](practical-session-2.html#cb78-2)spe <- spe[!is_mito, ]

Then, we apply methods from `scran` that give a list of HVGs, which can be used for further downstream analyses.

First we model the variance of the log-expression profiles for each gene, decomposing it into technical and biological components based on a fitted mean-variance trend.
    
    
    [](practical-session-2.html#cb79-1)## Fit mean-variance relationship
    [](practical-session-2.html#cb79-2)dec <- modelGeneVar(spe)
    [](practical-session-2.html#cb79-3)## Visualize mean-variance relationship
    [](practical-session-2.html#cb79-4)fit <- metadata(dec)
    [](practical-session-2.html#cb79-5)fit_df <- data.frame(mean = fit$mean,
    [](practical-session-2.html#cb79-6)                     var = fit$var,
    [](practical-session-2.html#cb79-7)                     trend = fit$trend(fit$mean))
    [](practical-session-2.html#cb79-8)
    [](practical-session-2.html#cb79-9)ggplot(data = fit_df, 
    [](practical-session-2.html#cb79-10)       aes(x = mean, y = var)) + 
    [](practical-session-2.html#cb79-11)  geom_point() + 
    [](practical-session-2.html#cb79-12)  geom_line(aes(y = trend), colour = "dodgerblue", linewidth = 1.5) + 
    [](practical-session-2.html#cb79-13)  labs(x = "mean of log-expression",
    [](practical-session-2.html#cb79-14)       y = "variance of log-expression") + 
    [](practical-session-2.html#cb79-15)  theme_classic()

![](_main_files/figure-html/02_features_FitModel-1.png)

The `trend` function that we used above is returned from the `modelGeneVar` function and returns the fitted value of the trend at any value of the mean. The âbiologicalâ variance of a gene is what remains when the fitted variance for a gene of that expression value is subtracted from the total variance (so genes above the blue trend line have a positive biological variance).

We select the top 10% of genes based on their biological variability The parameter `prop` defines how many HVGs we want. For example `prop = 0.1` returns the top 10% of genes. `prop = 1.0` would return all genes with a positive biological variability.
    
    
    [](practical-session-2.html#cb80-1)## Select top HVGs
    [](practical-session-2.html#cb80-2)top_hvgs <- getTopHVGs(dec, prop = 0.1)
    [](practical-session-2.html#cb80-3)
    [](practical-session-2.html#cb80-4)## How many HVGs?
    [](practical-session-2.html#cb80-5)length(top_hvgs)
    
    
    ## [1] 1429

**NOTE** \- we will return to feature selection in the next practical, as it is a complicated process with significant impacts on the chosen downstream analysis.

### 2.3.3 Spatially variable genes (SVGs)[](practical-session-2.html#spatially-variable-genes-svgs)

SVGs are genes with a highly spatially correlated pattern of expression, which varies along with the spatial distribution of a tissue structure of interest. This phenomenon is also called _spatial autocorrelation_ and underlies all types of spatial data, as we will discuss later.

The field of geography has developed some statistical measures to calculate spatial autocorrelation. Examples of these are Moranâs _I_ (âNotes on Continuous Stochastic Phenomena on JSTORâ 1950) and Gearyâs _C_ (âThe Contiguity Ratio and Statistical Mapping on JSTORâ 1954) that can be used to rank genes by the observed spatial autocorrelation to identify SVGs.

Several sophisticated new statistical methods to identify SVGs in STx data have also recently been developed. These include [SpatialDE](https://github.com/Teichlab/SpatialDE) (Svensson, Teichmann, and Stegle 2018), [SPARK](https://xzhoulab.github.io/SPARK/) (Sun, Zhu, and Zhou 2020), and [SPARK-X](https://xzhoulab.github.io/SPARK/) (Zhu, Sun, and Zhou 2021).

### 2.3.4 Integration of HVGs and SVGs[](practical-session-2.html#integration-of-hvgs-and-svgs)

A recent benchmark paper (Li et al. 2022) showed that integrating HVGs and SVGs to generate a combined set of features can improve downstream clustering performance in STx data. This confirms that SVGs contain additional biologically relevant information that is not captured by HVGs in these datasets. For example, a simple way to combine these features is to concatenate columns of principal components (PCs) calculated on the set of HVGs and the set of SVGs (excluding overlapping HVGs), and then using the combined set of features for further downstream analyses (Li et al. 2022).

## 2.4 Dimensionality reduction[](practical-session-2.html#dimensionality-reduction)

### 2.4.1 Background[](practical-session-2.html#background-2)

STx data, just like bulk and single-cell transcriptomics, is captured in high-dimensional space. The reduction of this complexity can be helpful for a number of applications. Principal Components Analysis (PCA) assumes linearity in the data and has historically been used for dimensionality reduction. More modern techniques, such as Uniform Manifold Approximation and Projection (UMAP) (McInnes, Healy, and Melville 2018) and t-Stochastic Neighbor Embedding (tSNE) (Maaten and Hinton 2008) do not assume linearity and provide some performance advantages.

The main practical difference between the output of these techniques is that the distances between the data points in PCA space are interpretable and can be used for clustering, while the distances in a UMAP/tSNE embedding are not interpretable in this way. As a result, we will be using PCA to reduce the dimensions of our dataset to assist clustering and UMAP to further reduce the principal components (PCs) in a two-dimensional space and produce better visualisations for the PCA.

Dimensionality reduction prior to clustering has two main advantages, firstly it reduces dataset noise from random variation. Secondly it improves the computational efficiency of downstream analyses such as clustering. In an STx experiment, like the one we are analysing here, we have more than 3,000 spots and almost 1,500 HVGs. As as result, each spot has 1,500 attributes (dimensions) as a basis for subsequent clustering. This large number of variables that differentiate or cluster together spots gives rise to the _curse of dimensionality_ (Keogh and Mueen 2017). This principle states that data points (spots) with a large number of features appear equidistant in attribute space resulting in poor clustering output.

### 2.4.2 PCA: Principal component analysis[](practical-session-2.html#pca-principal-component-analysis)

Here we will use an efficient implementation of PCA provided in the `scater` package (McCarthy et al. 2017) and retain the top 50 PCs for further downstream analyses. The random seed is required for reproducibility reasons because this implementation uses randomisation.
    
    
    [](practical-session-2.html#cb82-1)## Set seed
    [](practical-session-2.html#cb82-2)set.seed(987)
    [](practical-session-2.html#cb82-3)## Compute PCA
    [](practical-session-2.html#cb82-4)spe <- runPCA(spe, subset_row = top_hvgs)
    [](practical-session-2.html#cb82-5)## Check correctness - names
    [](practical-session-2.html#cb82-6)reducedDimNames(spe)
    
    
    ## [1] "PCA"
    
    
    [](practical-session-2.html#cb84-1)## Check correctness - dimensions
    [](practical-session-2.html#cb84-2)dim(reducedDim(spe, "PCA"))
    
    
    ## [1] 3511   50

### 2.4.3 UMAP: Uniform Manifold Approximation and Projection[](practical-session-2.html#umap-uniform-manifold-approximation-and-projection)

Here we will also run UMAP - using `scater`âs implementation - on the 50 PCs generated above and retain the top 2 UMAP components to visualise results.
    
    
    [](practical-session-2.html#cb86-1)## Set seed
    [](practical-session-2.html#cb86-2)set.seed(987)
    [](practical-session-2.html#cb86-3)## Compute UMAP on top 50 PCs
    [](practical-session-2.html#cb86-4)spe <- runUMAP(spe, dimred = "PCA")
    [](practical-session-2.html#cb86-5)## Check correctness - names
    [](practical-session-2.html#cb86-6)reducedDimNames(spe)
    
    
    ## [1] "PCA"  "UMAP"
    
    
    [](practical-session-2.html#cb88-1)## Check correctness - dimensions
    [](practical-session-2.html#cb88-2)dim(reducedDim(spe, "UMAP"))
    
    
    ## [1] 3511    2
    
    
    [](practical-session-2.html#cb90-1)## Update column names for easier plotting
    [](practical-session-2.html#cb90-2)colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

### 2.4.4 UMAP visualisations[](practical-session-2.html#umap-visualisations)

We can generate plots either using plotting functions from the [ggspavis](https://bioconductor.org/packages/ggspavis) package or [`ggplot2`](https://ggplot2.tidyverse.org/) package. When clustering later, we will add cluster labels to these reduced dimension plots for an off-tissue visualisation.
    
    
    [](practical-session-2.html#cb91-1)## Plot top 2 PCA dimensions
    [](practical-session-2.html#cb91-2)# plotDimRed(spe, type = "PCA")
    [](practical-session-2.html#cb91-3)
    [](practical-session-2.html#cb91-4)ggplot(data = as.data.frame(spe@int_colData@listData$reducedDims$PCA),
    [](practical-session-2.html#cb91-5)       aes(x = PC1, y = PC2, colour = spe@colData$ground_truth)) + 
    [](practical-session-2.html#cb91-6)  geom_point(size = 0.5) + 
    [](practical-session-2.html#cb91-7)  scale_colour_brewer(type = "qual") + 
    [](practical-session-2.html#cb91-8)  labs(title = "Reduced dimensions: PCA",
    [](practical-session-2.html#cb91-9)       x = "PC1",
    [](practical-session-2.html#cb91-10)       y = "PC2",
    [](practical-session-2.html#cb91-11)       colour = "Layers") +
    [](practical-session-2.html#cb91-12)  theme_classic()
    [](practical-session-2.html#cb91-13)
    [](practical-session-2.html#cb91-14)## Plot top 2 UMAP dimensions
    [](practical-session-2.html#cb91-15)# plotDimRed(spe, type = "UMAP")
    [](practical-session-2.html#cb91-16)
    [](practical-session-2.html#cb91-17)ggplot(data = as.data.frame(spe@int_colData@listData$reducedDims$UMAP),
    [](practical-session-2.html#cb91-18)       aes(x = UMAP1, y = UMAP2, colour = spe@colData$ground_truth)) + 
    [](practical-session-2.html#cb91-19)  geom_point(size = 0.5) + 
    [](practical-session-2.html#cb91-20)  scale_colour_brewer(type = "qual") + 
    [](practical-session-2.html#cb91-21)  labs(title = "Reduced dimensions: UMAP",
    [](practical-session-2.html#cb91-22)       x = "UMAP1",
    [](practical-session-2.html#cb91-23)       y = "UMAP2",
    [](practical-session-2.html#cb91-24)       colour = "Layers") +
    [](practical-session-2.html#cb91-25)  theme_classic()

![](_main_files/figure-html/02_dimRed_UMAP-vis-1.png)![](_main_files/figure-html/02_dimRed_UMAP-vis-2.png)

## 2.5 Clustering[](practical-session-2.html#clustering)

### 2.5.1 Background[](practical-session-2.html#background-3)

The clustering of observations into statistically similar groups is a well-established application in both bulk and single-cell RNA-Seq analysis. Clustering is a helpful tool because it structures and orders the data, allowing useful insights to be gained from complex, multivariate datasets and use those insights to classify the observed data or to generate hypotheses.

Common clustering methods are applied to ST data based on correlation or statistical distance of gene expression measurements. As we briefly mentioned above, the dimensionality of STx data means that sample distances in gene expression space tend to be small and not reliable for identifying clusters, so feature selection (gene selection) or dimensionality reduction approaches (i.e., PCA, UMAP) tend to be taken before clustering.

Common approaches to clustering gene expression data include k-means, hierarchical and Louvain algorithms, and all have been applied to the clustering of ST data. Some of these methods are implemented in some of the most popular single-cell analysis packages, such as `Seurat` (Hao et al. 2021) and `scran` (Lun, McCarthy, and Marioni 2016) and have been used for clustering in a number of ST studies.

### 2.5.2 Clustering on HVGs[](practical-session-2.html#clustering-on-hvgs)

Here, we apply graph-based clustering to the top 50 PCs calculated on the set of selected HVGs, using the Walktrap method implemented in `scran` (Lun, McCarthy, and Marioni 2016). To do so, we assume that (i) each spot is equal to a cell and (ii) we can detect from the gene expression the biologically informative spatial distribution patterns of cell types.
    
    
    [](practical-session-2.html#cb92-1)## Set seed
    [](practical-session-2.html#cb92-2)set.seed(987)
    [](practical-session-2.html#cb92-3)## Set number of Nearest-Neighbours (NNs)
    [](practical-session-2.html#cb92-4)k <- 10
    [](practical-session-2.html#cb92-5)## Build the k-NN graph
    [](practical-session-2.html#cb92-6)g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
    [](practical-session-2.html#cb92-7)## Run walktrap clustering
    [](practical-session-2.html#cb92-8)g_walk <- igraph::cluster_walktrap(g)
    [](practical-session-2.html#cb92-9)## Get the cluster labels
    [](practical-session-2.html#cb92-10)clus <- g_walk$membership
    [](practical-session-2.html#cb92-11)## Check how many
    [](practical-session-2.html#cb92-12)table(clus)
    
    
    ## clus
    ##   1   2   3   4   5   6 
    ## 350 354 661 895 366 885
    
    
    [](practical-session-2.html#cb94-1)## Store cluster labels in column 'label' in colData
    [](practical-session-2.html#cb94-2)colLabels(spe) <- factor(clus)

### 2.5.3 HVGs clustering visualisations[](practical-session-2.html#hvgs-clustering-visualisations)

We can visualise the clusters in two ways:

  1. plotting in spatial coordinates on the tissue map
  2. plotting in the UMAP/PCA embeddings.



We can use plotting functions either from the [ggspavis](https://bioconductor.org/packages/ggspavis) package.

For reference, we will also display the ground truth (manually annotated) labels available for this dataset.
    
    
    [](practical-session-2.html#cb95-1)## Plot in tissue map
    [](practical-session-2.html#cb95-2)plotSpots(spe, annotate = "label", 
    [](practical-session-2.html#cb95-3)          palette = "libd_layer_colors")
    [](practical-session-2.html#cb95-4)
    [](practical-session-2.html#cb95-5)## Plot ground truth in tissue map
    [](practical-session-2.html#cb95-6)plotSpots(spe, annotate = "ground_truth", 
    [](practical-session-2.html#cb95-7)          palette = "libd_layer_colors")

![](_main_files/figure-html/02_clust_vis-map-1.png)![](_main_files/figure-html/02_clust_vis-map-2.png)
    
    
    [](practical-session-2.html#cb96-1)## Plot clusters in PCA space
    [](practical-session-2.html#cb96-2)plotDimRed(spe, type = "PCA", 
    [](practical-session-2.html#cb96-3)           annotate = "label", palette = "libd_layer_colors")
    [](practical-session-2.html#cb96-4)
    [](practical-session-2.html#cb96-5)## Plot clusters in UMAP space
    [](practical-session-2.html#cb96-6)plotDimRed(spe, type = "UMAP", 
    [](practical-session-2.html#cb96-7)           annotate = "label", palette = "libd_layer_colors")

![](_main_files/figure-html/02_clust_vis-DimRed-1.png)![](_main_files/figure-html/02_clust_vis-DimRed-2.png)

From the visualizations, we can see that the clustering reproduces, to an extent, the known biological structure of the tissue, but not perfectly. One reason for this could be the fact that each spot may be comprised of a number different cells whose gene expression profiles are diluted in the overall profile of the spot, thus leading to low-quality clustering.

### 2.5.4 Spatially-aware clustering[](practical-session-2.html#spatially-aware-clustering)

In STx data, we can also perform clustering that takes spatial information into account, for example to identify spatially compact or spatially connected clusters.

A simple strategy is to perform graph-based clustering on a set of features (columns) that includes both molecular features (gene expression) and spatial features (x-y coordinates). In this case, a crucial tuning parameter is the relative amount of scaling between the two data modalities â if the scaling is chosen poorly, either the molecular or spatial features will dominate the clustering. Depending on data availability, further modalities could also be included. In this section, we will include some examples on this clustering approach.

## 2.6 Inter-cluster differentially expressed genes (DGEs)[](practical-session-2.html#inter-cluster-differentially-expressed-genes-dges)

### 2.6.1 Background[](practical-session-2.html#background-4)

Here, we will identify differentially expressed genes between clusters.

We will use the `findMarkers` implementation from the `scran` (Lun, McCarthy, and Marioni 2016). This implementation uses a binomial test, which tests for genes that differ in the proportion expressed vs.Â not expressed between clusters. This is a more stringent test than the default _t_ -tests, and tends to select genes that are easier to interpret and validate experimentally.

### 2.6.2 DGEs identification[](practical-session-2.html#dges-identification)
    
    
    [](practical-session-2.html#cb97-1)## Set gene names as row names ease of plotting
    [](practical-session-2.html#cb97-2)rownames(spe) <- rowData(spe)$gene_name
    [](practical-session-2.html#cb97-3)## Test for DGEs
    [](practical-session-2.html#cb97-4)markers <- findMarkers(spe, test = "binom", direction = "up")
    [](practical-session-2.html#cb97-5)## Check output
    [](practical-session-2.html#cb97-6)markers
    
    
    ## List of length 6
    ## names(6): 1 2 3 4 5 6

The output from the `findMarkers` implementation is a list of length equal to the number of clusters. Each element of the list contains the Log-Fold-Change (LogFC) of each gene between one cluster and all others.

### 2.6.3 DGEs visualisation[](practical-session-2.html#dges-visualisation)

Here we will plot LogFCs for cluster 1 against all other clusters
    
    
    [](practical-session-2.html#cb99-1)## Select cluster 1 genes
    [](practical-session-2.html#cb99-2)interesting <- markers[[1]]
    [](practical-session-2.html#cb99-3)## Get the top genes
    [](practical-session-2.html#cb99-4)best_set <- interesting[interesting$Top <= 5, ]
    [](practical-session-2.html#cb99-5)## Calculate the effect
    [](practical-session-2.html#cb99-6)logFCs <- getMarkerEffects(best_set)
    [](practical-session-2.html#cb99-7)## Plot a heat map
    [](practical-session-2.html#cb99-8)pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))

![](_main_files/figure-html/02_dges_vis-clst1-1.png)

Below we will plot the log-transformed normalised expression of the top genes for one cluster alongside their expression in the other clusters.
    
    
    [](practical-session-2.html#cb100-1)## Select genes
    [](practical-session-2.html#cb100-2)top_genes <- head(rownames(interesting))
    [](practical-session-2.html#cb100-3)## Plot expression
    [](practical-session-2.html#cb100-4)plotExpression(spe, x = "label", features = top_genes)

![](_main_files/figure-html/02_dges_vis-2-1.png)

## 2.7 Putting it all together[](practical-session-2.html#putting-it-all-together)
    
    
    [](practical-session-2.html#cb101-1)# clear workspace from previous chapters
    [](practical-session-2.html#cb101-2)rm(list = ls(all = TRUE))
    [](practical-session-2.html#cb101-3)
    [](practical-session-2.html#cb101-4)# LOAD DATA
    [](practical-session-2.html#cb101-5)
    [](practical-session-2.html#cb101-6)library(SpatialExperiment)
    [](practical-session-2.html#cb101-7)library(STexampleData)
    [](practical-session-2.html#cb101-8)spe <- Visium_humanDLPFC()
    [](practical-session-2.html#cb101-9)
    [](practical-session-2.html#cb101-10)# QUALITY CONTROL (QC)
    [](practical-session-2.html#cb101-11)
    [](practical-session-2.html#cb101-12)library(scater)
    [](practical-session-2.html#cb101-13)# subset to keep only spots over tissue
    [](practical-session-2.html#cb101-14)spe <- spe[, colData(spe)$in_tissue == 1]
    [](practical-session-2.html#cb101-15)# identify mitochondrial genes
    [](practical-session-2.html#cb101-16)is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    [](practical-session-2.html#cb101-17)# calculate per-spot QC metrics
    [](practical-session-2.html#cb101-18)spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
    [](practical-session-2.html#cb101-19)# select QC thresholds
    [](practical-session-2.html#cb101-20)qc_lib_size <- colData(spe)$sum < 600
    [](practical-session-2.html#cb101-21)qc_detected <- colData(spe)$detected < 400
    [](practical-session-2.html#cb101-22)qc_mito <- colData(spe)$subsets_mito_percent > 28
    [](practical-session-2.html#cb101-23)qc_cell_count <- colData(spe)$cell_count > 10
    [](practical-session-2.html#cb101-24)# combined set of discarded spots
    [](practical-session-2.html#cb101-25)discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
    [](practical-session-2.html#cb101-26)colData(spe)$discard <- discard
    [](practical-session-2.html#cb101-27)# filter low-quality spots
    [](practical-session-2.html#cb101-28)spe <- spe[, !colData(spe)$discard]
    [](practical-session-2.html#cb101-29)
    [](practical-session-2.html#cb101-30)# NORMALIZATION
    [](practical-session-2.html#cb101-31)
    [](practical-session-2.html#cb101-32)library(scran)
    [](practical-session-2.html#cb101-33)# calculate logcounts using library size factors
    [](practical-session-2.html#cb101-34)spe <- logNormCounts(spe)
    [](practical-session-2.html#cb101-35)
    [](practical-session-2.html#cb101-36)# FEATURE SELECTION
    [](practical-session-2.html#cb101-37)
    [](practical-session-2.html#cb101-38)# remove mitochondrial genes
    [](practical-session-2.html#cb101-39)spe <- spe[!is_mito, ]
    [](practical-session-2.html#cb101-40)# fit mean-variance relationship
    [](practical-session-2.html#cb101-41)dec <- modelGeneVar(spe)
    [](practical-session-2.html#cb101-42)# select top HVGs
    [](practical-session-2.html#cb101-43)top_hvgs <- getTopHVGs(dec, prop = 0.1)
    [](practical-session-2.html#cb101-44)
    [](practical-session-2.html#cb101-45)# DIMENSIONALITY REDUCTION
    [](practical-session-2.html#cb101-46)
    [](practical-session-2.html#cb101-47)# compute PCA
    [](practical-session-2.html#cb101-48)set.seed(123)
    [](practical-session-2.html#cb101-49)spe <- runPCA(spe, subset_row = top_hvgs)
    [](practical-session-2.html#cb101-50)# compute UMAP on top 50 PCs
    [](practical-session-2.html#cb101-51)set.seed(123)
    [](practical-session-2.html#cb101-52)spe <- runUMAP(spe, dimred = "PCA")
    [](practical-session-2.html#cb101-53)# update column names
    [](practical-session-2.html#cb101-54)colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
    [](practical-session-2.html#cb101-55)
    [](practical-session-2.html#cb101-56)# CLUSTERING
    [](practical-session-2.html#cb101-57)
    [](practical-session-2.html#cb101-58)# graph-based clustering
    [](practical-session-2.html#cb101-59)set.seed(123)
    [](practical-session-2.html#cb101-60)k <- 10
    [](practical-session-2.html#cb101-61)g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
    [](practical-session-2.html#cb101-62)g_walk <- igraph::cluster_walktrap(g)
    [](practical-session-2.html#cb101-63)clus <- g_walk$membership
    [](practical-session-2.html#cb101-64)colLabels(spe) <- factor(clus)
    [](practical-session-2.html#cb101-65)
    [](practical-session-2.html#cb101-66)# MARKER GENES
    [](practical-session-2.html#cb101-67)# test for marker genes
    [](practical-session-2.html#cb101-68)rownames(spe) <- rowData(spe)$gene_name
    [](practical-session-2.html#cb101-69)markers <- findMarkers(spe, test = "binom", direction = "up")

### References[](references.html#references)

Amezquita, Robert A., Aaron T. L. Lun, Etienne Becht, Vince J. Carey, Lindsay N. Carpp, Ludwig Geistlinger, Federico Marini, et al. 2020. âOrchestrating single-cell analysis with Bioconductor.â _Nat Methods_ 17 (February): 137â45. <https://doi.org/10.1038/s41592-019-0654-x>. 

Hao, Yuhan, Stephanie Hao, Erica Andersen-Nissen, William M. Mauck, Shiwei Zheng, Andrew Butler, Maddie J. Lee, et al. 2021. âIntegrated analysis of multimodal single-cell data.â _Cell_ 184 (13): 3573â3587.e29. <https://doi.org/10.1016/j.cell.2021.04.048>. 

Keogh, Eamonn, and Abdullah Mueen. 2017. âCurse of Dimensionality.â In _Encyclopedia of Machine Learning and Data Mining_ , 314â15. Boston, MA, USA: Springer, Boston, MA. <https://doi.org/10.1007/978-1-4899-7687-1_192>. 

Li, Yijun, Stefan Stanojevic, Bing He, Zheng Jing, Qianhui Huang, Jian Kang, and Lana X. Garmire. 2022. âBenchmarking Computational Integration Methods for Spatial Transcriptomics Data.â _bioRxiv_ , January, 2021.08.27.457741. <https://doi.org/10.1101/2021.08.27.457741>. 

Lun, Aaron T. L., Davis J. McCarthy, and John C. Marioni. 2016. âA step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.â _F1000Research_ 5 (2122): 2122. <https://doi.org/10.12688/f1000research.9501.2>. 

Maaten, Laurens van der, and Geoffrey Hinton. 2008. âVisualizing Data Using t-SNE.â _Journal of Machine Learning Research_ 9 (86): 2579â2605. <http://jmlr.org/papers/v9/vandermaaten08a.html>. 

Maynard, Kristen R., Leonardo Collado-Torres, Lukas M. Weber, Cedric Uytingco, Brianna K. Barry, Stephen R. Williams, Joseph L. Catallini, et al. 2021. âTranscriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex.â _Nat Neurosci_ 24 (March): 425â36. <https://doi.org/10.1038/s41593-020-00787-0>. 

McCarthy, Davis J., Kieran R. Campbell, Aaron T. L. Lun, and Quin F. Wills. 2017. âScater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R.â _Bioinformatics_ 33 (8): 1179â86. <https://doi.org/10.1093/bioinformatics/btw777>. 

McInnes, Leland, John Healy, and James Melville. 2018. âUMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.â _arXiv_ , February. <https://doi.org/10.48550/arXiv.1802.03426>. 

âNotes on Continuous Stochastic Phenomena on JSTOR.â 1950\. _Biometrika_. <https://www.jstor.org/stable/2332142>. 

Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 2020. âStatistical analysis of spatial expression patterns for spatially resolved transcriptomic studies.â _Nat Methods_ 17 (February): 193â200. <https://doi.org/10.1038/s41592-019-0701-7>. 

Svensson, Valentine, Sarah A. Teichmann, and Oliver Stegle. 2018. âSpatialDE: identification of spatially variable genes.â _Nat Methods_ 15 (May): 343â46. <https://doi.org/10.1038/nmeth.4636>. 

âThe Contiguity Ratio and Statistical Mapping on JSTOR.â 1954\. _Incorporated Statistician_. <https://www.jstor.org/stable/2986645>. 

Weber, Lukas M., and Helena L. Crowell. 2022. _Ggspavis: Visualization Functions for Spatially Resolved Transcriptomics Data_. <https://github.com/lmweber/ggspavis>. 

Zhu, Jiaqiang, Shiquan Sun, and Xiang Zhou. 2021. âSPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies.â _Genome Biol_ 22 (1): 1â25. <https://doi.org/10.1186/s13059-021-02404-0>. 


<!-- PAGE: practical-session-3.html -->

# Chapter 3 Practical session 3[](practical-session-3.html#practical-session-3)

This practical session will demonstrate the application of the most commonly used spatial analysis tools to STx data, and how we work with coordinate data alongside expression data.

## 3.1 Load packages[](practical-session-3.html#load-packages)

  * [`spdep`](https://cran.r-project.org/web/packages/spdep/index.html) is a collection of functions to create spatial weights matrix objects from polygon _contiguities_ , from point patterns by distance and tessellations. It is used for summarizing these objects, and for permitting their use in spatial data analysis like regional aggregation and tests for spatial _autocorrelation_.

  * [`sf`](https://cran.r-project.org/web/packages/sf/index.html) (_Simple Features for R_) is a package that offers support for simple features, a standardized way to encode spatial vector data.

  * [`GWmodel`](https://cran.r-project.org/web/packages/GWmodel/index.html) is a suite of models that fit situations when data are not described well by some global model, but where there are spatial regions where a suitably localised calibration provides a better description.




## 3.2 Background[](practical-session-3.html#background-5)

### 3.2.1 Main geocomputational data structures[](practical-session-3.html#main-geocomputational-data-structures)

There are three main data structures that we need to have ready before we undertake a geocomputational approach to STx data analysis. Namely these are; (1) geometries (point and polygon), (2) neighbours lists and (3) distance matrices.

  1. Spatial geometries can be points, lines, polygons and pixels. Polygons consist of a multitude of points connected by lines and can have many forms like circle, hexagon, non-canonical polygon etc.

  2. Neighbour lists are special types of lists that contain information about the neighbours of each polygon. The neighbours can be defined either by adjacency or by distance.

  3. Distance matrices contain the distances between different points and can be either weighted or un-weighted. The weighted distances are usually objective to each point and its neighbours. Meaning that the closer or farther a neighbour is from the point of focus, the weight of their distance changes according to an applied kernel. Usually in the case of STx data, like the ones generated by the 10X Visium platform, the un-weighted distance between two points is expressed in pixels and we acquire it from the `spaceranger` output.




### 3.2.2 The `sf` objects[](practical-session-3.html#the-sf-objects)

Package `sf` represents simple features as native R objects. All functions and methods in `sf` that operate on spatial data are prefixed by _st__ , which refers to _spatial type_. Simple features are implemented as R native data, using simple data structures (S3 classes, lists, matrix, vector). The typical use of `sf` involves reading, manipulating and writing of sets of features, with attributes and geometries.

As attributes are typically stored in `data.frame` objects (or the very similar `tbl_df`), we will also store feature geometries in a `data.frame` column. Since geometries are not single-valued, they are put in a list-column, a list of length equal to the number of records in the `data.frame`, with each list element holding the simple feature geometry of that feature. The three classes used to represent simple features are:

  * `sf`, the table (`data.frame`) with feature attributes and feature geometries, which contains
  * `sfc`, the list-column with the geometries for each feature (record), which is composed of
  * `sfg`, the feature geometry of an individual simple feature.



#### 3.2.2.1 Simple feature geometry types[](practical-session-3.html#simple-feature-geometry-types)

The following seven simple feature types are the most common:

type | description  
---|---  
`POINT` | zero-dimensional geometry containing a single point  
`LINESTRING` | sequence of points connected by straight, non-self intersecting line pieces; one-dimensional geometry  
`POLYGON` | geometry with a positive area (two-dimensional); sequence of points form a closed, non-self intersecting ring; the first ring denotes the exterior ring, zero or more subsequent rings denote holes in this exterior ring  
`MULTIPOINT` | set of points; a MULTIPOINT is simple if no two Points in the MULTIPOINT are equal  
`MULTILINESTRING` | set of linestrings  
`MULTIPOLYGON` | set of polygons  
`GEOMETRYCOLLECTION` | set of geometries of any type except GEOMETRYCOLLECTION  
  
Each of the geometry types can also be a (typed) empty set, containing zero coordinates (for `POINT` the standard is not clear how to represent the empty geometry). Empty geometries can be thought of as being the analogue to missing (`NA`) attributes, NULL values or empty lists.

#### 3.2.2.2 sf: objects with simple features[](practical-session-3.html#sf-objects-with-simple-features)

As we usually do not work with geometries of single `simple features`, but with datasets consisting of sets of features with attributes, the two are put together in `sf` (simple feature) objects. The following command reads a test dataset called `nc` from a file that is contained in the `sf` package:
    
    
    [](practical-session-3.html#cb102-1)nc <- st_read(system.file("shape/nc.shp", package = "sf"))
    
    
    ## Reading layer `nc' from data source 
    ##   `/home/sjcockell/R/x86_64-pc-linux-gnu-library/4.3/sf/shape/nc.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 100 features and 14 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
    ## Geodetic CRS:  NAD27

The short report printed gives the file name, the driver (ESRI Shapefile), mentions that there are 100 features (records, represented as rows) and 14 fields (attributes, represented as columns).

This object is of class:
    
    
    [](practical-session-3.html#cb104-1)class(nc)
    
    
    ## [1] "sf"         "data.frame"

meaning it extends (and âisâ a) `data.frame`, but with a single list-column with geometries, which is held in the column with name:
    
    
    [](practical-session-3.html#cb106-1)attr(nc, "sf_column")
    
    
    ## [1] "geometry"

If we print the first three features, we see their attribute values and an abridged version of the geometry
    
    
    [](practical-session-3.html#cb108-1)print(nc[9:15], n = 3)

which would give the following output:

![Overview of the `sf` object.](images/sf_xfig.png)

Figure 3.1: Overview of the `sf` object. 

In the output we see:

  * in green a simple feature: a single record, or `data.frame` row, consisting of attributes and geometry
  * in blue a single simple feature geometry (an object of class `sfg`)
  * in red a simple feature list-column (an object of class `sfc`, which is a column in the `data.frame`)
  * that although geometries are native R objects, they are printed as well-known text



It is also possible to create `data.frame` objects with geometry list-columns that are not of class `sf`, e.g.Â by:
    
    
    [](practical-session-3.html#cb109-1)nc.no_sf <- as.data.frame(nc)
    [](practical-session-3.html#cb109-2)class(nc.no_sf)
    
    
    ## [1] "data.frame"

However, such objects:

  * no longer register which column is the geometry list-column
  * no longer have a plot method, and
  * lack all of the other dedicated methods for class `sf`



#### 3.2.2.3 sfc: simple feature geometry list-column[](practical-session-3.html#sfc-simple-feature-geometry-list-column)

The column in the `sf` data.frame that contains the geometries is a list, of class `sfc`. We can retrieve the geometry list-column in this case by using standard `data.frame` notation like `nc$geom` or `nc[[15]]`, but the more general way uses `st_geometry`:
    
    
    [](practical-session-3.html#cb111-1)(nc_geom <- st_geometry(nc))
    
    
    ## Geometry set for 100 features 
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -84.32385 ymin: 33.88199 xmax: -75.45698 ymax: 36.58965
    ## Geodetic CRS:  NAD27
    ## First 5 geometries:
    
    
    ## MULTIPOLYGON (((-81.47276 36.23436, -81.54084 3...
    
    
    ## MULTIPOLYGON (((-81.23989 36.36536, -81.24069 3...
    
    
    ## MULTIPOLYGON (((-80.45634 36.24256, -80.47639 3...
    
    
    ## MULTIPOLYGON (((-76.00897 36.3196, -76.01735 36...
    
    
    ## MULTIPOLYGON (((-77.21767 36.24098, -77.23461 3...

Geometries are printed in abbreviated form, but we can view a complete geometry by selecting it, e.g.Â the first one by:
    
    
    [](practical-session-3.html#cb118-1)nc_geom[[1]]
    
    
    ## MULTIPOLYGON (((-81.47276 36.23436, -81.54084 36.27251, -81.56198 36.27359, -81.63306 36.34069, -81.74107 36.39178, -81.69828 36.47178, -81.7028 36.51934, -81.67 36.58965, -81.3453 36.57286, -81.34754 36.53791, -81.32478 36.51368, -81.31332 36.4807, -81.26624 36.43721, -81.26284 36.40504, -81.24069 36.37942, -81.23989 36.36536, -81.26424 36.35241, -81.32899 36.3635, -81.36137 36.35316, -81.36569 36.33905, -81.35413 36.29972, -81.36745 36.2787, -81.40639 36.28505, -81.41233 36.26729, -81.43104 36.26072, -81.45289 36.23959, -81.47276 36.23436)))

The way this is printed is called _well-known text_ , and is part of the standards. The word `MULTIPOLYGON` is followed by three parentheses, because it can consist of multiple polygons, in the form of `MULTIPOLYGON(POL1,POL2)`, where `POL1` might consist of an exterior ring and zero or more interior rings, as of `(EXT1,HOLE1,HOLE2)`. Sets of coordinates are held together with parentheses, so we get `((crds_ext)(crds_hole1)(crds_hole2))` where `crds_` is a comma-separated set of coordinates of a ring. This leads to the case above, where `MULTIPOLYGON(((crds_ext)))` refers to the exterior ring (1), without holes (2), of the first polygon (3) - hence three parentheses.

We can see there is a single polygon with no rings:
    
    
    [](practical-session-3.html#cb120-1)par(mar = c(0,0,1,0))
    [](practical-session-3.html#cb120-2)plot(nc[1], reset = FALSE) # reset = FALSE: we want to add to a plot with a legend
    [](practical-session-3.html#cb120-3)plot(nc[1,1], col = 'grey', add = TRUE)

![](_main_files/figure-html/03_sfc_Test3-1.png)

Following the `MULTIPOLYGON` data structure, in R we have a list of lists of lists of matrices. For instance, we get the first 3 coordinate pairs of the second exterior ring (first ring is always exterior) for the geometry of feature 4 by:
    
    
    [](practical-session-3.html#cb121-1)nc_geom[[4]][[2]][[1]][1:3,]
    
    
    ##           [,1]     [,2]
    ## [1,] -76.02717 36.55672
    ## [2,] -75.99866 36.55665
    ## [3,] -75.91192 36.54253

Geometry columns have their own class,
    
    
    [](practical-session-3.html#cb123-1)class(nc_geom)
    
    
    ## [1] "sfc_MULTIPOLYGON" "sfc"

#### 3.2.2.4 sfg: simple feature geometry[](practical-session-3.html#sfg-simple-feature-geometry)

Simple feature geometry (`sfg`) objects carry the geometry for a single feature, e.g.Â a point, linestring or polygon.

Simple feature geometries are implemented as R native data, using the following rules

  1. a single POINT is a numeric vector
  2. a set of points, e.g.Â in a LINESTRING or ring of a POLYGON is a `matrix`, each row containing a point
  3. any other set is a `list`



The below figure illustrates the different types of geometries:

![](_main_files/figure-html/03_sf_Test7-1.png)

Geometries can also be empty, as in
    
    
    [](practical-session-3.html#cb125-1)(x <- st_geometrycollection())
    [](practical-session-3.html#cb125-2)## GEOMETRYCOLLECTION EMPTY
    [](practical-session-3.html#cb125-3)length(x)
    [](practical-session-3.html#cb125-4)## [1] 0

_The above are taken from the very well written, well-descriptive and thorough`sf` package [vignette](https://cran.r-project.org/web/packages/sf/vignettes/sf1.html)._

## 3.3 Data structures preparation[](practical-session-3.html#data-structures-preparation)

For this practical we will be using a human steatotic kidney dataset from the [Liver Atlas](https://livercellatlas.org/index.php) (Guilliams et al. 2022). Specifically we will use the JBO019 sample.

### 3.3.1 Load new dataset[](practical-session-3.html#load-new-dataset)

**Note** \- Between now and section 3.8 nothing new is introduced (this is a repetition of the QC carried out in practical 2, but with this new liver dataset). Feel free to skip over this section and use the code block just before section 3.8 to run this QC quickly for this data.

First we generate the `SpatialFeaturesExperiment` object which is an extension of the `SpatialExperiment` (SPE) object that we used in the 2nd practical session. The difference is that the SFE object has incorporated the `sf` object structure and thus can accommodate the use of geocomputational tools.
    
    
    [](practical-session-3.html#cb126-1)sampleDir <- "./data/spaceranger_outs/Human_Liver_Steatotic/JBO019_Results"
    [](practical-session-3.html#cb126-2)sampleNames <- "JBO019"
    [](practical-session-3.html#cb126-3)sfe <- read10xVisiumSFE(samples = sampleDir, 
    [](practical-session-3.html#cb126-4)                        sample_id = sampleNames, 
    [](practical-session-3.html#cb126-5)                        type = "sparse", 
    [](practical-session-3.html#cb126-6)                        data = "filtered", 
    [](practical-session-3.html#cb126-7)                        images = "lowres", 
    [](practical-session-3.html#cb126-8)                        style = "W", 
    [](practical-session-3.html#cb126-9)                        zero.policy = TRUE)
    [](practical-session-3.html#cb126-10)
    [](practical-session-3.html#cb126-11)ground_truth <- read_table("./data/to_load/spotzonationGroup.txt")

## 3.4 Spot-level Quality Control[](practical-session-3.html#spot-level-quality-control-1)

### 3.4.1 Calculating QC metrics[](practical-session-3.html#calculating-qc-metrics-1)

In this section we are effectively recapitlating the spot- and gene-level QC from practical 2 for this new dataset, in order that we can use it for the practical exercises in the next session.
    
    
    [](practical-session-3.html#cb127-1)is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
    [](practical-session-3.html#cb127-2)sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
    [](practical-session-3.html#cb127-3)sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
    [](practical-session-3.html#cb127-4)sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)
    [](practical-session-3.html#cb127-5)
    [](practical-session-3.html#cb127-6)colData(sfe)
    
    
    ## DataFrame with 1185 rows and 15 columns
    ##                    in_tissue array_row array_col   sample_id            Barcode
    ##                    <logical> <integer> <integer> <character>        <character>
    ## AAACAAGTATCTCCCA-1      TRUE        50       102      JBO019 AAACAAGTATCTCCCA-1
    ## AAACATTTCCCGGATT-1      TRUE        61        97      JBO019 AAACATTTCCCGGATT-1
    ## AAACCCGAACGAAATC-1      TRUE        45       115      JBO019 AAACCCGAACGAAATC-1
    ## AAACGAGACGGTTGAT-1      TRUE        35        79      JBO019 AAACGAGACGGTTGAT-1
    ## AAACTAACGTGGCGAC-1      TRUE         8       110      JBO019 AAACTAACGTGGCGAC-1
    ## ...                      ...       ...       ...         ...                ...
    ## TTGTAATCCGTACTCG-1      TRUE        35        55      JBO019 TTGTAATCCGTACTCG-1
    ## TTGTGAACCTAATCCG-1      TRUE        56        90      JBO019 TTGTGAACCTAATCCG-1
    ## TTGTGCAGCCACGTCA-1      TRUE        60        74      JBO019 TTGTGCAGCCACGTCA-1
    ## TTGTGTTTCCCGAAAG-1      TRUE        51        59      JBO019 TTGTGTTTCCCGAAAG-1
    ## TTGTTGTGTGTCAAGA-1      TRUE        31        77      JBO019 TTGTTGTGTGTCAAGA-1
    ##                      Capt_area  annotation       index  sparsity       sum
    ##                    <character> <character> <character> <numeric> <numeric>
    ## AAACAAGTATCTCCCA-1           1          NA      spot_1  0.910410     13443
    ## AAACATTTCCCGGATT-1           1          NA      spot_2  0.967805      2648
    ## AAACCCGAACGAAATC-1           1         Mid      spot_3  0.864958     27733
    ## AAACGAGACGGTTGAT-1           1     Central      spot_4  0.835818     32973
    ## AAACTAACGTGGCGAC-1           1          NA      spot_5  0.995418       400
    ## ...                        ...         ...         ...       ...       ...
    ## TTGTAATCCGTACTCG-1           1          NA   spot_1181  0.933716      7612
    ## TTGTGAACCTAATCCG-1           1          NA   spot_1182  0.955831      4299
    ## TTGTGCAGCCACGTCA-1           1          NA   spot_1183  0.978252      1452
    ## TTGTGTTTCCCGAAAG-1           1          NA   spot_1184  0.956778      3831
    ## TTGTTGTGTGTCAAGA-1           1         Mid   spot_1185  0.852160     27755
    ##                     detected subsets_mito_sum subsets_mito_detected
    ##                    <integer>        <numeric>             <integer>
    ## AAACAAGTATCTCCCA-1      2933             1021                    12
    ## AAACATTTCCCGGATT-1      1054              285                    12
    ## AAACCCGAACGAAATC-1      4421             2087                    12
    ## AAACGAGACGGTTGAT-1      5375              821                    12
    ## AAACTAACGTGGCGAC-1       150              182                    11
    ## ...                      ...              ...                   ...
    ## TTGTAATCCGTACTCG-1      2170              733                    11
    ## TTGTGAACCTAATCCG-1      1446              515                    12
    ## TTGTGCAGCCACGTCA-1       712               54                    10
    ## TTGTGTTTCCCGAAAG-1      1415              422                    11
    ## TTGTTGTGTGTCAAGA-1      4840              906                    12
    ##                    subsets_mito_percent     total
    ##                               <numeric> <numeric>
    ## AAACAAGTATCTCCCA-1              7.59503     13443
    ## AAACATTTCCCGGATT-1             10.76284      2648
    ## AAACCCGAACGAAATC-1              7.52533     27733
    ## AAACGAGACGGTTGAT-1              2.48992     32973
    ## AAACTAACGTGGCGAC-1             45.50000       400
    ## ...                                 ...       ...
    ## TTGTAATCCGTACTCG-1              9.62953      7612
    ## TTGTGAACCTAATCCG-1             11.97953      4299
    ## TTGTGCAGCCACGTCA-1              3.71901      1452
    ## TTGTGTTTCCCGAAAG-1             11.01540      3831
    ## TTGTTGTGTGTCAAGA-1              3.26428     27755
    
    
    [](practical-session-3.html#cb129-1)rowData(sfe)
    
    
    ## DataFrame with 32738 rows and 18 columns
    ##                    gene_name              id       mean  detected     total
    ##                  <character>     <character>  <numeric> <numeric> <numeric>
    ## ENSG00000243485   MIR1302-10 ENSG00000243485 0.00000000  0.000000         0
    ## ENSG00000237613      FAM138A ENSG00000237613 0.00000000  0.000000         0
    ## ENSG00000186092        OR4F5 ENSG00000186092 0.00000000  0.000000         0
    ## ENSG00000238009 RP11-34P13.7 ENSG00000238009 0.00590717  0.590717         7
    ## ENSG00000239945 RP11-34P13.8 ENSG00000239945 0.00000000  0.000000         0
    ## ...                      ...             ...        ...       ...       ...
    ## ENSG00000215635   AC145205.1 ENSG00000215635          0         0         0
    ## ENSG00000268590        BAGE5 ENSG00000268590          0         0         0
    ## ENSG00000251180   CU459201.1 ENSG00000251180          0         0         0
    ## ENSG00000215616   AC002321.2 ENSG00000215616          0         0         0
    ## ENSG00000215611   AC002321.1 ENSG00000215611          0         0         0
    ##                 JBO019.sparsity JBO019.total JBO019.nLocations JBO019.s_min
    ##                       <numeric>    <numeric>         <integer>    <numeric>
    ## ENSG00000243485        1.000000            0                 0          Inf
    ## ENSG00000237613        1.000000            0                 0          Inf
    ## ENSG00000186092        1.000000            0                 0          Inf
    ## ENSG00000238009        0.994093            7                 7            1
    ## ENSG00000239945        1.000000            0                 0          Inf
    ## ...                         ...          ...               ...          ...
    ## ENSG00000215635               1            0                 0          Inf
    ## ENSG00000268590               1            0                 0          Inf
    ## ENSG00000251180               1            0                 0          Inf
    ## ENSG00000215616               1            0                 0          Inf
    ## ENSG00000215611               1            0                 0          Inf
    ##                 JBO019.max JBO019.s_mean JBO019.s_median JBO019.s_SD
    ##                  <numeric>     <numeric>       <numeric>   <numeric>
    ## ENSG00000243485          0           NaN              NA          NA
    ## ENSG00000237613          0           NaN              NA          NA
    ## ENSG00000186092          0           NaN              NA          NA
    ## ENSG00000238009          1             1               1           0
    ## ENSG00000239945          0           NaN              NA          NA
    ## ...                    ...           ...             ...         ...
    ## ENSG00000215635          0           NaN              NA          NA
    ## ENSG00000268590          0           NaN              NA          NA
    ## ENSG00000251180          0           NaN              NA          NA
    ## ENSG00000215616          0           NaN              NA          NA
    ## ENSG00000215611          0           NaN              NA          NA
    ##                 JBO019.p_mean JBO019.p_median JBO019.p_SD JBO019.s_CV
    ##                     <numeric>       <numeric>   <numeric>   <numeric>
    ## ENSG00000243485    0.00000000               0   0.0000000          NA
    ## ENSG00000237613    0.00000000               0   0.0000000          NA
    ## ENSG00000186092    0.00000000               0   0.0000000          NA
    ## ENSG00000238009    0.00590717               0   0.0766631           0
    ## ENSG00000239945    0.00000000               0   0.0000000          NA
    ## ...                       ...             ...         ...         ...
    ## ENSG00000215635             0               0           0          NA
    ## ENSG00000268590             0               0           0          NA
    ## ENSG00000251180             0               0           0          NA
    ## ENSG00000215616             0               0           0          NA
    ## ENSG00000215611             0               0           0          NA
    ##                 JBO019.p_CV
    ##                   <numeric>
    ## ENSG00000243485         NaN
    ## ENSG00000237613         NaN
    ## ENSG00000186092         NaN
    ## ENSG00000238009      1297.8
    ## ENSG00000239945         NaN
    ## ...                     ...
    ## ENSG00000215635         NaN
    ## ENSG00000268590         NaN
    ## ENSG00000251180         NaN
    ## ENSG00000215616         NaN
    ## ENSG00000215611         NaN
    
    
    [](practical-session-3.html#cb131-1)colGeometries(sfe)
    
    
    ## List of length 3
    ## names(3): spotPoly spotCntd spotHex

### 3.4.2 Plot manual annotation[](practical-session-3.html#plot-manual-annotation)
    
    
    [](practical-session-3.html#cb133-1)ggplot() + 
    [](practical-session-3.html#cb133-2)  geom_sf(aes(geometry = colGeometries(sfe)$spotHex$geometry, fill = colData(sfe)$annotation)) + 
    [](practical-session-3.html#cb133-3)  theme_void() + 
    [](practical-session-3.html#cb133-4)  theme(legend.position = "right") + 
    [](practical-session-3.html#cb133-5)  labs(fill = "Annotation")

![](_main_files/figure-html/03_QC_sfe2-1.png)

### 3.4.3 Library size threshold[](practical-session-3.html#library-size-threshold)
    
    
    [](practical-session-3.html#cb134-1)# ----------------------------------------------- #
    [](practical-session-3.html#cb134-2)## Density and histogram of library sizes
    [](practical-session-3.html#cb134-3)ggplot(data = as.data.frame(colData(sfe)),
    [](practical-session-3.html#cb134-4)       aes(x = sum)) +
    [](practical-session-3.html#cb134-5)    geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-3.html#cb134-6)                   colour = "black", 
    [](practical-session-3.html#cb134-7)                   fill = "grey",
    [](practical-session-3.html#cb134-8)                   bins = 50) +
    [](practical-session-3.html#cb134-9)    geom_density(alpha = 0.5,
    [](practical-session-3.html#cb134-10)                 adjust = 0.5,
    [](practical-session-3.html#cb134-11)                 fill = "#A0CBE8",
    [](practical-session-3.html#cb134-12)                 colour = "#4E79A7") +
    [](practical-session-3.html#cb134-13)    geom_vline(xintercept = c(1000, NA),
    [](practical-session-3.html#cb134-14)               colour = "red", 
    [](practical-session-3.html#cb134-15)               linetype = "dashed") + 
    [](practical-session-3.html#cb134-16)    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb134-17)    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb134-18)    xlab("Library size") + 
    [](practical-session-3.html#cb134-19)    ylab("Density") + 
    [](practical-session-3.html#cb134-20)    theme_classic()
    
    
    ## Warning: Removed 1 rows containing missing values (`geom_vline()`).

![](_main_files/figure-html/03_QC_sfe3-1.png)
    
    
    [](practical-session-3.html#cb136-1)## Select library size threshold
    [](practical-session-3.html#cb136-2)qc_lib_size <- colData(sfe)$sum < 1000 #| colData(sfe)$sum > 45000
    [](practical-session-3.html#cb136-3)## Check how many spots are filtered out
    [](practical-session-3.html#cb136-4)table(qc_lib_size)
    
    
    ## qc_lib_size
    ## FALSE  TRUE 
    ##  1166    19
    
    
    [](practical-session-3.html#cb138-1)## Add threshold in colData
    [](practical-session-3.html#cb138-2)colData(sfe)$qc_lib_size <- qc_lib_size
    [](practical-session-3.html#cb138-3)## Check putative spatial patterns of removed spots
    [](practical-session-3.html#cb138-4)ggplot() + 
    [](practical-session-3.html#cb138-5)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb138-6)            aes(geometry = geometry)) + 
    [](practical-session-3.html#cb138-7)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb138-8)            aes(geometry = geometry, fill = colData(sfe)$qc_lib_size)) +
    [](practical-session-3.html#cb138-9)    scale_fill_manual(values = c("grey95", "red")) + 
    [](practical-session-3.html#cb138-10)    labs(fill = "Discarded") + 
    [](practical-session-3.html#cb138-11)    theme_bw()

![](_main_files/figure-html/03_QC_sfe3-2.png)

### 3.4.4 Number of expressed genes[](practical-session-3.html#number-of-expressed-genes-1)
    
    
    [](practical-session-3.html#cb139-1)# ----------------------------------------------- #
    [](practical-session-3.html#cb139-2)## Density and histogram of expressed genes
    [](practical-session-3.html#cb139-3)ggplot(data = as.data.frame(colData(sfe)),
    [](practical-session-3.html#cb139-4)       aes(x = detected)) +
    [](practical-session-3.html#cb139-5)    geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-3.html#cb139-6)                   colour = "black", 
    [](practical-session-3.html#cb139-7)                   fill = "grey",
    [](practical-session-3.html#cb139-8)                   bins = 50) +
    [](practical-session-3.html#cb139-9)    geom_density(alpha = 0.5,
    [](practical-session-3.html#cb139-10)                 adjust = 0.5,
    [](practical-session-3.html#cb139-11)                 fill = "#A0CBE8",
    [](practical-session-3.html#cb139-12)                 colour = "#4E79A7") + 
    [](practical-session-3.html#cb139-13)    geom_vline(xintercept = c(550, NA),
    [](practical-session-3.html#cb139-14)               colour = "red", 
    [](practical-session-3.html#cb139-15)               linetype = "dashed") +
    [](practical-session-3.html#cb139-16)    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb139-17)    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb139-18)    xlab("Genes expressed in each spot") + 
    [](practical-session-3.html#cb139-19)    ylab("Density") + 
    [](practical-session-3.html#cb139-20)    theme_classic()
    
    
    ## Warning: Removed 1 rows containing missing values (`geom_vline()`).

![](_main_files/figure-html/03_QC_sfe4-1.png)
    
    
    [](practical-session-3.html#cb141-1)## Select expressed genes threshold
    [](practical-session-3.html#cb141-2)qc_detected <- colData(sfe)$detected < 550 #| colData(sfe)$detected > 6000
    [](practical-session-3.html#cb141-3)## Check how many spots are filtered out
    [](practical-session-3.html#cb141-4)table(qc_detected)
    
    
    ## qc_detected
    ## FALSE  TRUE 
    ##  1165    20
    
    
    [](practical-session-3.html#cb143-1)## Add threshold in colData
    [](practical-session-3.html#cb143-2)colData(sfe)$qc_detected <- qc_detected
    [](practical-session-3.html#cb143-3)## Check for putative spatial pattern of removed spots
    [](practical-session-3.html#cb143-4)ggplot() + 
    [](practical-session-3.html#cb143-5)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb143-6)            aes(geometry = geometry)) + 
    [](practical-session-3.html#cb143-7)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb143-8)            aes(geometry = geometry, fill = colData(sfe)$qc_detected)) +
    [](practical-session-3.html#cb143-9)    scale_fill_manual(values = c("grey95", "red")) + 
    [](practical-session-3.html#cb143-10)    labs(fill = "Discarded") + 
    [](practical-session-3.html#cb143-11)    theme_bw()

![](_main_files/figure-html/03_QC_sfe4-2.png)

### 3.4.5 Percentage of mitochondrial expression[](practical-session-3.html#percentage-of-mitochondrial-expression-1)
    
    
    [](practical-session-3.html#cb144-1)# ----------------------------------------------- #
    [](practical-session-3.html#cb144-2)## Density and histogram of percentage of mitochondrial expression
    [](practical-session-3.html#cb144-3)ggplot(data = as.data.frame(colData(sfe)),
    [](practical-session-3.html#cb144-4)       aes(x = subsets_mito_percent)) +
    [](practical-session-3.html#cb144-5)    geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-3.html#cb144-6)                   colour = "black", 
    [](practical-session-3.html#cb144-7)                   fill = "grey",
    [](practical-session-3.html#cb144-8)                   bins = 50) +
    [](practical-session-3.html#cb144-9)    geom_density(alpha = 0.5,
    [](practical-session-3.html#cb144-10)                 adjust = 0.5,
    [](practical-session-3.html#cb144-11)                 fill = "#A0CBE8",
    [](practical-session-3.html#cb144-12)                 colour = "#4E79A7") + 
    [](practical-session-3.html#cb144-13)    geom_vline(xintercept = c(22, NA),
    [](practical-session-3.html#cb144-14)               colour = "red", 
    [](practical-session-3.html#cb144-15)               linetype = "dashed") +
    [](practical-session-3.html#cb144-16)    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb144-17)    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb144-18)    xlab("Percentage of mitochondrial expression") + 
    [](practical-session-3.html#cb144-19)    ylab("Density") + 
    [](practical-session-3.html#cb144-20)    theme_classic()
    
    
    ## Warning: Removed 1 rows containing missing values (`geom_vline()`).

![](_main_files/figure-html/03_QC_sfe5-1.png)
    
    
    [](practical-session-3.html#cb146-1)## Select mitochondrial percentage threshold
    [](practical-session-3.html#cb146-2)qc_mito <- colData(sfe)$subsets_mito_percent > 22
    [](practical-session-3.html#cb146-3)## Check how many spots are filtered out
    [](practical-session-3.html#cb146-4)table(qc_mito)
    
    
    ## qc_mito
    ## FALSE  TRUE 
    ##  1180     5
    
    
    [](practical-session-3.html#cb148-1)## Add threshold in colData
    [](practical-session-3.html#cb148-2)colData(sfe)$qc_mito <- qc_mito
    [](practical-session-3.html#cb148-3)## Check for putative spatial pattern of removed spots
    [](practical-session-3.html#cb148-4)ggplot() + 
    [](practical-session-3.html#cb148-5)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb148-6)            aes(geometry = geometry)) + 
    [](practical-session-3.html#cb148-7)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb148-8)            aes(geometry = geometry, fill = colData(sfe)$qc_mito)) +
    [](practical-session-3.html#cb148-9)    scale_fill_manual(values = c("grey95", "red")) + 
    [](practical-session-3.html#cb148-10)    labs(fill = "Discarded") + 
    [](practical-session-3.html#cb148-11)    theme_bw()

![](_main_files/figure-html/03_QC_sfe5-2.png)

### 3.4.6 Remove low-quality spots[](practical-session-3.html#remove-low-quality-spots-1)
    
    
    [](practical-session-3.html#cb149-1)# ----------------------------------------------- #
    [](practical-session-3.html#cb149-2)## Check the number of discarded spots for each metric
    [](practical-session-3.html#cb149-3)apply(cbind(qc_lib_size, qc_detected, qc_mito), 2, sum)
    
    
    ## qc_lib_size qc_detected     qc_mito 
    ##          19          20           5
    
    
    [](practical-session-3.html#cb151-1)## Combine together the set of discarded spots
    [](practical-session-3.html#cb151-2)discard <- qc_lib_size | qc_detected | qc_mito
    [](practical-session-3.html#cb151-3)table(discard)
    
    
    ## discard
    ## FALSE  TRUE 
    ##  1161    24
    
    
    [](practical-session-3.html#cb153-1)## Store the set in the object
    [](practical-session-3.html#cb153-2)colData(sfe)$discard <- discard
    [](practical-session-3.html#cb153-3)## Check for putative spatial pattern of removed spots
    [](practical-session-3.html#cb153-4)ggplot() + 
    [](practical-session-3.html#cb153-5)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb153-6)            aes(geometry = geometry)) + 
    [](practical-session-3.html#cb153-7)    geom_sf(data = colGeometry(sfe, "spotHex"),
    [](practical-session-3.html#cb153-8)            aes(geometry = geometry, fill = colData(sfe)$discard)) +
    [](practical-session-3.html#cb153-9)    scale_fill_manual(values = c("grey95", "red")) + 
    [](practical-session-3.html#cb153-10)    labs(fill = "Discarded") + 
    [](practical-session-3.html#cb153-11)    theme_bw()

![](_main_files/figure-html/03_QC_sfe6-1.png)
    
    
    [](practical-session-3.html#cb154-1)# ----------------------------------------------- #
    [](practical-session-3.html#cb154-2)## remove combined set of low-quality spots
    [](practical-session-3.html#cb154-3)sfe <- sfe[, !colData(sfe)$discard]

## 3.5 Normalisation of counts[](practical-session-3.html#normalisation-of-counts-1)

### 3.5.1 Log-tranformation of counts[](practical-session-3.html#log-tranformation-of-counts-1)
    
    
    [](practical-session-3.html#cb155-1)## Calculate library size factors
    [](practical-session-3.html#cb155-2)sfe <- computeLibraryFactors(sfe)
    [](practical-session-3.html#cb155-3)## Have a look at the size factors
    [](practical-session-3.html#cb155-4)summary(sizeFactors(sfe))
    
    
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.07961 0.36902 0.95469 1.00000 1.54936 2.77256
    
    
    [](practical-session-3.html#cb157-1)## Density and histogram of library sizes
    [](practical-session-3.html#cb157-2)ggplot(data = data.frame(sFact = sizeFactors(sfe)), 
    [](practical-session-3.html#cb157-3)       aes(x = sFact)) +
    [](practical-session-3.html#cb157-4)    geom_histogram(aes(y = after_stat(density)), 
    [](practical-session-3.html#cb157-5)                   colour = "black", 
    [](practical-session-3.html#cb157-6)                   fill = "grey",
    [](practical-session-3.html#cb157-7)                   bins = 40) +
    [](practical-session-3.html#cb157-8)    geom_density(alpha = 0.5,
    [](practical-session-3.html#cb157-9)                 adjust = 0.5,
    [](practical-session-3.html#cb157-10)                 fill = "#A0CBE8",
    [](practical-session-3.html#cb157-11)                 colour = "#4E79A7") +
    [](practical-session-3.html#cb157-12)    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb157-13)    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    [](practical-session-3.html#cb157-14)    xlab("Library size") + 
    [](practical-session-3.html#cb157-15)    ylab("Density") + 
    [](practical-session-3.html#cb157-16)    theme_classic()

![](_main_files/figure-html/03_LogNorm_sfe-1.png)
    
    
    [](practical-session-3.html#cb158-1)# calculate logcounts using library size factors
    [](practical-session-3.html#cb158-2)sfe <- logNormCounts(sfe)

## 3.6 Gene-level Quality Control[](practical-session-3.html#gene-level-quality-control)

### 3.6.1 Calculating extra QC metrics[](practical-session-3.html#calculating-extra-qc-metrics)
    
    
    [](practical-session-3.html#cb159-1)rowData(sfe)[["JBO019.s_logMean"]] <- rowSums(assay(sfe, "logcounts")) / rowData(sfe)[["JBO019.nLocations"]]

### 3.6.2 Set and apply filters[](practical-session-3.html#set-and-apply-filters)
    
    
    [](practical-session-3.html#cb160-1)is_zero <- rowData(sfe)$total == 0
    [](practical-session-3.html#cb160-2)is_logLow <- rowData(sfe)[["JBO019.s_logMean"]] <= 1
    [](practical-session-3.html#cb160-3)discard_gs <- is_zero | is_mito | is_logLow
    [](practical-session-3.html#cb160-4)table(discard_gs)
    
    
    ## discard_gs
    ## FALSE  TRUE 
    ##  8535 24203
    
    
    [](practical-session-3.html#cb162-1)rowData(sfe)$discard <- discard_gs
    [](practical-session-3.html#cb162-2)
    [](practical-session-3.html#cb162-3)## FEATURE SELECTION
    [](practical-session-3.html#cb162-4)## remove mitochondrial and other genes
    [](practical-session-3.html#cb162-5)sfe <- sfe[!rowData(sfe)$discard, ]

## 3.7 Selecting genes[](practical-session-3.html#selecting-genes-1)

### 3.7.1 Highly Variable Genes (HVGs)[](practical-session-3.html#highly-variable-genes-hvgs-1)
    
    
    [](practical-session-3.html#cb163-1)## Fit mean-variance relationship
    [](practical-session-3.html#cb163-2)dec <- modelGeneVar(sfe,
    [](practical-session-3.html#cb163-3)                    assay.type = "logcounts")
    [](practical-session-3.html#cb163-4)
    [](practical-session-3.html#cb163-5)## Visualize mean-variance relationship
    [](practical-session-3.html#cb163-6)fit <- metadata(dec)
    [](practical-session-3.html#cb163-7)fit_df <- data.frame(mean = fit$mean,
    [](practical-session-3.html#cb163-8)                     var = fit$var,
    [](practical-session-3.html#cb163-9)                     trend = fit$trend(fit$mean))
    [](practical-session-3.html#cb163-10)
    [](practical-session-3.html#cb163-11)ggplot(data = fit_df, 
    [](practical-session-3.html#cb163-12)       aes(x = mean, y = var)) + 
    [](practical-session-3.html#cb163-13)    geom_point() + 
    [](practical-session-3.html#cb163-14)    geom_line(aes(y = trend), colour = "dodgerblue", linewidth = 1.5) + 
    [](practical-session-3.html#cb163-15)    labs(x = "mean of log-expression",
    [](practical-session-3.html#cb163-16)         y = "variance of log-expression") + 
    [](practical-session-3.html#cb163-17)    theme_classic()

![](_main_files/figure-html/03_HVGs_sfe-1.png)
    
    
    [](practical-session-3.html#cb164-1)## Select top HVGs
    [](practical-session-3.html#cb164-2)top_hvgs <- getTopHVGs(dec, 
    [](practical-session-3.html#cb164-3)                       var.field = "bio", 
    [](practical-session-3.html#cb164-4)                       prop = 0.5,
    [](practical-session-3.html#cb164-5)                       var.threshold = 0,
    [](practical-session-3.html#cb164-6)                       fdr.threshold = 0.1)

### 3.7.2 Code for 3.3 to 3.7[](practical-session-3.html#code-for-3.3-to-3.7)
    
    
    [](practical-session-3.html#cb165-1)## Import data
    [](practical-session-3.html#cb165-2)sampleDir <- "./data/spaceranger_outs/Human_Liver_Steatotic/JBO019_Results"
    [](practical-session-3.html#cb165-3)sampleNames <- "JBO019"
    [](practical-session-3.html#cb165-4)sfe <- read10xVisiumSFE(samples = sampleDir, 
    [](practical-session-3.html#cb165-5)                        sample_id = sampleNames, 
    [](practical-session-3.html#cb165-6)                        type = "sparse", 
    [](practical-session-3.html#cb165-7)                        data = "filtered", 
    [](practical-session-3.html#cb165-8)                        images = "lowres", 
    [](practical-session-3.html#cb165-9)                        style = "W", 
    [](practical-session-3.html#cb165-10)                        zero.policy = TRUE)
    [](practical-session-3.html#cb165-11)# ----------------------------------------------- #
    [](practical-session-3.html#cb165-12)ground_truth <- read_table("./data/to_load/spotzonationGroup.txt")
    [](practical-session-3.html#cb165-13)## Add QC metrics
    [](practical-session-3.html#cb165-14)is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
    [](practical-session-3.html#cb165-15)sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
    [](practical-session-3.html#cb165-16)sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
    [](practical-session-3.html#cb165-17)sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)
    [](practical-session-3.html#cb165-18)# ----------------------------------------------- #
    [](practical-session-3.html#cb165-19)## SPOT SELECTION
    [](practical-session-3.html#cb165-20)## Select library size threshold
    [](practical-session-3.html#cb165-21)qc_lib_size <- colData(sfe)$sum < 1000
    [](practical-session-3.html#cb165-22)## Add threshold in colData
    [](practical-session-3.html#cb165-23)colData(sfe)$qc_lib_size <- qc_lib_size
    [](practical-session-3.html#cb165-24)## Select expressed genes threshold
    [](practical-session-3.html#cb165-25)qc_detected <- colData(sfe)$detected < 550
    [](practical-session-3.html#cb165-26)## Add threshold in colData
    [](practical-session-3.html#cb165-27)colData(sfe)$qc_detected <- qc_detected
    [](practical-session-3.html#cb165-28)## Select mitochondrial percentage threshold
    [](practical-session-3.html#cb165-29)qc_mito <- colData(sfe)$subsets_mito_percent > 22
    [](practical-session-3.html#cb165-30)## Add threshold in colData
    [](practical-session-3.html#cb165-31)colData(sfe)$qc_mito <- qc_mito
    [](practical-session-3.html#cb165-32)## Combine together the set of discarded spots
    [](practical-session-3.html#cb165-33)discard <- qc_lib_size | qc_detected | qc_mito
    [](practical-session-3.html#cb165-34)## Store the set in the object
    [](practical-session-3.html#cb165-35)colData(sfe)$discard <- discard
    [](practical-session-3.html#cb165-36)## Remove combined set of low-quality spots
    [](practical-session-3.html#cb165-37)sfe <- sfe[, !colData(sfe)$discard]
    [](practical-session-3.html#cb165-38)# ----------------------------------------------- #
    [](practical-session-3.html#cb165-39)## FEATURE SELECTION
    [](practical-session-3.html#cb165-40)## Calculate library size factors
    [](practical-session-3.html#cb165-41)sfe <- computeLibraryFactors(sfe)
    [](practical-session-3.html#cb165-42)## Calculate logcounts using library size factors
    [](practical-session-3.html#cb165-43)sfe <- logNormCounts(sfe)
    [](practical-session-3.html#cb165-44)## Calculate log-counts sample mean
    [](practical-session-3.html#cb165-45)rowData(sfe)[["JBO019.s_logMean"]] <- rowSums(assay(sfe, "logcounts")) / rowData(sfe)[["JBO019.nLocations"]]
    [](practical-session-3.html#cb165-46)## Set and apply filters
    [](practical-session-3.html#cb165-47)is_zero <- rowData(sfe)$total == 0
    [](practical-session-3.html#cb165-48)is_logLow <- rowData(sfe)[["JBO019.s_logMean"]] <= 1
    [](practical-session-3.html#cb165-49)discard_gs <- is_zero | is_mito | is_logLow
    [](practical-session-3.html#cb165-50)rowData(sfe)$discard <- discard_gs
    [](practical-session-3.html#cb165-51)## Remove mitochondrial and other genes
    [](practical-session-3.html#cb165-52)sfe <- sfe[!rowData(sfe)$discard, ]
    [](practical-session-3.html#cb165-53)## Fit mean-variance relationship
    [](practical-session-3.html#cb165-54)dec <- modelGeneVar(sfe,
    [](practical-session-3.html#cb165-55)                    assay.type = "logcounts")
    [](practical-session-3.html#cb165-56)## Select top HVGs
    [](practical-session-3.html#cb165-57)top_hvgs <- getTopHVGs(dec, 
    [](practical-session-3.html#cb165-58)                       var.field = "bio", 
    [](practical-session-3.html#cb165-59)                       prop = 0.5,
    [](practical-session-3.html#cb165-60)                       var.threshold = 0,
    [](practical-session-3.html#cb165-61)                       fdr.threshold = 0.05)

## 3.8 Neighbour graph and distance matrix[](practical-session-3.html#neighbour-graph-and-distance-matrix)

### 3.8.1 Adding spatial weights[](practical-session-3.html#adding-spatial-weights)

The neighbour lists can be supplemented with spatial weights using the `nb2listw` and `nb2listwdist` function from `spdep` package for the chosen type and coding scheme style. There are 6 different coding scheme styles that can be used to weigh neighbour relationships:

  1. **B** : is the basic binary coding (1 for neighbour, 0 for no neighbour).
  2. **W** : is row standardised (sums over all links to n).
  3. **C** : is globally standardised (sums over all links to n).
  4. **U** : is equal to C divided by the number of neighbours (sums over all links to unity).
  5. **S** : is the variance-stabilizing coding scheme (sums over all links to n).
  6. **minmax** : divides the weights by the minimum of the maximum row sums and maximum column sums of the input weights; It is similar to the C and U styles.



The coding scheme style is practically the value each neighbour will get. For example, in a binary coding scheme style (**B**) if a spot is a neighbour of the spot in focus then gets the value of **1** , else gets **0**. Another example, in a row standardised coding scheme style (**W**) if the spot in focus has a total of 10 neighbours and each neighbour has a weight of 1, then the sum of all neighbour weights is 10, and each neighbour will get a normalised weight of 1/10 = 0.1. As a result, in the row standardised coding scheme, spots with many neighbours will have neighbours with lower weights and thus will not be over-emphasised.

Starting from a binary neighbours list, in which regions are either listed as neighbours or are absent (thus not in the set of neighbours for some definition), we can add a distance-based weights list. The `nb2listwdist` function supplements a neighbours list with spatial weights for the chosen types of distance modelling and coding scheme. While the offered coding schemes parallel those of the `nb2listw` function above, three distance-based types of weights are available: inverse distance weighting (IDW), double-power distance weights (DPD), and exponential distance decay (EXP). The three types of distance weight calculations are based on pairwise distances ððð, all of which are controlled by parameter _âalphaâ_ (ð¼ below):

  1. **idw** : ð¤ðð=ðâð¼ðð,
  2. **exp** : ð¤ðð=exp(âð¼â ððð),
  3. **dpd** : ð¤ðð=[1â(ððð/ðmax)ð¼]ð¼,



the latter of which leads to ð¤ðð=0 for all ððð>ðmax. Note that _IDW_ weights show extreme behaviour close to 0 and can take on the value infinity. In such cases, the infinite values are replaced by the largest finite weight present in the weights list.

### 3.8.2 Generate distance matrices[](practical-session-3.html#generate-distance-matrices)

A distance matrix is a mirrored matrix that contains the distance between a spot and every other spot. This distance can be a simple Euclidean distance based on the coordinates of the spots or a weighted distance according to a bandwidth around each spot using a kernel that gives higher scores to distances between spots that are closer together compared to the ones that are farther away. These weighted distance matrices are later used to run geographically weighted (GW) models.

There are 6 different kernels that can be used to weight the distances between spots. The next two figures are from the `GWmodel` publication (Gollini et al. 2015) and illustrate the mathematical application of these kernals, and show graphically how they weight by distance.

![The math equations that define the kernels.](images/gwmodel_kernel_math.png)

Figure 3.2: The math equations that define the kernels. 

![Examples from using each kernel.](images/gwmodel_kernel_graphs.png)

Figure 3.3: Examples from using each kernel. 

In the below we choose one of the many possible ways of building a neighbour graph for the steatotic liver data set. In this example we are using a k-nearest neighbours approach with row-standardised distance-based weights.
    
    
    [](practical-session-3.html#cb166-1)## add a neighbour graph using a weighted distance matrix
    [](practical-session-3.html#cb166-2)sfe <- addSpatialNeighGraphs(sfe, "JBO019", type = "knearneigh", style = "W", distMod = "raw", k = 6)
    [](practical-session-3.html#cb166-3)
    [](practical-session-3.html#cb166-4)colGraphs(sfe)
    
    
    ## $col
    ## Characteristics of weights list object:
    ## Neighbour list object:
    ## Number of regions: 1161 
    ## Number of nonzero links: 6966 
    ## Percentage nonzero weights: 0.5167959 
    ## Average number of links: 6 
    ## Non-symmetric neighbours list
    ## 
    ## Weights style: W 
    ## Weights constants summary:
    ##      n      nn   S0       S1       S2
    ## W 1161 1347921 1161 376.8333 4674.667
    
    
    [](practical-session-3.html#cb168-1)## Calculate a simple distance matrix
    [](practical-session-3.html#cb168-2)sfe <- addDistMat(sfe, p = 2)

We can use a `geom` from the `tidyterra` package (commonly used for map visualisations) to plot the neighbour graph we generated in the previous step.
    
    
    [](practical-session-3.html#cb169-1)## Retrieve the tissue image
    [](practical-session-3.html#cb169-2)sfei <- getImg(sfe, image_id = "lowres")
    [](practical-session-3.html#cb169-3)## Extract the spot locations
    [](practical-session-3.html#cb169-4)spot_coords <- spatialCoords(sfe) %>% as.data.frame()
    [](practical-session-3.html#cb169-5)
    [](practical-session-3.html#cb169-6)## Set limits
    [](practical-session-3.html#cb169-7)xlim <- c(min(spot_coords$pxl_col_in_fullres) - 100, 
    [](practical-session-3.html#cb169-8)          max(spot_coords$pxl_col_in_fullres) + 100)
    [](practical-session-3.html#cb169-9)ylim <- c(min(spot_coords$pxl_row_in_fullres) - 100, 
    [](practical-session-3.html#cb169-10)          max(spot_coords$pxl_row_in_fullres) + 100)
    [](practical-session-3.html#cb169-11)nbs <- colGraph(sfe)
    [](practical-session-3.html#cb169-12)ggplot() + 
    [](practical-session-3.html#cb169-13)    geom_spatraster_rgb(data = imgRaster(sfei)) + 
    [](practical-session-3.html#cb169-14)    geom_sf(data = as(nb2lines(nbs$neighbours, coords = spatialCoords(sfe)), "sf")) + 
    [](practical-session-3.html#cb169-15)    lims(x = xlim, y = ylim) +
    [](practical-session-3.html#cb169-16)    coord_sf() + 
    [](practical-session-3.html#cb169-17)    theme_void()
    
    
    ## Warning in CRS(proj4string): CRS: projargs should not be NULL; set to NA

![](_main_files/figure-html/03_visualise_neighbours-1.png)

Now that we have a fully QC-ed dataset with spatial weights and a neighbour graph applied, we have prepared our data fully for the application of geospatial methods - specifically in practical 4, geogrpahically weighted principal components analysis (GWPCA).

## 3.9 Putting it all together[](practical-session-3.html#putting-it-all-together-1)

The below code puts all these steps in order by selecting one of the options at each step.
    
    
    [](practical-session-3.html#cb171-1)## Import data
    [](practical-session-3.html#cb171-2)sampleDir <- "./data/spaceranger_outs/Human_Liver_Steatotic/JBO019_Results"
    [](practical-session-3.html#cb171-3)sampleNames <- "JBO019"
    [](practical-session-3.html#cb171-4)sfe <- read10xVisiumSFE(samples = sampleDir, 
    [](practical-session-3.html#cb171-5)                        sample_id = sampleNames, 
    [](practical-session-3.html#cb171-6)                        type = "sparse", 
    [](practical-session-3.html#cb171-7)                        data = "filtered", 
    [](practical-session-3.html#cb171-8)                        images = "lowres", 
    [](practical-session-3.html#cb171-9)                        style = "W", 
    [](practical-session-3.html#cb171-10)                        zero.policy = TRUE)
    [](practical-session-3.html#cb171-11)# ----------------------------------------------- #
    [](practical-session-3.html#cb171-12)ground_truth <- read_table("./data/to_load/spotzonationGroup.txt")
    [](practical-session-3.html#cb171-13)## Add QC metrics
    [](practical-session-3.html#cb171-14)is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
    [](practical-session-3.html#cb171-15)sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
    [](practical-session-3.html#cb171-16)sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
    [](practical-session-3.html#cb171-17)sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)
    [](practical-session-3.html#cb171-18)# ----------------------------------------------- #
    [](practical-session-3.html#cb171-19)## SPOT SELECTION
    [](practical-session-3.html#cb171-20)## Select library size threshold
    [](practical-session-3.html#cb171-21)qc_lib_size <- colData(sfe)$sum < 1000
    [](practical-session-3.html#cb171-22)## Add threshold in colData
    [](practical-session-3.html#cb171-23)colData(sfe)$qc_lib_size <- qc_lib_size
    [](practical-session-3.html#cb171-24)## Select expressed genes threshold
    [](practical-session-3.html#cb171-25)qc_detected <- colData(sfe)$detected < 550
    [](practical-session-3.html#cb171-26)## Add threshold in colData
    [](practical-session-3.html#cb171-27)colData(sfe)$qc_detected <- qc_detected
    [](practical-session-3.html#cb171-28)## Select mitochondrial percentage threshold
    [](practical-session-3.html#cb171-29)qc_mito <- colData(sfe)$subsets_mito_percent > 22
    [](practical-session-3.html#cb171-30)## Add threshold in colData
    [](practical-session-3.html#cb171-31)colData(sfe)$qc_mito <- qc_mito
    [](practical-session-3.html#cb171-32)## Combine together the set of discarded spots
    [](practical-session-3.html#cb171-33)discard <- qc_lib_size | qc_detected | qc_mito
    [](practical-session-3.html#cb171-34)## Store the set in the object
    [](practical-session-3.html#cb171-35)colData(sfe)$discard <- discard
    [](practical-session-3.html#cb171-36)## Remove combined set of low-quality spots
    [](practical-session-3.html#cb171-37)sfe <- sfe[, !colData(sfe)$discard]
    [](practical-session-3.html#cb171-38)# ----------------------------------------------- #
    [](practical-session-3.html#cb171-39)## FEATURE SELECTION
    [](practical-session-3.html#cb171-40)## Calculate library size factors
    [](practical-session-3.html#cb171-41)sfe <- computeLibraryFactors(sfe)
    [](practical-session-3.html#cb171-42)## Calculate logcounts using library size factors
    [](practical-session-3.html#cb171-43)sfe <- logNormCounts(sfe)
    [](practical-session-3.html#cb171-44)## Calculate log-counts sample mean
    [](practical-session-3.html#cb171-45)rowData(sfe)[["JBO019.s_logMean"]] <- rowSums(assay(sfe, "logcounts")) / rowData(sfe)[["JBO019.nLocations"]]
    [](practical-session-3.html#cb171-46)## Set and apply filters
    [](practical-session-3.html#cb171-47)is_zero <- rowData(sfe)$total == 0
    [](practical-session-3.html#cb171-48)is_logLow <- rowData(sfe)[["JBO019.s_logMean"]] <= 1
    [](practical-session-3.html#cb171-49)discard_gs <- is_zero | is_mito | is_logLow
    [](practical-session-3.html#cb171-50)rowData(sfe)$discard <- discard_gs
    [](practical-session-3.html#cb171-51)## Remove mitochondrial and other genes
    [](practical-session-3.html#cb171-52)sfe <- sfe[!rowData(sfe)$discard, ]
    [](practical-session-3.html#cb171-53)## Fit mean-variance relationship
    [](practical-session-3.html#cb171-54)dec <- modelGeneVar(sfe,
    [](practical-session-3.html#cb171-55)                    assay.type = "logcounts")
    [](practical-session-3.html#cb171-56)## Select top HVGs
    [](practical-session-3.html#cb171-57)top_hvgs <- getTopHVGs(dec, 
    [](practical-session-3.html#cb171-58)                       var.field = "bio", 
    [](practical-session-3.html#cb171-59)                       prop = 0.5,
    [](practical-session-3.html#cb171-60)                       var.threshold = 0,
    [](practical-session-3.html#cb171-61)                       fdr.threshold = 0.05)
    [](practical-session-3.html#cb171-62)# ----------------------------------------------- #
    [](practical-session-3.html#cb171-63)## ADD GEOGRAPHY
    [](practical-session-3.html#cb171-64)## Add a neighbour graph using a weighted distance matrix
    [](practical-session-3.html#cb171-65)sfe <- addSpatialNeighGraphs(sfe, "JBO019", type = "knearneigh", style = "W", distMod = "raw", k = 6)
    [](practical-session-3.html#cb171-66)## Calculate a simple distance matrix
    [](practical-session-3.html#cb171-67)sfe <- addDistMat(sfe, p = 2)

### References[](references.html#references)

Gollini, Isabella, Binbin Lu, Martin Charlton, Christopher Brunsdon, and Paul Harris. 2015. âGWmodel: An R Package for Exploring Spatial Heterogeneity Using Geographically Weighted Models.â _J Stat Soft_ 63 (February): 1â50. <https://doi.org/10.18637/jss.v063.i17>. 

Guilliams, Martin, Johnny Bonnardel, Birthe Haest, Bart Vanderborght, Camille Wagner, Anneleen Remmerie, Anna Bujko, et al. 2022. âSpatial Proteogenomics Reveals Distinct and Evolutionarily Conserved Hepatic Macrophage Niches.â _Cell_ 185 (2): 379â396.e38. https://doi.org/<https://doi.org/10.1016/j.cell.2021.12.018>. 


<!-- PAGE: practical-session-4.html -->

# Chapter 4 Practical session 4[](practical-session-4.html#practical-session-4)

In this session we will have a hands-on exploration of GW-PCA and its application to STx data. What can we learn from this novel technique?

## 4.1 Geographically Weighted Principal Components Analysis (GWPCA)[](practical-session-4.html#geographically-weighted-principal-components-analysis-gwpca)

A standard PCA can pick out the key multivariate modes of variability in the data. Looking at outlying values of the principal components of these data gives us an idea of unusual sites (in terms of combinations of gene expression profiles - and to a certain extent of combinations of cell types in each spot). Next, geographically weighted PCA can be used to find spatial multivariate outliers. Sounds complicated, but really all this means is it identifies sites that have an unusual multi-way combination of gene expression in relation to their immediate geographical neighbours. It might be that the values observed at these sites as a combination is not uncommon in the tissue as a whole - but is very unusual in its locality.

To find such outliers the procedure is relatively simple - instead of doing a PCA on the tissue as a whole, for each sample we do a PCA on data falling into a window centred on the location of that spot. In that way we can check whether the spot is like its neighbours or not, from a multivariate viewpoint.

The procedure we will follow in this practical carries out a geographically weighted PCA. In short, it runs a âwindowedâ PCA around each of the spots.

## 4.2 Load packages[](practical-session-4.html#load-packages-1)

## 4.3 Load Quality Controled and Normalised data[](practical-session-4.html#load-quality-controled-and-normalised-data)

First of all, we need to load the data we prepared in the previous practical.
    
    
    [](practical-session-4.html#cb172-1)sfe <- readRDS(file = "./data/to_load/practical03_sfe.rds")
    [](practical-session-4.html#cb172-2)top_hvgs <- readRDS(file = "./data/to_load/practical03_topHVGs.rds")

## 4.4 Parameter prearation for GWPCA[](practical-session-4.html#parameter-prearation-for-gwpca)

The `gwpca` method uses `princomp` internally to run the PCAs - this function does not allow the number of variables (genes) to be greater than the number of samples (spots). This imposes a hard requirement on the data pre-processing. We have, however, already identified the highly variable genes in our sample, and for this case, there are fewer genes than spots.

Some other parameterisation is neccessary and these required parameters (as we have used for this dataset) are illustrated here:
    
    
    [](practical-session-4.html#cb173-1)## Get the gene names that are going to be evaluated
    [](practical-session-4.html#cb173-2)vars = top_hvgs
    [](practical-session-4.html#cb173-3)## Set a fixed bandwidth
    [](practical-session-4.html#cb173-4)bw = 6*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]
    [](practical-session-4.html#cb173-5)## Set the number of components to be retained
    [](practical-session-4.html#cb173-6)k = 20
    [](practical-session-4.html#cb173-7)## Set the kernel to be used
    [](practical-session-4.html#cb173-8)kernel = "gaussian"
    [](practical-session-4.html#cb173-9)## Set the Minkowski distance power: p = 2 --> Euclidean
    [](practical-session-4.html#cb173-10)p = 2
    [](practical-session-4.html#cb173-11)## Is the bandwidth adaptive?: No because spots are fixed
    [](practical-session-4.html#cb173-12)adaptive = FALSE
    [](practical-session-4.html#cb173-13)## Cross-Validate GWPCA?
    [](practical-session-4.html#cb173-14)cv = TRUE
    [](practical-session-4.html#cb173-15)## Calculate PCA scores?
    [](practical-session-4.html#cb173-16)scores = FALSE
    [](practical-session-4.html#cb173-17)## Run a robust GWPCA?
    [](practical-session-4.html#cb173-18)robust = FALSE
    [](practical-session-4.html#cb173-19)## Make a cluster for parallel computing (otherwise GWPCA is slow!)
    [](practical-session-4.html#cb173-20)my.cl <- parallel::makeCluster(parallelly::availableCores() - 1, type = 'FORK')

The bandwidth defines a radius around each spot - every spot that falls inside this radius is considered a neighbour. We can set bandwidth as a fixed value (as here) or we can select the bandwidth automatically. Without going into detail here, this is achieved by a form of cross validation, where each observation is omitted, and it is attempted to reconstruct the values on the basis of principal components, derived from the other observations. The bandwidth achieving the optimal results is the one selected. For a complete explanation, see Harris, Brunsdon, and Charlton (2011). The function `bw.gwpca` from `GWmodel` can be used to computes this.

  * **NOTE** : Larger bandwidths imply bigger moving spatial windows, which in turn imply smoother spatially varying outputs.



## 4.5 Run GWPCA[](practical-session-4.html#run-gwpca)

Here we present the invocation to run GWPCA, however because this process is computationally intensive and time-consuming, we do not suggest running it on posit.cloud. We have pre-computed the result and provide it for you to load.
    
    
    [](practical-session-4.html#cb174-1)# DO NOT RUN THIS CHUNK
    [](practical-session-4.html#cb174-2)
    [](practical-session-4.html#cb174-3)pcagw <- gwpcaSTE(sfe = sfe, 
    [](practical-session-4.html#cb174-4)                  assay = "logcounts",
    [](practical-session-4.html#cb174-5)                  vars = vars, 
    [](practical-session-4.html#cb174-6)                  p = p, 
    [](practical-session-4.html#cb174-7)                  k = k, 
    [](practical-session-4.html#cb174-8)                  bw = bw, 
    [](practical-session-4.html#cb174-9)                  kernel = kernel,
    [](practical-session-4.html#cb174-10)                  adaptive = adaptive, 
    [](practical-session-4.html#cb174-11)                  scores = scores, 
    [](practical-session-4.html#cb174-12)                  robust = robust,
    [](practical-session-4.html#cb174-13)                  cv = cv,
    [](practical-session-4.html#cb174-14)                  future = FALSE,
    [](practical-session-4.html#cb174-15)                  strategy = "cluster",
    [](practical-session-4.html#cb174-16)                  workers = my.cl,
    [](practical-session-4.html#cb174-17)                  verbose = FALSE)
    [](practical-session-4.html#cb174-18)saveRDS(pcagw, file = "./data/to_load/practical04_pcagw.rds")

Because GWPCA can take some time to run, we ran it for you and below you can load the output:
    
    
    [](practical-session-4.html#cb175-1)pcagw <- readRDS(file = "./data/to_load/practical04_pcagw.rds")

## 4.6 Plot global PCA results[](practical-session-4.html#plot-global-pca-results)

In the next steps we will take a look inside the output from the `gwpca` function and we are going to extract some basic information. Since GWPCA consists of multiple local PCAs, it is good to know how many PCs makes sense to look at. We can do so by running a global PCA and plotting a scree plot:
    
    
    [](practical-session-4.html#cb176-1)plotGWPCA_global(gwpca = pcagw,
    [](practical-session-4.html#cb176-2)                 comps = 1:10,
    [](practical-session-4.html#cb176-3)                 type = "scree",
    [](practical-session-4.html#cb176-4)                 point_args = list(size = 3, colour = "red"),
    [](practical-session-4.html#cb176-5)                 line_args = list(linewidth = 1, colour = "dodgerblue"))

![](_main_files/figure-html/04_scree_plot-1.png)

In a Principal Component Analysis (PCA), the first three principal components may explain less than 15% of the variance in the data if the data is highly dispersed or if there is a large amount of noise in the data. This means that the first three principal components are not capturing a significant portion of the variability in the data. This could be due to a lack of clear structure in the data or a lack of meaningful patterns that can be captured by the PCA. Alternatively, it could be due to the presence of many irrelevant features or variables in the data that are not contributing to the overall variance. This is one more of the reasons why GWPCA is more appropriate for STx data. Because, it may be true that the global PCs are not strong but locally this can change.

## 4.7 Identify the leading genes in each location[](practical-session-4.html#identify-the-leading-genes-in-each-location)

The genes with the highest loading scores (where loading score = correlation between variable and component) at each location can be thought of as the âleading genesâ - i.e.Â those with the most explanatory power with respect to the variability of gene expression at that location. These leading genes can be a local indicator of relevant biology.

Here we look at leading genes in 2 ways - (1) by finding the single gene with the highest loading at each location; (2) by finding sets of the top 4 genes by loading score, where the order of those genes does not matter (so the ordered set A,B,C,D is considered the same as D,B,A,C).
    
    
    [](practical-session-4.html#cb177-1)## Extract leading genes
    [](practical-session-4.html#cb177-2)pcagw <- gwpca_LeadingGene(gwpca = pcagw, 
    [](practical-session-4.html#cb177-3)                           sfe = sfe, 
    [](practical-session-4.html#cb177-4)                           pc_nos = 1:4, 
    [](practical-session-4.html#cb177-5)                           type = "single", 
    [](practical-session-4.html#cb177-6)                           names = "gene_names")
    
    
    ## 16  leading genes found for  PC1
    ## The leading genes in  PC1  are:
    ##     ADH1A        C7       CRP    CYP3A4      GLUL     GSTA2      HAMP      HBA2 
    ##         2        11         4       365         7         1        13        33 
    ##     IGLL5    MALAT1 MTRNR2L12  MTRNR2L8      NNMT     PTGDS      SAA1       SDS 
    ##        87        39       153       181        23        73        36       133 
    ## 21  leading genes found for  PC2
    ## The leading genes in  PC2  are:
    ##        C7       CAT     CFHR1       CRP    CYP3A4      GLUL      HBA2       HBB 
    ##         3         6        38        39       149        83         2        37 
    ##    IGFBP3    IGFBP7       IGJ     IGLL5    MALAT1 MTRNR2L10 MTRNR2L12  MTRNR2L8 
    ##        49        39        34       246        80        10        78       124 
    ##      NNMT      SAA1       SDS     TAGLN    UGT2B7 
    ##        42        12        69        20         1 
    ## 24  leading genes found for  PC3
    ## The leading genes in  PC3  are:
    ##     AEBP1        C7       CAT     CFHR1       CRP    CYP3A4      GLUL      HBA2 
    ##         2         2        27        20         5        20        17        27 
    ##       HBB    IGFBP3    IGFBP7       IGJ     IGLL5    MALAT1 MTRNR2L10 MTRNR2L12 
    ##       150        41        77         6       399       136         6        61 
    ##  MTRNR2L8      MYL9      NNMT      SAA1   SCGB3A1       SDS     TAGLN    UGT2B7 
    ##        25         9        24         6        56        15        26         4 
    ## 25  leading genes found for  PC4
    ## The leading genes in  PC4  are:
    ##     AEBP1       CAT     CFHR1       CRP     FXYD2      GLUL     GSTA2      HBA2 
    ##         1        53        15         7         7        33         3         2 
    ##       HBB    IGFBP3    IGFBP7       IGJ     IGLL5    MALAT1 MTRNR2L10 MTRNR2L12 
    ##       181       100        51        60       281       201         5        16 
    ##  MTRNR2L8      MYL9      NNMT      ORM2      SAA1       SDS    SPINK1     TAGLN 
    ##        16         5        55         6         6        37        12         4 
    ##    UGT2B7 
    ##         4
    
    
    [](practical-session-4.html#cb179-1)pcagw <- gwpca_LeadingGene(gwpca = pcagw, 
    [](practical-session-4.html#cb179-2)                           sfe = sfe, 
    [](practical-session-4.html#cb179-3)                           pc_nos = 1:4, 
    [](practical-session-4.html#cb179-4)                           genes_n = 4, 
    [](practical-session-4.html#cb179-5)                           type = "multi", 
    [](practical-session-4.html#cb179-6)                           method = "membership", 
    [](practical-session-4.html#cb179-7)                           names = "gene_names")
    
    
    ## The number of individual leading genes groups found for PC1 is: 110 
    ## These groups are: Too many to print them!
    ## The number of individual leading genes groups found for PC2 is: 240 
    ## These groups are: Too many to print them!
    ## The number of individual leading genes groups found for PC3 is: 310 
    ## These groups are: Too many to print them!
    ## The number of individual leading genes groups found for PC4 is: 421 
    ## These groups are: Too many to print them!

We can also plot these leading genes on the spot map - as each location by definition has (potentially) a different leading gene.
    
    
    [](practical-session-4.html#cb181-1)## Plot leading genes
    [](practical-session-4.html#cb181-2)plotGWPCA_leadingG(gwpca = pcagw,
    [](practical-session-4.html#cb181-3)                   comps = 1:2,
    [](practical-session-4.html#cb181-4)                   type = "single",
    [](practical-session-4.html#cb181-5)                   arrange = FALSE)
    [](practical-session-4.html#cb181-6)
    [](practical-session-4.html#cb181-7)plotGWPCA_leadingG(gwpca = pcagw,
    [](practical-session-4.html#cb181-8)                   comps = 1,
    [](practical-session-4.html#cb181-9)                   type = "multi",
    [](practical-session-4.html#cb181-10)                   arrange = FALSE)

![](_main_files/figure-html/leading_genes2-1.png)![](_main_files/figure-html/leading_genes2-2.png)![](_main_files/figure-html/leading_genes2-3.png) The âmultiâ plot here is problematic, because there are too many groups of genes to be able to print a legible legend. The alternative below is provided to highlight gene groups that are found in at least 12 spots.
    
    
    [](practical-session-4.html#cb182-1)### Plot multi type (alternative)
    [](practical-session-4.html#cb182-2)## The data
    [](practical-session-4.html#cb182-3)leadingGsMulti <- pcagw$leadingGeneMulti
    [](practical-session-4.html#cb182-4)## The Legend labels
    [](practical-session-4.html#cb182-5)spot_labels <- data.frame(table(leadingGsMulti[1])) %>%
    [](practical-session-4.html#cb182-6)    dplyr::rename(LeadingGs = colnames(leadingGsMulti)[1], 
    [](practical-session-4.html#cb182-7)                  count = Freq) %>%
    [](practical-session-4.html#cb182-8)    dplyr::arrange(desc(count)) %>%
    [](practical-session-4.html#cb182-9)    mutate(show = ifelse(count > 12, TRUE, FALSE))
    [](practical-session-4.html#cb182-10)    
    [](practical-session-4.html#cb182-11)## The legend breaks:
    [](practical-session-4.html#cb182-12)spot_breaks <- spot_labels %>%
    [](practical-session-4.html#cb182-13)    dplyr::filter(show == TRUE) %>% 
    [](practical-session-4.html#cb182-14)    dplyr::arrange(LeadingGs) %>%
    [](practical-session-4.html#cb182-15)    dplyr::select(LeadingGs) %>% 
    [](practical-session-4.html#cb182-16)    .[["LeadingGs"]] %>% 
    [](practical-session-4.html#cb182-17)    as.vector()
    [](practical-session-4.html#cb182-18)    
    [](practical-session-4.html#cb182-19)## The colours:
    [](practical-session-4.html#cb182-20)col_No <- sum(spot_labels$show)
    [](practical-session-4.html#cb182-21)colour_values <- getColours(col_No)
    [](practical-session-4.html#cb182-22)names(colour_values) <- spot_labels$LeadingGs[spot_labels$show]
    [](practical-session-4.html#cb182-23)pc <- "PC1"
    [](practical-session-4.html#cb182-24)    
    [](practical-session-4.html#cb182-25)## The Plot:
    [](practical-session-4.html#cb182-26)ggplot() + 
    [](practical-session-4.html#cb182-27)    geom_sf(data = leadingGsMulti, 
    [](practical-session-4.html#cb182-28)            aes(geometry = geometry$geometry,
    [](practical-session-4.html#cb182-29)                fill = .data[[pc]]),
    [](practical-session-4.html#cb182-30)            colour = "grey30", 
    [](practical-session-4.html#cb182-31)            show.legend = TRUE) + 
    [](practical-session-4.html#cb182-32)    scale_fill_manual(values = colour_values,
    [](practical-session-4.html#cb182-33)                      breaks = spot_breaks,
    [](practical-session-4.html#cb182-34)                      na.value = "gray95") +
    [](practical-session-4.html#cb182-35)    labs(title = NULL,
    [](practical-session-4.html#cb182-36)         fill = "Group of\nLeading\nGenes") + 
    [](practical-session-4.html#cb182-37)    theme_void() +
    [](practical-session-4.html#cb182-38)    theme(legend.position = "bottom", legend.text = element_text(size=6)) +
    [](practical-session-4.html#cb182-39)    guides(fill = guide_legend(ncol = 3, byrow = TRUE))

![](_main_files/figure-html/leading_genes3-1.png)

## 4.8 Percentage of Total Variation (PTV)[](practical-session-4.html#percentage-of-total-variation-ptv)

Another useful diagnostic for PCA is the percentage of variability in the data explained by each of the components. Locally, this can be achieved by looking at the `local.PV` component of `pcagw`; this is written as `pcagw$local.PV`. This is an 1161 by 20 matrix - where 1161 is the number of observations and 20 is the number of components (`k`). For each location, the 20 columns correspond to the percentage of the total variance explained by each of the principal components at that location. If, say, the first two components contributed 90% of the total variance, then it is reasonable to assume that much of the variability in the data can be seen by just looking at these two components. Because this is geographically weighted PCA, this quantity varies across the map.
    
    
    [](practical-session-4.html#cb183-1)## Calculate the PTV for multiple Components
    [](practical-session-4.html#cb183-2)pcagw <- gwpca_PropVar(gwpca = pcagw, n_comp = 2:10, sfe = sfe)
    
    
    ##     Comps_01         Comps_02        Comps_03        Comps_04    
    ##  Min.   : 6.279   Min.   :11.67   Min.   :16.43   Min.   :20.69  
    ##  1st Qu.: 9.483   1st Qu.:16.13   1st Qu.:21.24   1st Qu.:25.69  
    ##  Median :16.782   Median :25.54   Median :30.37   Median :34.46  
    ##  Mean   :17.370   Mean   :25.92   Mean   :31.35   Mean   :35.49  
    ##  3rd Qu.:22.534   3rd Qu.:32.87   3rd Qu.:39.42   3rd Qu.:43.81  
    ##  Max.   :38.254   Max.   :46.50   Max.   :54.25   Max.   :57.51  
    ##     Comps_05        Comps_06        Comps_07        Comps_08    
    ##  Min.   :24.64   Min.   :28.28   Min.   :31.49   Min.   :34.26  
    ##  1st Qu.:29.65   1st Qu.:33.13   1st Qu.:36.54   1st Qu.:39.53  
    ##  Median :37.79   Median :40.86   Median :43.53   Median :46.17  
    ##  Mean   :38.98   Mean   :42.07   Mean   :44.84   Mean   :47.38  
    ##  3rd Qu.:47.17   3rd Qu.:49.78   3rd Qu.:52.16   3rd Qu.:54.19  
    ##  Max.   :60.60   Max.   :62.97   Max.   :65.04   Max.   :67.03  
    ##     Comps_09        Comps_10    
    ##  Min.   :36.76   Min.   :39.15  
    ##  1st Qu.:42.34   1st Qu.:45.05  
    ##  Median :48.60   Median :50.96  
    ##  Mean   :49.73   Mean   :51.91  
    ##  3rd Qu.:56.07   3rd Qu.:57.77  
    ##  Max.   :68.83   Max.   :70.39
    
    
    [](practical-session-4.html#cb185-1)## Plot PTV
    [](practical-session-4.html#cb185-2)plotGWPCA_ptv(gwpca = pcagw,
    [](practical-session-4.html#cb185-3)              comps = 1:10,
    [](practical-session-4.html#cb185-4)              type = "violin")
    [](practical-session-4.html#cb185-5)
    [](practical-session-4.html#cb185-6)## Map PTV
    [](practical-session-4.html#cb185-7)plotGWPCA_ptv(gwpca = pcagw,
    [](practical-session-4.html#cb185-8)              comps = 1:6,
    [](practical-session-4.html#cb185-9)              type = "map")

![](_main_files/figure-html/04_ptv-1.png)![](_main_files/figure-html/04_ptv-2.png)

## 4.9 Identify discrepancies[](practical-session-4.html#identify-discrepancies)

Global PCA can be used to identify multivariate outliers. Extending this, it is also possible to use local PCA (i.e., GWPCA) to identify local outliers. One way of doing this links back to the cross-validation idea that can be used to select a bandwidth. Recall that this is based on a score of how well each observation can be reconstructed on the basis of local PCs. The score measures the total discrepancies of true data values from the reconstructed ones - and the bandwidth chosen is the one minimising this. However, the total discrepancy score is the sum of the individual discrepancies. A very large individual discrepancy associated with an observation suggests it is very different - in a multidimensional way, to the observations near to it.
    
    
    [](practical-session-4.html#cb186-1)## Plot the discrepancies as boxplot
    [](practical-session-4.html#cb186-2)plotGWPCA_discr(pcagw, type = "box")

![](_main_files/figure-html/04_discrep1-1.png)
    
    
    [](practical-session-4.html#cb187-1)## Plot the discrepancies map
    [](practical-session-4.html#cb187-2)plotGWPCA_discr(pcagw, type = "map")

![](_main_files/figure-html/04_discrep2-1.png)
    
    
    [](practical-session-4.html#cb188-1)## Get location data for the discrepancies
    [](practical-session-4.html#cb188-2)discrepancy_loc_dt <- getDiscrepancyLocData(sfe = sfe, 
    [](practical-session-4.html#cb188-3)                                            gwpca = pcagw, 
    [](practical-session-4.html#cb188-4)                                            sample_id = "JBO019")

Another possibility to understand the nature of the outlier is a parallel coordinates heatmap. Here, each observation neighbouring the location that has been found to be an outlier is shown as a column with the genes in rows. Since here we are investigating local outliers, one particular observation is highlighted in red - the outlier - and the remaining ones in grey, but with the intensity of the grey fading according to their distance from the red observation. This enables you to see what characteristic the red observation has that means it as outlying from its neighbours. The plot can be created using `STExplorerDev::plotGWPCA_discrHeatmap`:
    
    
    [](practical-session-4.html#cb189-1)head(discrepancy_loc_dt)
    
    
    ##                              barcodes coords.pxl_col_in_fullres
    ## AAGTGCCTTGACTGTA-1 AAGTGCCTTGACTGTA-1                     11086
    ## ACCCGGATGACGCATC-1 ACCCGGATGACGCATC-1                      9908
    ## ACCTCCGTTATTCACC-1 ACCTCCGTTATTCACC-1                      9113
    ## AGATGATGGAGTCTGG-1 AGATGATGGAGTCTGG-1                      9117
    ## AGGTATAATTGATAGT-1 AGGTATAATTGATAGT-1                      9312
    ## AGTGAACAAACTTCTC-1 AGTGAACAAACTTCTC-1                     11088
    ##                    coords.pxl_row_in_fullres discScore
    ## AAGTGCCTTGACTGTA-1                      5148  22793.56
    ## ACCCGGATGACGCATC-1                      5607  24165.26
    ## ACCTCCGTTATTCACC-1                      4255  27035.77
    ## AGATGATGGAGTCTGG-1                      4933  27423.00
    ## AGGTATAATTGATAGT-1                      4593  23839.25
    ## AGTGAACAAACTTCTC-1                      5374  25111.79
    ##                                          geometry
    ## AAGTGCCTTGACTGTA-1 POLYGON ((11020.03 5034.788...
    ## ACCCGGATGACGCATC-1 POLYGON ((9841.741 5494.291...
    ## ACCTCCGTTATTCACC-1 POLYGON ((9046.741 4142.291...
    ## AGATGATGGAGTCTGG-1 POLYGON ((9050.741 4820.291...
    ## AGGTATAATTGATAGT-1 POLYGON ((9245.076 4480.294...
    ## AGTGAACAAACTTCTC-1 POLYGON ((11020.91 5261.585...
    
    
    [](practical-session-4.html#cb191-1)focus <- discrepancy_loc_dt$barcodes[1:2]
    [](practical-session-4.html#cb191-2)bw = 3*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]
    [](practical-session-4.html#cb191-3)
    [](practical-session-4.html#cb191-4)# Plot the heatmap to visualise the genes that make this location an outlier
    [](practical-session-4.html#cb191-5)plotGWPCA_discrHeatmap(sfe = sfe,
    [](practical-session-4.html#cb191-6)                       assay = "logcounts",
    [](practical-session-4.html#cb191-7)                       vars = NULL,
    [](practical-session-4.html#cb191-8)                       focus = focus,
    [](practical-session-4.html#cb191-9)                       dMetric = "euclidean", 
    [](practical-session-4.html#cb191-10)                       sample_id = "JBO019",
    [](practical-session-4.html#cb191-11)                       bw = bw, 
    [](practical-session-4.html#cb191-12)                       mean.diff = 1, 
    [](practical-session-4.html#cb191-13)                       show.vars = "top", 
    [](practical-session-4.html#cb191-14)                       scale = "row", 
    [](practical-session-4.html#cb191-15)                       gene.names = TRUE,
    [](practical-session-4.html#cb191-16)                       color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000)),
    [](practical-session-4.html#cb191-17)                       fontsize_row = 3)

![](_main_files/figure-html/04_discrep4-1.png)![](_main_files/figure-html/04_discrep4-2.png)
    
    
    [](practical-session-4.html#cb192-1)discrepancy_gene_dt <- getDiscrepancyGeneData(sfe = sfe,
    [](practical-session-4.html#cb192-2)                                              assay = "logcounts",
    [](practical-session-4.html#cb192-3)                                              vars = NULL,
    [](practical-session-4.html#cb192-4)                                              focus = focus[2],
    [](practical-session-4.html#cb192-5)                                              dMetric = "euclidean", 
    [](practical-session-4.html#cb192-6)                                              sample_id = "JBO019",
    [](practical-session-4.html#cb192-7)                                              bw = bw, 
    [](practical-session-4.html#cb192-8)                                              mean.diff = 1, 
    [](practical-session-4.html#cb192-9)                                              show.vars = "top",
    [](practical-session-4.html#cb192-10)                                              exportExpression = TRUE)
    [](practical-session-4.html#cb192-11)head(discrepancy_gene_dt)
    
    
    ##                 AACCCTACTGTCAATA-1 ACATGGCGCCAAAGTA-1 ACCCGGATGACGCATC-1
    ## ENSG00000078808          1.9800764          1.7845077           0.000000
    ## ENSG00000157916          1.1230410          2.3434121           0.000000
    ## ENSG00000171603          0.6681445          0.8599719           1.667949
    ## ENSG00000162496          1.1230410          0.8599719           0.000000
    ## ENSG00000074964          1.4683664          0.0000000           1.964768
    ## ENSG00000158828          1.1230410          0.8599719           0.000000
    ##                 ACGATCATCTTGTAAA-1 AGGGTCGATGCGAACT-1 ATAGTTCCACCCACTC-1
    ## ENSG00000078808           1.434549          0.0000000          0.0000000
    ## ENSG00000157916           1.094451          1.1491034          0.7197148
    ## ENSG00000171603           0.000000          0.0000000          0.0000000
    ## ENSG00000162496           1.094451          2.2182483          1.5560955
    ## ENSG00000074964           1.709571          0.6860569          0.7197148
    ## ENSG00000158828           1.434549          2.0159260          1.1976848
    ##                 ATATCAACCTACAGAG-1 CATCTTACACCACCTC-1 CCATCTCACCAGTGAA-1
    ## ENSG00000078808           1.848204           2.383137          0.7918769
    ## ENSG00000157916           1.356318           1.489110          2.2192853
    ## ENSG00000171603           0.000000           0.000000          1.3002121
    ## ENSG00000162496           2.214350           0.000000          2.2192853
    ## ENSG00000074964           1.028715           1.140629          0.7918769
    ## ENSG00000158828           2.632509           1.140629          0.7918769
    ##                 CCGATCTCAACCTTAT-1 CGCACGTGCGCTATCA-1 CGCTAGAGACCGCTGC-1
    ## ENSG00000078808          1.1456336          2.0109592           1.070644
    ## ENSG00000157916          1.7760473          2.4713829           1.070644
    ## ENSG00000171603          0.0000000          0.0000000           0.000000
    ## ENSG00000162496          0.6836662          1.7105588           1.678406
    ## ENSG00000074964          1.1456336          0.8136275           1.070644
    ## ENSG00000158828          0.6836662          1.3307021           1.070644
    ##                 CTAGGTCTGAAGGAAT-1 GAAAGAACAGCGTTAT-1 GCAGACCCAGCACGTA-1
    ## ENSG00000078808          1.1275395          2.0533150          2.1052992
    ## ENSG00000157916          0.0000000          1.5312450          0.8684635
    ## ENSG00000171603          0.6712291          0.7049757          0.0000000
    ## ENSG00000162496          1.1275395          2.5947724          1.4067734
    ## ENSG00000074964          1.4736763          1.5312450          1.4067734
    ## ENSG00000158828          1.1275395          1.5312450          0.8684635
    ##                 GCGCAAGAGCGCGCTG-1 GCTCGCTCATGTCCAA-1 GCTGTTGCTACCGAAC-1
    ## ENSG00000078808           2.505794          1.4756527           2.004239
    ## ENSG00000157916           1.117631          1.8764253           1.588731
    ## ENSG00000171603           0.000000          0.0000000           1.325320
    ## ENSG00000162496           2.173035          2.1897262           1.588731
    ## ENSG00000074964           0.000000          0.9188039           1.002827
    ## ENSG00000158828           1.739765          0.9188039           1.811392
    ##                 GGTTCTACTCGTCTGA-1 GTGCGACAGGGAGTGT-1 TAGACTACCTAGCGTT-1
    ## ENSG00000078808          1.7977169           1.633005           1.421361
    ## ENSG00000157916          1.4066078           1.633005           0.000000
    ## ENSG00000171603          0.8683433           0.000000           0.000000
    ## ENSG00000162496          2.3583475           1.263669           2.123256
    ## ENSG00000074964          0.8683433           2.170771           0.000000
    ## ENSG00000158828          1.4066078           1.633005           1.695032
    ##                 TAGGGAGCTTGGGATG-1 TATAGATGGTCGCAGT-1 TATATCCCTGGGAGGA-1
    ## ENSG00000078808          2.2316954          0.0000000          0.7102708
    ## ENSG00000157916          0.7985659          1.9335702          2.6061791
    ## ENSG00000171603          0.0000000          0.0000000          0.0000000
    ## ENSG00000162496          0.7985659          1.9335702          1.1841040
    ## ENSG00000074964          0.7985659          0.9560551          0.7102708
    ## ENSG00000158828          0.0000000          1.5260667          1.1841040
    ##                 TATTACCATCCTGCTT-1 TCATTTAAGTCTCCGA-1 TCTGGGAACCTTTGAA-1
    ## ENSG00000078808          1.9382158          2.1150669          2.3926955
    ## ENSG00000157916          1.0928336          1.2220876          1.4967274
    ## ENSG00000171603          0.6475105          0.0000000          0.6846746
    ## ENSG00000162496          2.4677235          1.5846069          2.0131754
    ## ENSG00000074964          0.6475105          1.2220876          0.6846746
    ## ENSG00000158828          1.0928336          0.7367522          1.4967274
    ##                 TGACATCGAGCGGACC-1 TGACGATGCACTAGAA-1 TTGAATTCACGTGAGG-1
    ## ENSG00000078808          1.9076262          2.4843308          1.8133018
    ## ENSG00000157916          1.2485083          1.3402294          2.0507006
    ## ENSG00000171603          0.7552962          1.0152672          0.0000000
    ## ENSG00000162496          2.7012885          1.3402294          2.0507006
    ## ENSG00000074964          1.2485083          0.5951586          0.7036449
    ## ENSG00000158828          1.6153854          0.0000000          1.1745524
    ##                 TTGACGCTCCATGAGC-1 gene_name
    ## ENSG00000078808           2.296259      SDF4
    ## ENSG00000157916           1.851153      RER1
    ## ENSG00000171603           0.000000    CLSTN1
    ## ENSG00000162496           1.851153     DHRS3
    ## ENSG00000074964           0.000000 ARHGEF10L
    ## ENSG00000158828           1.851153     PINK1

## 4.10 Final Summary[](practical-session-4.html#final-summary)

In this practical we have shown the utility of a geospatial method, GWPCA, to explore the variability of an STx dataset at the local level. By assessing features of the output of this method, we can learn things about the spatial distribution of biologically relevant gene expression.

Hopefully this, alongside the other practicals today, have given you a basic grounding in how to work with STx data and some of the practical considerations of doing so. Although we have demonstrated all of these methods with 10X Genomics [Visium](https://www.10xgenomics.com/products/spatial-gene-expression) data, there is no reason why they are not applicable to any STx method, such as [Slide-Seq](https://curiobioscience.com/) or [Stereo-Seq](https://bgi-australia.com.au/stomics).

You can learn more about the application of geospatial methods to this liver dataset by coming to see our poster: **B-122** in Poster Session B - Tuesday, July 25, between 18:00 CEST and 19:00 CEST.

### References[](references.html#references)

Harris, Paul, Chris Brunsdon, and Martin Charlton. 2011. âGeographically weighted principal components analysis.â _International Journal of Geographical Information Science_ 25 (10): 1717â36. <https://doi.org/10.1080/13658816.2011.554838>. 


<!-- PAGE: references.html -->

# References[](references.html#references)

Amezquita, Robert A., Aaron T. L. Lun, Etienne Becht, Vince J. Carey, Lindsay N. Carpp, Ludwig Geistlinger, Federico Marini, et al. 2020. âOrchestrating single-cell analysis with Bioconductor.â _Nat Methods_ 17 (February): 137â45. <https://doi.org/10.1038/s41592-019-0654-x>. 

Gollini, Isabella, Binbin Lu, Martin Charlton, Christopher Brunsdon, and Paul Harris. 2015. âGWmodel: An R Package for Exploring Spatial Heterogeneity Using Geographically Weighted Models.â _J Stat Soft_ 63 (February): 1â50. <https://doi.org/10.18637/jss.v063.i17>. 

Guilliams, Martin, Johnny Bonnardel, Birthe Haest, Bart Vanderborght, Camille Wagner, Anneleen Remmerie, Anna Bujko, et al. 2022. âSpatial Proteogenomics Reveals Distinct and Evolutionarily Conserved Hepatic Macrophage Niches.â _Cell_ 185 (2): 379â396.e38. https://doi.org/<https://doi.org/10.1016/j.cell.2021.12.018>. 

Hao, Yuhan, Stephanie Hao, Erica Andersen-Nissen, William M. Mauck, Shiwei Zheng, Andrew Butler, Maddie J. Lee, et al. 2021. âIntegrated analysis of multimodal single-cell data.â _Cell_ 184 (13): 3573â3587.e29. <https://doi.org/10.1016/j.cell.2021.04.048>. 

Harris, Paul, Chris Brunsdon, and Martin Charlton. 2011. âGeographically weighted principal components analysis.â _International Journal of Geographical Information Science_ 25 (10): 1717â36. <https://doi.org/10.1080/13658816.2011.554838>. 

Keogh, Eamonn, and Abdullah Mueen. 2017. âCurse of Dimensionality.â In _Encyclopedia of Machine Learning and Data Mining_ , 314â15. Boston, MA, USA: Springer, Boston, MA. <https://doi.org/10.1007/978-1-4899-7687-1_192>. 

Li, Yijun, Stefan Stanojevic, Bing He, Zheng Jing, Qianhui Huang, Jian Kang, and Lana X. Garmire. 2022. âBenchmarking Computational Integration Methods for Spatial Transcriptomics Data.â _bioRxiv_ , January, 2021.08.27.457741. <https://doi.org/10.1101/2021.08.27.457741>. 

Lun, Aaron T. L., Davis J. McCarthy, and John C. Marioni. 2016. âA step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.â _F1000Research_ 5 (2122): 2122. <https://doi.org/10.12688/f1000research.9501.2>. 

Maaten, Laurens van der, and Geoffrey Hinton. 2008. âVisualizing Data Using t-SNE.â _Journal of Machine Learning Research_ 9 (86): 2579â2605. <http://jmlr.org/papers/v9/vandermaaten08a.html>. 

Maynard, Kristen R., Leonardo Collado-Torres, Lukas M. Weber, Cedric Uytingco, Brianna K. Barry, Stephen R. Williams, Joseph L. Catallini, et al. 2021. âTranscriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex.â _Nat Neurosci_ 24 (March): 425â36. <https://doi.org/10.1038/s41593-020-00787-0>. 

McCarthy, Davis J., Kieran R. Campbell, Aaron T. L. Lun, and Quin F. Wills. 2017. âScater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R.â _Bioinformatics_ 33 (8): 1179â86. <https://doi.org/10.1093/bioinformatics/btw777>. 

McInnes, Leland, John Healy, and James Melville. 2018. âUMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.â _arXiv_ , February. <https://doi.org/10.48550/arXiv.1802.03426>. 

âNotes on Continuous Stochastic Phenomena on JSTOR.â 1950\. _Biometrika_. <https://www.jstor.org/stable/2332142>. 

Righelli, Dario, Lukas M. Weber, Helena L. Crowell, Brenda Pardo, Leonardo Collado-Torres, Shila Ghazanfar, Aaron T. L. Lun, Stephanie C. Hicks, and Davide Risso. 2022. âSpatialExperiment: infrastructure for spatially-resolved transcriptomics data in R using Bioconductor.â _Bioinformatics_ 38 (11): 3128â31. <https://doi.org/10.1093/bioinformatics/btac299>. 

Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 2020. âStatistical analysis of spatial expression patterns for spatially resolved transcriptomic studies.â _Nat Methods_ 17 (February): 193â200. <https://doi.org/10.1038/s41592-019-0701-7>. 

Svensson, Valentine, Sarah A. Teichmann, and Oliver Stegle. 2018. âSpatialDE: identification of spatially variable genes.â _Nat Methods_ 15 (May): 343â46. <https://doi.org/10.1038/nmeth.4636>. 

âThe Contiguity Ratio and Statistical Mapping on JSTOR.â 1954\. _Incorporated Statistician_. <https://www.jstor.org/stable/2986645>. 

Weber, Lukas M., and Helena L. Crowell. 2022. _Ggspavis: Visualization Functions for Spatially Resolved Transcriptomics Data_. <https://github.com/lmweber/ggspavis>. 

Zhu, Jiaqiang, Shiquan Sun, and Xiang Zhou. 2021. âSPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies.â _Genome Biol_ 22 (1): 1â25. <https://doi.org/10.1186/s13059-021-02404-0>. 
