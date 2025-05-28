## InterSpatial

InterSpatial is a R package for augmenting cell-cell communication (CCC) analysis for scRNA-seq/snRNA-seq with a low resolution spatial transcriptomic data. Using gene expressions from these different data modalities, InterSpatial is able to detect CCC through ligand-receptor interaction



## Installation

The InterSpatial R package can be installed using devtools:

```r
devtools::install_github("jichunxie/InterSpatial")

```
### Installation of Dependencies

Make sure you have installed the following dependencies with correct version

- The package is built and tested using [Seurat(<= 4.4.0)](https://satijalab.org/seurat/). Users are recommended to use this versions, the output may be altered if using newer Seurat versions.

- Install [HGC](https://www.bioconductor.org/packages/devel/bioc/html/HGC.html) package for Hierarchical Graph Clustering using the following command:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install the HGC package from Bioconductor
BiocManager::install("HGC")
```

## Tutorials

Tutorials for using InterSpatial can be found [here](https://htmlpreview.github.io/?https://github.com/tuhin-majumder/InterSpatial/blob/main/vignettes/interspatial_tutorial.html).
