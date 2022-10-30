# Unique DUOX2+ACE2+ small cholangiocytes are pathogenic targets for primary biliary cholangitis


## System requirements

These source code written manily by R. And these should be run on Systems support R. R 3.5.1 was use for test these code. R packages denpendencies include: Seurat(version 3.1.2), harmony(version 1.0),clusterProfiler(version 4.4.4). cellphoneDB (version 2.1.0) was used for the cell-cell interaction analysis.

## Installation guide

Please install these software or R package following the guide in their home page or github.

* R: https://cran.r-project.org/
* Seurat: https://satijalab.org/seurat/
* clusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
* harmony: https://github.com/immunogenomics/harmony
* cellphoneDB: https://github.com/Teichlab/cellphonedb

## Demo

* Seurat is used for single-cell RNA sequencing data analysis. Mainly work flow contains QC, sample aggregating, dimension reduction, clustering and visualization. The expected outputs include cell cluster, gene expression matrix and different type visualization. You can find tutorial in the home page of Seurat: https://satijalab.org/seurat/vignettes.html
* harmony was mainly used to process batch effect. The result of harmony analysis will be used in the following dimension reduction step. You can find tutorial in https://github.com/immunogenomics/harmony/blob/master/vignettes/.
* clusterProfiler performs enrichment analysis based on GO or KEGG database. The output of this package including enriched functions and related statistics. You can follow the demo: http://yulab-smu.top/clusterProfiler-book/
* cellphoneDB performs cellullar communication analysis and the visualization of the result. You can get the strength of interactions and the related statistics by cellphoneDB. The understanding of result can refer to https://github.com/Teichlab/cellphonedb/blob/master/Docs/RESULTS-DOCUMENTATION.md. And you can follow the tutorial: https://github.com/Teichlab/cellphonedb.

## Instructions for use

### `main_analysis.R`

This script was main script for data analysis including integrated, cell-cell interaction and enrichment analysis.

#### Data integration

##### Input

* `exp_list` The path of 10X expression matrix corresponding to the sample.


##### Output

* `umap.cluster.pdf` UMAP plot of the integrated data.

* `.diffgenes.xls` The differential expressioned genes list of each cluster.

* `integrated_cluster.rds` The Seurat object after integrated.

#### Annotation

##### Input

* `celltypes` A vector storing cell type infomation.

##### Output

* `umap.label.pdf` Umap plot after add cell type label.

* `integrated_add_celllabel.rds` The Seurat object after annotation and add cell label.

#### Visializtion

##### Input

* `features` Genes need to display.

#### Cell-cell interaction

Cell-cell interaction analysis performed by cellphoneDB. 

##### Input

* `inputdata` The path of Seurat object rds used for cellphoneDB.

#### Enrichment analysis

Enrichment analysis was performed by R package clusterprofiler.

##### Input 

* `diffgene` Gene list used for enrichment analysis.

* `GO/KEGG_term` A vector storing description of terms you want to display.

##### Output

* `Enrichment_GO/KEGG.xls` Enrichment analysis result.

* `Enrichment_GO/KEGG.pdf` The dot plot of the term you want to display in result.

### `sig_score.R`

This script was used to calculate Sig score and display the result.

#### Input
* `clusterrds` Integrated Seruat object rds before annotation.

* `marker` Marker gene list with the first colum was celltype and the second cloumn was the marker genes corresponding to celltype.

* `cluster` A table whit the first column was celltype, and the second column was cluster number corresponding to celltype.

#### Output
* `_sig_score_dotplot.pdf` The dot plot of the result.

* `_data.xls` Sig Score result.

