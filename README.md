# BioCluster
Identify biologically meaningful subclusters

==========
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Tutorial](#Tutorial)
* [Citation](#Citation)

<a name="Introduction"/>

# Introduction

BioCluster is a R package to identify biologically relevant substructures from single-cell RNAseq data
<a name="Installation"/>

# Installation

```
devtools::install_github("liuqivandy/BioCluster")
```



<a name="Tutorial"/>

# Tutorial
```R
library(BioCluster)
library(data.table)
## load the processed seurat object lungobj1 generated from lung single-cell RNAseq data (GSE)
load(url("https://www.dropbox.com/s/c3veuuuu9dk73nw/lung_dropseq1.Rdata?dl=1"))
## read the bulk lung RNAseq from GTEx
lung_bulk<-data.frame(fread("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_tpm/gene_tpm_2017-06-05_v8_lung.gct.gz"),row.names=1)
lung_bulk<-lung_bulk[!duplicated(lung_bulk[,2]),]
rownames(lung_bulk)<-lung_bulk[,2]
lung_bulk<-lung_bulk[,c(-1,-2)]

## evaluate the biological reproducibility of each subcluster of lungobj1 in the GTEx bulk lung expression
lungobj1<-BioCluster(lungobj1,lung_bulk)
## Plot the BioCluster result
PlotBioCluster(lungobj1)
```
<p align="center">
  <img width="800"  src="https://github.com/liuqivandy/BioCluster/blob/master/dimplot.png">
</p>

```R
## Prune the clustering by removing splits that are not supported by the external dataset
pruneres<-BioCluster_prune(lungobj1)
##plot the final clustering result
DimPlot(pruneres$obj,label=T)+NoLegend()
```
<p align="center">
  <img width="800"  src="https://github.com/liuqivandy/BioCluster/dimplot.png">
</p>

```R
## plot the pruned tree
pruneres$ggplot
```
<p align="center">
  <img width="800"  src="https://github.com/liuqivandy/BioCluster/pruneresult">
</p>

<a name="Citation"/>

# Citation
Qi Liu, Yu Shyr. Identification of biologically relevant substructures from single-cell transcriptomics.
