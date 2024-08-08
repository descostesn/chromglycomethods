# Geneset Enrichment of the three clusters presenting different O-GlcNac occupancy outcomes upon RNA Polymerase II degradation

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  


## Description

Geneset enrichment analysis of the 5 clusters defined in [fig4C](../C/README.md). Genes overlapping one or more O-GlcNac peaks were retrieved as input. Cluster 1 and 3 did not give any enrichment. Cluster 2 (239 genes/584 peaks) is mainly enriched for genes coding for histone methyltransferases and transcription co-activators. Cluster 4 (347 genes/1,200 peaks) gained O-GlcNac on genes coding for proteins linked to the RNA Polymerase II activity. Cluster 5 (751 genes/2,986 peaks) gave an enrichment for general transcription initiation factor activity. The alteration of O-GlcNac occupancy on hundreds of genes coding for proteins participating to the transcriptional activity of the cells is in line with our observation of the relocalization and alteration in numbers of transcription factories (fig3-FG). These observations suggest a role of O-GlcNac in a feedback loop mechanism important for regulating the transcriptional architecture of the genome.

## Data

```
#!/bin/bash

mkdir data

wget https://zenodo.org/records/12793186/files/fig4D_genes_cluster1.gff  -P data/
wget https://zenodo.org/records/12793186/files/fig4D_genes_cluster2.gff  -P data/
wget https://zenodo.org/records/12793186/files/fig4D_genes_cluster3.gff  -P data/
wget https://zenodo.org/records/12793186/files/fig4D_genes_cluster4.gff -P data/
wget https://zenodo.org/records/12793186/files/fig4D_genes_cluster5.gff  -P data/
```

## Installation


