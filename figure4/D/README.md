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

## The coordinates of the peaks sorted in 5 groups
wget https://zenodo.org/records/12793186/files/peakscoord-fig4C.bed

## The file of the union of the peaks
wget https://zenodo.org/records/12793186/files/union_OGlcNac_noauxaux-fig4C.bed -P data

```

## Installation

***!!!!!! TO FINISH !!!!! ***
Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig4D.yml](fig4D.yml), run:

```
conda env create -n fig4d --file ./fig4D.yml
conda activate fig4d
```

## Figure generation

Run the script:

```
Rscript scripts/enrichmentanalysis.R
```

The script should output:

***!!!!!!!!!!!!!!!!!!! TO MODIFY NB OF ANNOTATIONS REMOVED !!!!!!!!!!!!!!!!!!!!!!!***
```
Reading gff files and return conversion table
Defining background
### No background provided
Retrieving info from biomart
Connecting to biomart
# Attempt 1/5 # Connection to Ensembl ... 
Connected with success.
	 Retrieving gene info
# Attempt 1/5 # Retrieving information about genes from biomaRt ...
Information retrieved with success.
Filtering out non-canonical chromosomes from genesinfo

!!!!!!!!!!!!!!!!!!! TO MODIFY !!!!!!!!!!!!!!!!!!!!!!!


	 Removing 70/57186 annotations with non canonical chromosomes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Creating the input list of entrezID
Creating the entrezID-symbol table
Performing clusters comparison on molecular function
Output the dotplot of the comparison into results/
```

You should obtain the raw figure:

***!!!!!!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!***


### Pre-processing


Perform the genomic repartition of the union of the glc peaks for each replicate:

```
Rscript genomicRepartition.R
```

The script should output:

```
Filtering database chromosomes
Building list of repeats (this might take a while)
     Processing data/LINE.gff
	 Processing data/LTR.gff
	 Processing data/SINE.gff
	 Processing data/Satellite.gff
	 Processing data/DNA.gff
	 Processing data/Simple_repeat.gff
	 Processing data/RNA.gff
	 Processing data/Low_complexity.gff

```


The steps are:

src/R/conf_generation/genomicRepartition_RepClasses_Grouping/glcsept2023_polII_human_enhancers_unionofpeaks.R
src/R/tools/detailsGroupsHeatmapHumanSept2023.R
src/R/conf_generation/venndiagrams_twoexp/glcpolIIhumansept2023/humanpolIIglcnac_heatmapclusters.R

The files with the ordered peaks coordinates `peakscoord-fig4C.bed` was obtained following the methods of [fig4C](../C/README.md).