# Geneset enrichment analysis of differentially expressed genes in siOGT

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  


## Description

On the 867 down regulated and 594 upregulated genes (see [fig2B](../B/README.md)), 779 down- and 545 up-regulated genes were found to be enriched in functional categories. On the 33 down-regulated genes having an O-GlcNac peak, 31 were found in functional categories. Among the unexpected categories, down-regulated genes appear to be coding for the structural componants of ribosomes. Strikingly, altered genes of the three groups (down, up, down-glc) also appear to be coding for RNA Polymerase II-specific DNA-binding transcription factors.


## Data


```
#!/bin/bash

mkdir data
mkdir results

## The down- and up-regulated genes obtained with DESeq2
wget https://zenodo.org/records/12793186/files/log0_siogtdown-ensembl.gff -P data/
wget https://zenodo.org/records/12793186/files/log0_siogtup-ensembl.gff -P data/

## The down-regulated genes having an O-GlcNac peak
wget https://zenodo.org/records/12793186/files/siogtdown_withOGlcNac.gff -P data/

## The results given by cluster profiler (not needed to perform the analysis)
wget https://zenodo.org/records/12793186/files/fig2C_clusterprofiler.txt -P results/
```

## Figure generation

The script should output:

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
	 Removing 70/57186 annotations with non canonical chromosomes

```