# Expression of O-GlcNac bound promoters

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  

## Description

Violin plot showing the mouse ESC RNA-seq expression levels of O-GlcNAc enriched gene promoters (236 genes). The median expression of genes having O-GlcNac at promoters is higher when compared to 236 randomly selected promoters without O-GlcNac, 236 randomly selected promoters, and 21,085 promoters. The three groups to which O-GlcNac promoters are compared to have at least 1 RNA-seq read.

The statistical test used is a two-sided Mann-Whitney test (mu = 0, paired = FALSE) giving the p-value:

glcprom vs noglcprom: 0.000692837371518593
glcprom vs randomprom: 0.0356250808138736
glcprom vs allprom: 0.000158394211256996

## Data

Download the following data:

```
#!/bin/bash

mkdir data

# O-GlcNac peaks
wget XXX/ESCHGGlcNAc_rep1.gff -P data/

# RNA-seq counts and feature lengths
wget https://zenodo.org/records/12793186/files/ESCRNAseq_SRR11294181counts.txt -P data/
wget https://zenodo.org/records/12793186/files/ESCRNAseq_SRR11294181countslength.txt -P data/
```

## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1C.yml](fig1C.yml), run:

```
conda env create -n fig1c --file ./fig1C.yml
conda activate fig1c
```

## Figure Generation

Run the script [promexpression.R](promexpression.R) from the current folder, it uses the files downloaded in the subfolders `data/`. From the terminal:

```
Rscript promexpression.R
```

It should output:

```
Reading gene lengths
Reading counts and computing TPM
Retrieve the promoters coordinates
Filtering chromosomes
Reading the OGlcNac peaks
         Processing data/ESCHGGlcNAc_rep1.gff
Perform overlap between O-GlcNac peaks and promoters
# Attempt 1/5 # Connection to Ensembl ...
Connected with success.
# Attempt 1/5 # Retrieving information about genes from biomaRt ...
Information retrieved with success.
Making counts list with only expressed one
Processing rep1
Writing file
Saving 7 x 7 in image
Computing mann-whitney-wilcoxson on each violin
Writing file
```

You should obtain the raw figure:

