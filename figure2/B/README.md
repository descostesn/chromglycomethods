# Differentially expressed genes upon degradation of OGT by siRNA

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [CutnRun](#cutnrun)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.II. [RNA-seq](#rna-seq)  
&nbsp;&nbsp; V.II. [DEseq2](#deseq2)  
&nbsp;&nbsp; V.III. [Merging replicates](#merging-replicates)  
&nbsp;&nbsp; V.IV. [DEG overlap O-GlcNac](#deg-overlap-o-glcnac)  


## Description

    Volcano plot indicating the differentially expressed genes upon degradation of OGT by siRNA in mouse ES cells.
We found 867 down regulated (dark blue) and 594 upregulated genes (dark red). Among these, 44 down- and 12 up-regulated genes have a fold-change higher than two. 33 down- (light blue) and 2 up-regulated (orange) genes had an O-GlcNac peaks. Note that the differences with the numbers obtained with the script below are due to the conversion from ensembl to symbol.

## Data

```
#!/bin/bash

mkdir data

# The tabular file given by DEseq2
wget https://zenodo.org/records/12793186/files/resultDeseq2_siogt.txt  -P data/

## The down- and up-regulated genes having an O-GlcNac peak
wget https://zenodo.org/records/12793186/files/siogtdown_withOGlcNac.gff -P data/
wget https://zenodo.org/records/12793186/files/siogtup_withOGlcNac.gff  -P data/
```

## Installation

TO DO

## Figure generation

Run the command:

```
Rscript volcanosiogt.R
```

The script should output:

```
Reading files
Defining status
The number of genes per category is:
down    downwithglc     nodiff      up      upwithglc
836          31         53940       592           2
Plotting volcano
```

You should obtain the raw figure:

<img src="volcano-siogt.png" alt="volcano plot DEG siOgt" width="400"/>


## Pre-processing

### Workflows
#### CutnRun
#### RNA-seq
### DEseq2
### Merging replicates
### DEG overlap O-GlcNac

