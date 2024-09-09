# RNA Polymerase II and O-GlcNac co-localize at gene promoters in mouse ES cells

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Data](#data-1)  
&nbsp;&nbsp; V.II. [Workflow](#workflow)  

## Description

The heatmap of binding values of RNA Polymerase II and O-GlcNac shows a global colocalization at the 55,634 Refseq genes.

## Data

```
#!/bin/bash

mkdir data

## Bigwig files used to generate the heatmap
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep1.bw  -P data/
wget https://zenodo.org/records/12793186/files/RNApolymeraseII_SRX8556273.bw  -P data/

## Refseq annotations
wget https://zenodo.org/records/12793186/files/refGeneUCSC-mm10-March2021.gff -P annotations/
```

## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1E.yml](fig1E.yml), run:

```
conda env create -n fig1E --file ./fig1E.yml
conda activate fig1E
```
