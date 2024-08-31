# Genome browser view of genes from cluster 4 and 5

I. [Description](#description)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  

## Description

Examples of RNA Pol II and O-GlcNAc enrichment upon Dox-Aux treatment by genome browser view. Two examples of genes at cluster 4 (ARMC5) and cluster 5 (TAF8) are reported to show respectively the enrichment and loss of O-GlcNAc signal upon RNA Pol II degradation.

## Data

```
#!/bin/bash

mkdir data

## The bigwigs of O-GlcNca before and after treatment by (Dox)/Auxin
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcNoDoxAux_rep2.bw -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcDoxAux_rep2.bw -P data/

## The bigwigs of RNAPol II before and after treatment by (Dox)/Auxin
wget https://zenodo.org/records/12793186/files/RNAPolII_SRX11070611_control.bw -P data/
wget https://zenodo.org/records/12793186/files/RNAPolII_SRX11070613_auxin.bw -P data/

## The replicate 1 is also available
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcNoDoxAux_rep1.bw -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcDoxAux_rep1.bw -P data/
```

## Figure Generation

The bw files were uploaded to [IGV](https://igv.org/) v2.13.0 selecting the mouse mm10 genome.


## Pre-processing

See [fig. 4D](../D/README.md#pre-processing)