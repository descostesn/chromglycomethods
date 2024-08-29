# Genome Browser Screenshot

I. [Description](#description)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  


## Description

Genome browser view of O-GlcNac and RNA-seq signal in control and siOgt Knock-Down (KD) conditions. Both genes show downregulation following Ogt KD.

## Data

The bigwig files used in this figure can be obtained with:

```
#!/bin/bash

## The screenshots were performed with rep 2
wget https://zenodo.org/records/13444099/files/TPMbw_sictrl_rep2.bw
wget https://zenodo.org/records/13444099/files/TPMbw_siogt_rep2.bw

## Replicate 1 is also provided
wget https://zenodo.org/records/13444099/files/TPMbw_sictrl_rep1.bw
wget https://zenodo.org/records/13444099/files/TPMbw_siogt_rep1.bw
```


