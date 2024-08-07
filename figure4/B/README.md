# RNA Polymerase II and O-GlcNac co-localize at gene promoters in human DLD-1 cells

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [ChIP-seq and CutnRun](#cutnrun)  

## Description

In [fig1](../../figure1/E/README.md), we observed that the potential candidates underlying the O-GlcNac signal are related to the RNA Polymerase II (RNA Pol II). We thought investigating further the relationship between O-GlcNac and the transcriptional machinery by removing RNA Pol II. We took advantage of an existing system in Human colon adenocarcinoma DLD-1 cells. These cells express OsTIR and a cassette encoding mini-AID (mAID) and fluorescent protein mClover (mAID+mClover) at the initiation site of the endogenous Rpb1 gene locus (POLR2A) (Nagashima 2019). Stimulation by doxocyclin enables to remove RNAPol II.

The heatmap of binding values of RNA Polymerase II and O-GlcNac shows a global colocalization at the 21,519 GRCh38 Ensembl genes. Note that upon RNAPol II removal by doxocyclin induction, O-GlcNac signal remains globally unaffected. As conveyed in the following analysis, we hypothetize that O-GlcNac might be maintained by the recruitment of other O-GlcNacylated proteins at gene promoters.

## Data

```
#!/bin/bash

mkdir data

## The bigwig files of RNAPol II and O-GlcNac before and after induction
wget XXX/DLD1GlcNAcDoxAux_rep1.bw -P data/
wget XXX/DLD1GlcNAcDoxAux_rep2.bw -P data/
wget XXX/DLD1GlcNAcNoDoxAux_rep1.bw -P data/
wget XXX/DLD1GlcNAcNoDoxAux_rep2.bw -P data/
wget https://zenodo.org/records/12793186/files/RNAPolII_SRX11070611_control.bw -P data/
wget https://zenodo.org/records/12793186/files/RNAPolII_SRX11070613_auxin.bw -P data/

## The gene annotations
wget https://zenodo.org/records/12793186/files/Homo_sapiens.GRCh38.100.chrfiltered.tar.gz
```

## Installation
