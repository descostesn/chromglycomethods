# Heatmap of O-GlcNac peaks along with ChIP-Atlas candidates

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [ChIP-seq and CutnRun](#cutnrun)  
&nbsp;&nbsp; V.II. [Union](#union)  

!! Perform the union of peaks before computing the heatmap, show a venn diagram in the doc
!! Delete union_sept2023mouse_HG.gff from zenodo and upload the new one called union_sept2023mouse_HG1-2.gff, replace the name in the data section
!! precise number of peaks in the description
!! Define a single heatmap with 3 groups for pol II
!! Delete mouseESC_fig1E_peakorder.bed on zenodo and upload the new file with same name

## Description

Heatmap of CuntRun O-GlcNac signal and ChIP-Seq ChIP-Atlas selected candidates signal 1 kb around the start of XXX O-GlcNac peaks. The list of peaks was obtained computing the union of the replicate peaks. Three groups were defined and highlight the RNA Polymerase II (RNAPol II) before, at, and after the O-GlcNac peaks. For the groups I and II, RNAPol II describes a double peak signal overlapping O-GlcNac independently of its position. Other selected candidates follow a similar pattern. It indicates that a sub-category of these proteins might carry O-GlcNac in a transient state. The description of the function of each protein can be found in [fig1D](../D/README.md).

## Data

Download the following data:

```
#!/bin/bash

mkdir data

# The bed file containing the sorted peak coordinates:
wget https://zenodo.org/records/12793186/files/mouseESC_fig1E_peakorder.bed -P data/

# The union of the O-GlcNac peaks
wget https://zenodo.org/records/12793186/files/union_sept2023mouse_HG.gff -P data/

# The bigwig files of the experiments
wget XXX/ESCHGGlcNAc_rep1.bw  -P data/
wget https://zenodo.org/records/12793186/files/RNApolymeraseII_SRX8556273.bw  -P data/
wget https://zenodo.org/records/12793186/files/Tbp_SRX9195301.bw  -P data/
wget https://zenodo.org/records/12793186/files/Taf12_SRX11221932.bw  -P data/
wget https://zenodo.org/records/12793186/files/Nelfa_SRX017058.bw  -P data/
wget https://zenodo.org/records/12793186/files/Med1_SRX9195310.bw  -P data/
wget https://zenodo.org/records/12793186/files/Med12_SRX1670201.bw  -P data/
wget https://zenodo.org/records/12793186/files/Med24_SRX5926394.bw  -P data/
wget https://zenodo.org/records/12793186/files/Med26_SRX4167136.bw  -P data/
wget https://zenodo.org/records/12793186/files/Dr1_SRX2894853.bw  -P data/

## The peak files of each replicates
wget XXX/ESCHGGlcNAc_rep1_peaks.gff -P data/
wget XXX/ESCHGGlcNAc_rep2_peaks.gff -P data/
```



