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
wget XXX/DLD1GlcNAcNoDoxAux_rep1.bw -P data/
wget https://zenodo.org/records/12793186/files/RNApolymeraseII_SRX10580013.bw -P data/

## The gene annotations
wget https://zenodo.org/records/12793186/files/Homo_sapiens.GRCh38.110.chr_march2024_filtered.tar.gz
tar -xvzf Homo_sapiens.GRCh38.110.chr_march2024_filtered.tar.gz

## The peaks coordinates after decreasing sorting of RNAPol II
wget https://zenodo.org/records/12793186/files/peakscoord-fig4B.bed
```

## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig4B.yml](fig4B.yml), run:

```
conda env create -n fig4b --file ./fig4B.yml
conda activate fig4b
```

## Figure generation

Generate a matrix with the RNAPol II bigwig using the coordinates of the Ensembl genes and plot a heatmap in descending order:

```
#!/bin/bash

mkdir results

## Define the number of CPUs
NBCPU=1

## Build the deeptools matrix
computeMatrix scale-regions --regionsFileName Homo_sapiens.GRCh38.110.chr_march2024_filtered.bed --scoreFileName RNApolymeraseII_SRX10580013.bw --outFileName results/polII.mat --samplesLabel RNAPol_II --numberOfProcessors $NBCPU --regionBodyLength 2000 --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --unscaled5prime 0 --unscaled3prime 0

## Plot the RNAPol II signal using decreasing sorting
FILENAME="heatmap_polII.png"

plotHeatmap --matrixFile results/polII.mat --outFileName results/$FILENAME --plotFileFormat 'png' --outFileSortedRegions results/peakscoord-fig4B.bed --dpi '200' --sortRegions 'descend' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --xAxisLabel 'distance from TSS (bp)' --yAxisLabel 'genes' --heatmapWidth 7.5 --heatmapHeight 25.0 --whatToShow 'plot, heatmap and colorbar' --startLabel 'TSS' --endLabel 'TES' --refPointLabel 'TSS' --samplesLabel RNAPolII --legendLocation 'best' --labelRotation '0'
```

Using the sorted peak coordinates `peakscoord-fig4B.bed`, generate a matrix of O-GlcNac signal before and after (Dox)/Auxin treatment:

```
computeMatrix scale-regions --regionsFileName results/peakscoord-fig4B.bed --scoreFileName data/DLD1GlcNAcNoDoxAux_rep1.bw data/DLD1GlcNAcDoxAux_rep1.bw --outFileName results/OGlcNac.mat --samplesLabel GlcNAcNoDox GlcNAcDox --numberOfProcessors $NBCPU --regionBodyLength 2000 --beforeRegionStartLength 2000 --afterRegionStartLength 2000  --unscaled5prime 0 --unscaled3prime 0
```

The matrix is already sorted because it followed the order of `peakscoord-fig4B.bed`. Remains plotting the signal without performing sorting:

```
FILENAME="heatmap_OGlcNac.png"

plotHeatmap --matrixFile results/OGlcNac.mat --outFileName results/$FILENAME  --plotFileFormat 'png' --dpi '200' --sortRegions 'no' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --xAxisLabel 'distance from TSS (bp)' --yAxisLabel 'genes' --heatmapWidth 7.5 --heatmapHeight 25.0 --whatToShow 'plot, heatmap and colorbar' --startLabel 'TSS' --endLabel 'TES' --refPointLabel 'TSS' --samplesLabel NoDox Dox --legendLocation 'best' --labelRotation '0'
```

You should obtain the following heatmaps:

!!!!!!! TO DO
