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

## Ensembl annotations
wget https://github.com/descostesn/chromglycomethods/blob/main/figure1/E/annotations/ensemblmm10.tar.gz -P data/
tar -xvzf data/ensemblmm10.tar.gz
```

## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1E.yml](fig1E.yml), run:

```
conda env create -n fig1E --file ./fig1E.yml
conda activate fig1E
```

## Figure generation

Generate a matrix with the RNAPol II bigwig using the coordinates of the Refseq genes and plot a heatmap in descending order (deeptools v3.5.5):

```
#!/bin/bash

mkdir results

## Define the number of CPUs
NBCPU=1

## Build the deeptools matrix
computeMatrix scale-regions --regionsFileName data/ensembl_mm10_Sept2020-chrfiltered.gff --scoreFileName data/RNApolymeraseII_SRX8556273.bw --outFileName results/polII.mat --samplesLabel RNAPol_II --numberOfProcessors $NBCPU --regionBodyLength 2000 --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --unscaled5prime 0 --unscaled3prime 0

## Plot the RNAPol II signal using decreasing sorting
FILENAME="heatmap_polII.png"

plotHeatmap --matrixFile results/polII.mat --outFileName results/$FILENAME --plotFileFormat 'png' --outFileSortedRegions results/peakscoord-fig1E.bed --dpi '200' --sortRegions 'descend' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --xAxisLabel 'distance from TSS (bp)' --yAxisLabel 'genes' --heatmapWidth 7.5 --heatmapHeight 25.0 --whatToShow 'plot, heatmap and colorbar' --startLabel 'TSS' --endLabel 'TES' --refPointLabel 'TSS' --samplesLabel RNAPolII --legendLocation 'best' --labelRotation '0'
```

You should obtain the raw figure:

<img src="pictures/heatmap_polII.png" alt="PolII heatmap" width="200" heigth="200"/>

Using the sorted peak coordinates `peakscoord-fig1E.bed`, generate a matrix of O-GlcNac signal:

```
computeMatrix scale-regions --regionsFileName results/peakscoord-fig1E.bed --scoreFileName data/ESCHGGlcNAc_rep1.bw --outFileName results/OGlcNac.mat --samplesLabel Glc --numberOfProcessors $NBCPU --regionBodyLength 2000 --beforeRegionStartLength 2000 --afterRegionStartLength 2000  --unscaled5prime 0 --unscaled3prime 0
```

The matrix is already sorted because it followed the order of `peakscoord-fig1E.bed`. Remains plotting the signal without performing sorting:

```
FILENAME="heatmap_OGlcNac.png"

plotHeatmap --matrixFile results/OGlcNac.mat --outFileName results/$FILENAME  --plotFileFormat 'png' --dpi '200' --sortRegions 'no' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --xAxisLabel 'distance from TSS (bp)' --yAxisLabel 'genes' --heatmapWidth 7.5 --heatmapHeight 25.0 --whatToShow 'plot, heatmap and colorbar' --startLabel 'TSS' --endLabel 'TES' --refPointLabel 'TSS' --samplesLabel NoDox Dox --legendLocation 'best' --labelRotation '0'
```

You should obtain the following heatmaps:

<img src="pictures/heatmap_OGlcNac.png" alt="OGlcNac heatmap" width="400" heigth="400"/>


