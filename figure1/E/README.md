# Heatmap of O-GlcNac peaks along with ChIP-Atlas candidates

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
&nbsp;&nbsp; IV.II. [Union of peaks](#union-of-peaks)  
&nbsp;&nbsp; IV.II. [Matrix Generation](#matrix-generation)
&nbsp;&nbsp; IV.III. [Heatmaps](#heatmaps)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [ChIP-seq and CutnRun](#cutnrun)  



!! Define a single heatmap with 3 groups for pol II
!! Delete mouseESC_fig1E_peakorder.bed on zenodo and upload the new file with same name

## Description

Heatmap of CuntRun O-GlcNac signal and ChIP-Seq ChIP-Atlas selected candidates signal 1 kb around the start of 702 O-GlcNac peaks. The list of peaks was obtained computing the union of the replicate peaks. Three groups were defined and highlight the RNA Polymerase II (RNAPol II) before, at, and after the O-GlcNac peaks. For the groups I and II, RNAPol II describes a double peak signal overlapping O-GlcNac independently of its position. Other selected candidates follow a similar pattern. It indicates that a sub-category of these proteins might carry O-GlcNac in a transient state. The description of the function of each protein can be found in [fig1D](../D/README.md).

## Data

Download the following data:

```
#!/bin/bash

mkdir data

# The bed file containing the sorted peak coordinates:
wget https://zenodo.org/records/12793186/files/mouseESC_fig1E_peakorder.bed -P data/

# The union of the O-GlcNac peaks
wget https://zenodo.org/records/12793186/files/union_sept2023mouse_HG1-2.gff -P data/

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

## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1E.yml](fig1E.yml), run:

```
conda env create -n fig1e --file ./fig1E.yml
conda activate fig1e
```

## Figure Generation

Generate the gff file of the union of peaks by running:

```
Rscript union.R
```

The script should give the output:

```
Reading peak files
Reducing intervals
The union returned 702 peaks
Writing results/union_sept2023mouse_HG1-2.bed
```

Generate a matrix with the RNAPolII bigwig using the coordinates of the union of peaks:

```
## Define the number of CPUs
NBCPU=1

## Build the deeptools matrix
computeMatrix  reference-point --regionsFileName results/union_sept2023mouse_HG1-2.bed --scoreFileName data/RNApolymeraseII_SRX8556273.bw --outFileName results/polIImatrix.mat --samplesLabel RNAPolII --numberOfProcessors $NBCPU --referencePoint TSS --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean'   --binSize 50  --transcriptID transcript --exonID exon --transcript_id_designator transcript_id
```

Generate a heatmap of RNAPol II performing a k-means clustering on 3 groups to retrieve the peak coordinates:

```
plotHeatmap --matrixFile results/polIImatrix.mat --outFileName results/rnapolII-3groups-initial.pdf --plotFileFormat 'pdf'   --outFileSortedRegions results/peakscoord.bed --dpi '200' --sortRegions 'descend' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --xAxisLabel 'distance from peak start (bp)' --yAxisLabel 'peaks' --heatmapWidth 7.5 --heatmapHeight 25.0  --whatToShow 'plot, heatmap and colorbar' --startLabel 'peak' --endLabel 'TES' --refPointLabel 'peak'     --legendLocation 'best'  --labelRotation '0'
```


