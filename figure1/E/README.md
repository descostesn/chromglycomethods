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
wget https://zenodo.org/records/12793186/files/peakscoord-fig1E.bed -P data/

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
#!/bin/bash

## Define the number of CPUs
NBCPU=1

## Build the deeptools matrix
computeMatrix  reference-point --regionsFileName results/union_sept2023mouse_HG1-2.bed --scoreFileName data/RNApolymeraseII_SRX8556273.bw --outFileName results/polIImatrix.mat --samplesLabel RNAPolII --numberOfProcessors $NBCPU --referencePoint TSS --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean' --binSize 50 --transcriptID transcript --exonID exon --transcript_id_designator transcript_id
```

Generate a heatmap of RNAPol II performing a k-means clustering on 3 groups to retrieve the peak coordinates:

```
#!/bin/bash

FORMAT='png'
FILENAME='rnapolII-3groups-initial.png'
plotHeatmap --matrixFile results/polIImatrix.mat --outFileName results/$FILENAME --plotFileFormat $FORMAT --outFileSortedRegions results/peakscoord-fig1E.bed --dpi '200' --sortRegions 'descend'  --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --xAxisLabel 'distance from peak start (bp)' --yAxisLabel 'peaks' --heatmapWidth 7.5 --heatmapHeight 25.0 --whatToShow 'plot, heatmap and colorbar' --startLabel 'peak' --endLabel 'TES' --refPointLabel 'peak' --legendLocation 'best' --labelRotation '0' --kmeans 3
```

Replace the groups 'cluster_2' and 'cluster_3' in results/peakscoord-fig1E.bed to avoid visual separation of the groups:

```
#!/bin/bash

sed "s/cluster_[2-3]/cluster_1/" results/peakscoord-fig1E.bed > results/tmp.bed
rm results/peakscoord-fig1E.bed
mv results/tmp.bed results/peakscoord-fig1E.bed
```

Regenerate the RNAPol II matrix with the new coordinates and plot the signal without changing the coordinate order:

```
#!/bin/bash

## Matrix
computeMatrix  reference-point --regionsFileName results/peakscoord-fig1E.bed --scoreFileName data/RNApolymeraseII_SRX8556273.bw --outFileName results/polIImatrixsorted.mat --samplesLabel RNAPolII --numberOfProcessors $NBCPU --referencePoint TSS --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean' --binSize 50 --transcriptID transcript --exonID exon --transcript_id_designator transcript_id

## Heatmap
FORMAT='png'
FILENAME='rnapolII-3groups.png'
plotHeatmap --matrixFile results/polIImatrixsorted.mat --outFileName results/$FILENAME --plotFileFormat $FORMAT    --dpi '200' --sortRegions 'no' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines'  --missingDataColor 'black' --alpha '1.0' --colorList white,blue --zMin 0 --zMax 180 --yMin 0.0 --yMax 80.0  --xAxisLabel 'distance from peak start (bp)' --yAxisLabel 'peaks' --heatmapWidth 7.5 --heatmapHeight 25.0  --whatToShow 'plot, heatmap and colorbar' --startLabel 'peak' --endLabel 'TES' --refPointLabel 'peak'     --legendLocation 'best' --labelRotation '0'
```

Performing the same procedure for the other datasets:

```
FILEVEC=( 'ESCHGGlcNAc_rep1.bw' 'Tbp_SRX9195301.bw' 'Taf12_SRX11221932.bw' 'Nelfa_SRX017058.bw' 'Med1_SRX9195310.bw' 'Med12_SRX1670201.bw' 'Med24_SRX5926394.bw' 'Med26_SRX4167136.bw' 'Dr1_SRX2894853.bw' )
HEATMAPMAXVEC=( 100 70 80 250 50 40 40 60 140 )
PROFILEMAXVEC=( 60 30 40 100 20 18 18 30 50 )
FORMAT='png'

for i in "${!FILEVEC[@]}"
do
    filename=${FILEVEC[$i]}
    echo "## " $filename
    rootname="$(basename $filename .bw)"
    
    echo "-- matrix"
    computeMatrix  reference-point --regionsFileName results/peakscoord-fig1E.bed --scoreFileName data/$filename --outFileName results/$rootname-matrixsorted.mat --samplesLabel $rootname --numberOfProcessors $NBCPU --referencePoint TSS --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --sortRegions 'keep' --sortUsing 'mean' --averageTypeBins 'mean' --binSize 50 --transcriptID transcript --exonID exon --transcript_id_designator transcript_id

    echo "-- heatmap"
    zmaxi=${HEATMAPMAXVEC[$i]}
    ymaxi=${PROFILEMAXVEC[$i]}
    plotHeatmap --matrixFile results/$rootname-matrixsorted.mat --outFileName results/$rootname.$FORMAT --plotFileFormat $FORMAT --dpi '200' --sortRegions 'no' --sortUsing 'mean' --averageTypeSummaryPlot 'mean' --plotType 'lines' --missingDataColor 'black' --alpha '1.0' --colorList white,blue --zMin 0 --zMax $zmaxi --yMin 0.0 --yMax $ymaxi --xAxisLabel 'distance from peak start (bp)' --yAxisLabel 'peaks' --heatmapWidth 7.5 --heatmapHeight 25.0 --whatToShow 'plot, heatmap and colorbar' --startLabel 'peak' --endLabel 'TES' --refPointLabel 'peak' --legendLocation 'best' --labelRotation '0'

    rm results/$rootname-matrixsorted.mat
done
```

