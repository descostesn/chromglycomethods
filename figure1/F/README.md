# Heatmap of O-GlcNac peaks along with ChIP-Atlas candidates

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
&nbsp;&nbsp; IV.II. [Union of peaks](#union-of-peaks)  
&nbsp;&nbsp; IV.II. [RNAPol II clustering](#rnapol-ii-clustering)  
&nbsp;&nbsp; IV.III. [Heatmaps](#heatmaps)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Data](#data-1)  
&nbsp;&nbsp; V.II. [Workflow](#workflow)  
&nbsp;&nbsp; V.I. [Peak detection](#peak-detection)  


## Description

Heatmap of CuntRun O-GlcNac signal and ChIP-Seq ChIP-Atlas selected candidates signal 1 kb around the start of 702 O-GlcNac peaks. The list of peaks was obtained computing the union of the replicate peaks. Three groups were defined and highlight the RNA Polymerase II (RNAPol II) before, at, and after the O-GlcNac peaks. For the groups I and II, RNAPol II describes a double peak signal overlapping O-GlcNac independently of its position. Other selected candidates follow a similar pattern. It indicates that a sub-category of these proteins might carry O-GlcNac in a transient state. The description of the function of each protein can be found in [fig1D](../D/README.md).

## Data

Download the following data:

```
#!/bin/bash

mkdir data

# The bed file containing the sorted peak coordinates (file is correct, panel letter changed):
wget https://zenodo.org/records/12793186/files/peakscoord-fig1E.bed -P data/

# The union of the O-GlcNac peaks
wget https://zenodo.org/records/12793186/files/union_sept2023mouse_HG1-2.gff -P data/

# The bigwig files of the experiments
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep1.bw  -P data/
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
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep1_peaks.gff -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep2_peaks.gff -P data/
```

## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1F.yml](fig1F.yml), run:

```
conda env create -n fig1f --file ./fig1f.yml
conda activate fig1f
```

## Figure Generation

### Union of peaks

Generate the bed file of the union of peaks by running:

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

### RNAPol II clustering

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

### Heatmaps


Performing the same procedure for the other datasets:

```
FILEVEC=( 'ESCHGGlcNAc_rep1.bw' 'Tbp_SRX9195301.bw' 'Taf12_SRX11221932.bw' 'Nelfa_SRX017058.bw' 'Med1_SRX9195310.bw' 'Med12_SRX1670201.bw' 'Med24_SRX5926394.bw' 'Med26_SRX4167136.bw' 'Dr1_SRX2894853.bw' )
HEATMAPMAXVEC=( 25 80 150 340 50 28 38 85 170 )
PROFILEMAXVEC=( 15 25 40 105 20 12 15 25 55 )
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

The script should output:

```
##  ESCHGGlcNAc_rep1.bw
-- matrix
-- heatmap
##  Tbp_SRX9195301.bw
-- matrix
-- heatmap
##  Taf12_SRX11221932.bw
-- matrix
-- heatmap
##  Nelfa_SRX017058.bw
-- matrix
-- heatmap
##  Med1_SRX9195310.bw
-- matrix
-- heatmap
##  Med12_SRX1670201.bw
-- matrix
-- heatmap
##  Med24_SRX5926394.bw
-- matrix
-- heatmap
##  Med26_SRX4167136.bw
-- matrix
-- heatmap
##  Dr1_SRX2894853.bw
-- matrix
-- heatmap
```


## Pre-processing

### Data

| Target | ID | library layout | link |
|--------|----|----------------|------|
| O-GlcNAc rep1 | E-MTAB-14308 | single | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/033/ERR13430733/ERR13430733.fastq.gz |
| RNApolymeraseII | SRX8556273 | single | see the provided [snakemake](snakemake/Snakefile) (technical replicates) | 
| Tbp | SRR12716556 | single | see the provided [snakemake](snakemake/Snakefile) |
| Taf12 | SRR14907447 | single | see the provided [snakemake](snakemake/Snakefile) |
| Nelfa | SRR036737 | single | see the provided [snakemake](snakemake/Snakefile) |
| Med1 | SRR12716565 | single | see the provided [snakemake](snakemake/Snakefile) |
| Med12 | SRR3313261 | single | see the provided [snakemake](snakemake/Snakefile) |
| Med24 | SRR9153306 | single | see the provided [snakemake](snakemake/Snakefile) |
| Med26 | SRR7262971 | single | see the provided [snakemake](snakemake/Snakefile) |
| Dr1 | SRR5658441 | single | see the provided [snakemake](snakemake/Snakefile) |


### Workflow

The pre-processing was performed with the Galaxy workflows [OGlcNac_ChIP-SeqPEmm10](../A/galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqPEmm10.ga) and [OGlcNac_ChIP-SeqSEmm10](../A/galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqSEmm10.ga) for paired-end and single-end data respectively. The .ga files can be imported in one own galaxy account.

Quality control was done with FastQC v0.11.9: `fastqc --outdir $outputfolder --threads $nbcpu --quiet --extract --kmers 7 -f 'fastq' $input.fastq.gz`.

Adapters and low quality reads were removed with trim-galore v0.4.3: `trim_galore --phred33 --quality 20 --stringency 1 -e 0.1 --length 20 --output_dir ./ $input.fastq.gz`.

Reads were aligned to mm10 with Bowtie 2.3.4.1 and the bam were sorted using samtools v1.9:
single: `bowtie2 -p $nbcpu -x m.musculus/mm10/mm10 -U $input.fastq.gz --sensitive --no-unal 2> $log |  samtools sort -@$nbcpu -O bam -o $output.bam`
paired: `bowtie2 -p $nbcpu -x m.musculus/mm10/mm10 -1 $input1.fastq.gz -2 $input2.fastq.gz -I 0 -X 500 --fr --dovetail --sensitive --no-unal 2> $log  | samtools sort -@$nbcpu -O bam -o $output.bam`.

Only primary alignments were kept using samtools v1.9: `samtools view -o $output.bam -h -b -q 20 -F 0x800 $input.bam`.

Reads not aligned to consensus chromosomes were excluded with samtools v1.9: `samtools view -o $output.bam -h -b $input.bam 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chrX' 'chrY'`.

Duplicates were removed with picard v2.18.2: `picard MarkDuplicates INPUT=$input.bam OUTPUT=$output.bam METRICS_FILE=$metrics.txt REMOVE_DUPLICATES='true' ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR`

Bigwig files normalized by the genome size were generated with deeptools v3.0.2:
single: `bamCoverage --numberOfProcessors $NBCPU --bam $input.bam --outFileName $output.bw --outFileFormat 'bigwig' --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 --scaleFactor 1.0  --extendReads 150 --minMappingQuality '1'`
paired: `bamCoverage --numberOfProcessors $NBCPU --bam $input.bam --outFileName $output.bw --outFileFormat 'bigwig' --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 --scaleFactor 1.0  --extendReads --minMappingQuality '1'`

### Peak detection

| Target | Broad | q-value | Duplicates Thres. | Tag size |
|--------|-------|---------|-------------------|----------|
| ESCHGGlcNAc_rep1 | NO | 0.04 | 7 | 82 |
| ESCHGGlcNAc_rep2 | NO | 0.04 | 7 | 82 |

* Macs2 v2.2.7.1 Narrow: `macs2 callpeak -t $input.bam -c $control.bam -n $expname --outdir $outfold -f BAM -g 1.87e9 -s $tagsize -q $qvalue --nomodel --extsize 150 --keep-dup $dupthres`
