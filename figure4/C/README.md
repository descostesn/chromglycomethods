# Clustering of O-GlcNac before and after RNA Polymerase II removal reveals different functional categories

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [ChIP-seq and CutnRun](#chip-seq-and-cutnrun)  
&nbsp;&nbsp; V.II. [Peak Detection](#peak-detection)  
&nbsp;&nbsp; V.I. [Peak Union](#peak-union)  




## Description

The comparison of O-GlcNac signal at 6,544 loci (union of replicate peaks) revealed 3 categories of loci (5 groups) in which the reponse of O-GlcNac occupancy to the removal of RNA Polymerase II (RNAPol II, see [fig4B](../B/README.md#description)) is different. Cluster 1-2-3 (284, 584, and 1490 peaks respectively) do not see major changes in O-GlcNac occupancy. The observed O-GlcNac signal is likely coming from other proteins than RNAPol II. Similarly, cluster 4 (1,200 peaks) is RNAPol II independent as poor signal is observed prior removal. However, the absence of RNAPol II enables the binding of proteins carrying O-GlcNac as an increase of signal is observed. It is tempting to speculate that this cluster has an actively transcribing RNAPol II, the absence of signal being explained by the lack of time for the antibody to catch O-GlcNac. Upon RNAPol II removal, the binding of repressing factors should decrease the production of RNAs. Strikingly, the figS3 shows a decrease in nascent rna signal at these loci. Cluster 5 (2,986 peaks) is the RNAPol II dependent one since it sees a decrease of O-GlcNac signal upon its removal. Overall, we can envision O-GlcNac as being a key component of transcription whether it is active or inactive. Perturbation of RNAPol II activity reveals a re-localization of O-GlcNac occupancy at thousands of loci.


## Data

```
#!/bin/bash

mkdir data

## The bigwigs of O-GlcNca before and after treatment by (Dox)/Auxin
wget XXX/DLD1GlcNAcNoDoxAux_rep1.bw -P data/
wget XXX/DLD1GlcNAcDoxAux_rep1.bw -P data/

## The bam files to perform the peak detection
wget XXX/DLD1GlcNAcDoxAux_rep1.bam -P data/
wget XXX/DLD1GlcNAcDoxAux_rep2.bam -P data/
wget XXX/DLD1GlcNAcNoDoxAux_rep1.bam -P data/
wget XXX/DLD1GlcNAcNoDoxAux_rep2.bam -P data/

## The peaks to perform the union
wget XXX/DLD1GlcNAcNoDoxAux_rep1.gff -P data/
wget XXX/DLD1GlcNAcNoDoxAux_rep2.gff -P data/
wget XXX/DLD1GlcNAcDoxAux_rep1_peaks.gff -P data/
wget XXX/DLD1GlcNAcDoxAux_rep2_peaks.gff -P data/

## The bigwigs of RNAPol II used as control
wget https://zenodo.org/records/12793186/files/RNAPolII_SRX11070611_control.bw -P data/
wget https://zenodo.org/records/12793186/files/RNAPolII_SRX11070613_auxin.bw -P data/

## The file of the union of the peaks
wget https://zenodo.org/records/12793186/files/union_OGlcNac_noauxaux-fig4C.bed -P data

## The coordinates of the peaks sorted in 5 groups
wget https://zenodo.org/records/12793186/files/peakscoord-fig4C.bed
```


## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig4C.yml](fig4C.yml), run:

```
conda env create -n fig4c --file ./fig4C.yml
conda activate fig4c
```

## Figure Generation

Using the union of peaks as the reference loci, compute a matrix of O-GlcNac without and with Auxin treatment. This matrix is used to plot a heatmap with K-means clustering with 5 groups (deeptools v3.5.5):

```
#!/bin/bash

mkdir results

NBCPU=1

computeMatrix reference-point --regionsFileName data/union_OGlcNac_noauxaux-fig4C.bed --scoreFileName DLD1GlcNAcNoDoxAux_rep1.bw DLD1GlcNAcDoxAux_rep1.bw --outFileName results/OGlcNacnoauxaux.mat --samplesLabel DLD1GlcNAcNoAux DLD1GlcNAcAux  --numberOfProcessors $NBCPU --referencePoint TSS  --beforeRegionStartLength 1000 --afterRegionStartLength 1000
```






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


To retrieve the genes associated to promoters in each cluster of the heatmap, the script src\R\tools\vennDiagram_twoExp.R has been updated. The conf for overlap is src\R\conf_generation\venndiagrams_twoexp\glcpolIIhumansept2023\humanpolIIglcnac_heatmapclusters.R

The results are in /g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_*-compartmentsgff/promoters_vs_ensemblgenes/

‌The files to perform the union of peaks are: /home/descostes/vscode/o-n-acetylglucosamine/src/R/conf_generation/make_union/sept2023_glcPolII_human.R

For the heatmaps, the galaxy histories are:
- human glc polII heatmaps

The matrices, heatmaps, and coordinates are in:
- Human: /g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/5groups

With the division in 5 groups, of 6544 peaks, the number of peaks per cluster is:

cluster_1: 284
cluster_2: 584
cluster_3: 1490
cluster_4: 1200
cluster_5 : 2986

The script to compute the genomic repartition of each cluster of the human heatmap is src/R/tools/detailsGroupsHeatmapHumanSept2023.R and the results are written in the heatmap folder (/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!