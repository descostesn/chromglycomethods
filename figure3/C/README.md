# Differential ATAC-seq in mouse Embryonic Stem Cells


I. [Dwget XXX/ESCription](#dwget XXX/ESCription)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; IV.I. [ATAC-seq](#atac-seq)  
&nbsp;&nbsp; IV.II. [Peak Detection](#peak-detection)  
&nbsp;&nbsp; V.IV. [Peak union](#peak-union)  


## DESCription

This experiment uses a transgenic cell line expressing bacterial OGA BtGH84 fused to a localization peptide (NLS) and regulated by Tet-On system. OGA is a glycosidase that removes O-GlcNAc modifications. We evaluated the changes in chromatin accessibility before and after O-GlcNac removal by OGA. No differences were found which tells us that O-GlcNac does not influence the chromatin accessibility.

## Data


```
#!/bin/bash

mkdir data

## The bam files of ATAC-seq before and after doxocyclin induction
wget XXX/ESC_1b_atac_NoDox_rep1.bam -P data/
wget XXX/ESC_1b_atac_NoDox_rep2.bam -P data/
wget XXX/ESC_1b_atac_NoDox_rep3.bam -P data/
wget XXX/ESC_1b_atac_Dox_rep1.bam -P data/
wget XXX/ESC_1b_atac_Dox_rep2.bam -P data/
wget XXX/ESC_1b_atac_Dox_rep3.bam -P data/

## The peak annotations (not necessary to generate the MA plot, see union file hereafter)
wget XXX/ESC_1b_atac_NoDox_rep1_peaks.bed -P data/
wget XXX/ESC_1b_atac_NoDox_rep2_peaks.bed -P data/
wget XXX/ESC_1b_atac_Dox_rep3_peaks.bed -P data/
wget XXX/ESC_1b_atac_Dox_rep1_peaks.bed -P data/
wget XXX/ESC_1b_atac_Dox_rep2_peaks.bed -P data/
wget XXX/ESC_1b_atac_NoDox_rep3_peaks.bed -P data/

## Annotations of the union of ATAC-seq peaks
wget https://zenodo.org/api/records/12793186/files/ESC1bNoDox_vs_Dox.gtf -P data/
```


## Figure Generation

The MA plot was obtained with DESeq2 1.22.1. The workflow [OGlcNac_deseq2_atacseqPE_peaks](galaxy-workflow/Galaxy-Workflow-OGlcNac_deseq2_atacseqPE_peaks.ga). Two lists of bam files for each condition (NoDox/Dox) are given as input. The counts are retrieved on the union of peaks (ESC1bNoDox_vs_Dox.gtf) with the featureCounts function of subread v2.0.3: `featureCounts -a ESC1bNoDox_vs_Dox.gtf -F "GTF" -o $output -T $nbcpu -s 0 -Q 0 -t 'peak' -g 'peak_id' --minOverlap 1 --fracOverlap 0 --fracOverlapFeature 0 -p -B -C --countReadPairs $input.bam`. The lists of counts are used as input of DESeq2 1.22.1 that is launched with the script [deseq2.R](../../figure2/B/others/deseq2.R). The deseq2 table was then filtered on columns 3 (fold-change) and 7 (adj p-val) with `abs(c3)>0 and c7<0.05` to keep significant fold-changes that are not equal to zero. The lists of down- and up-regulated genes were separated to their respective files using the third column: `c3<0` for down and `c3>0` for up.

## Pre-processing

### ATAC-seq
### Peak Detection
### Peak union
