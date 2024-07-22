# Pan-O-GlcNac Cut&Run in non-repetitive DNA in mouse ES cells

## Description

Upset plot of the location of O-GlcNac Cut&Run peaks (replicates 1/2) in functional genomic compartments: Active promoters (79/84 peaks), Transcription Initiation (64/71 peaks), Heterochromatin (42/52 peaks), Bivalent Promoters (29/30 peaks), PcG (polycombs) Domains (24/23 peaks), Transcription Elongation (23/12 peaks), and Transcription Termination (9/11 peaks). O-GlcNac preferentially occupies active promoters at the transcription initiation sites.

## Details

The genomic compartments were defined as follows:

1) Active promoter: K27ac peaks in TSS-/+1Kb.
2) Transcription initiation: TSS-1/+1kb overlapping with Ser5P peaks.
3) Heterochromatin: H3K9me3 peak intervals.
4) Bivalent promoters: H3K4me3/H3K27me3 peaks overlapping TSS-/+1Kb. 
5) Polycomb domain: Suz12 and RING1B peaks overlapping each other.
6) Transcription elongation: TSS+1kb to TES overlapping Ser2P peaks.
7) Transcription termination: TES+50bp intervals.

## Data

The processed data to use in the script can be downloaded from:

```
#!/bin/bash

mkdir annotations
mkdir data

# gencode, refgene, refseq annotations
wget XXX/gencode.vM25.annotation.gff -P annotations/
wget XXX/refGeneUCSC-mm10-March2021.gff -P annotations/
wget XXX/refseqNCBI-mm10-March2021.gff -P annotations/

# Histone marks
wget XXX/H3K27ac_SRX19148013_peaks_broadPeak.gff -P data/
wget XXX/H3K4me1_SRR5466745_peaks_broadPeak.gff -P data/
wget XXX/H3K27me3_SRR10032683_peaks_hiddendomains.gff -P data/
wget XXX/H3K4me3_SRX5382140_peaks_broadPeak.gff -P data/
wget XXX/H3K9me3_SRR925652_peaks_broadPeak.gff -P data/

# Polycomb marks
wget XXX/Suz12_SRR034190_peaks_broadPeak.gff -P data/
wget XXX/Ring1B_SRR10095137_peaks_narrowPeak.gff -P data/

# CTD marks
wget XXX/Ser5P_SRR391050_peaks_broadPeak.gff -P data/
wget XXX/Ser2P_SRR391039_peaks_broadPeak.gff -P data/

# ATAC-seq
wget XXX/ATAC_SRR5466767_peaks_narrow.gff -P data/

# O-GlcNac peak replicates
wget XXX/ESCHGGlcNAc_rep1.gff -P data/
wget XXX/ESCHGGlcNAc_rep2.gff -P data/

# O-GlcNac bigwig replicates
wget XXX/ESCHGGlcNAc_rep1.bw -P data/
wget XXX/ESCHGGlcNAc_rep2.bw -P data/
```


## Installation

## Figure generation

The script should output the following numbers:

```
# Replicate 1
Performing overlap with each compartment
		 activeProm
		76/5371 compartments contain a peak
		 PcGDomain
		22/10515 compartments contain a peak
		 heteroChrom
		28/48747 compartments contain a peak
		 bivalentProm
		29/13237 compartments contain a peak
		 initiation
		62/3364 compartments contain a peak
		 elongation
		20/2455 compartments contain a peak
		 termination
		9/13962 compartments contain a peak

# Replicate 2
Performing overlap with each compartment
		 activeProm
		82/5371 compartments contain a peak
		 PcGDomain
		21/10515 compartments contain a peak
		 heteroChrom
		37/48747 compartments contain a peak
		 bivalentProm
		29/13237 compartments contain a peak
		 initiation
		70/3364 compartments contain a peak
		 elongation
		11/2455 compartments contain a peak
		 termination
		11/13962 compartments contain a peak
```

You should obtain the raw figure:

[]()

## Pre-processing

### Raw to bam

**Workflows:**
**ChIP-seq and CutnRun:**
**ATAC-seq:**

### Peak detection

**Scripts**: