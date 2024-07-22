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

- gencode annotation: gencode.vM25.annotation.gff
- refgene annotation: refGeneUCSC-mm10-March2021.gff
- refseq annotation: refseqNCBI-mm10-March2021.gff
- H3K27ac peaks: H3K27ac_SRX19148013_peaks_broadPeak.gff
- H3K4me1 peaks: H3K4me1_SRR5466745_peaks_broadPeak.gff
- H3K27me3 peaks: H3K27me3_SRR10032683_peaks_hiddendomains.gff
- H3K4me3 peaks: H3K4me3_SRX5382140_peaks_broadPeak.gff
- Suz12 peaks: Suz12_SRR034190_peaks_broadPeak.gff
- RING1B peaks: Ring1B_SRR10095137_peaks_narrowPeak.gff
- H3K9me3 peaks: H3K9me3_SRR925652_peaks_broadPeak.gff
- Ser5P peaks: Ser5P_SRR391050_peaks_broadPeak.gff
- Ser2P peaks: Ser2P_SRR391039_peaks_broadPeak.gff
- ATACSeq peaks: ATAC_SRR5466767_peaks_narrow.gff
- O-GlcNac replicate 1 peaks: ESCHGGlcNAc_rep1.gff
- O-GlcNac replicate 2 peaks: ESCHGGlcNAc_rep2.gff

## Installation

## Figure generation

## Pre-processing

### Raw to bam

**Workflows:**
**ChIP-seq and CutnRun:**
**ATAC-seq:**

### Peak detection

**Scripts**: