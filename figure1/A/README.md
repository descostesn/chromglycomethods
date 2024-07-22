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
- Suz12 peaks:
- RING1B peaks:
- H3K9me3 peaks:
- Ser5P peaks: 
- Ser2P peaks:
- ATACSeq peaks:


- O-GlcNac Cut&Run replicate 1:






peakspathqueryvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc2_lane1sample13_peaks_narrowPeak.gff") # nolint

querynamevec <- c("ESCHGGlcNAc1", "ESCHGGlcNAc2")

repeatsannovec <- c(LINE = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/LINE.gff",  # nolint
        LTR = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/LTR.gff",  # nolint
        SINE = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/SINE.gff") # nolint

outputfolder <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/extendedGenomicCompartments/sept2023_mouse_HG/rep1-2" # nolint

glcnacbwvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_GlcPolIIGlucose_Sept2023/glc_glucoce_mouseESC/data/sequencing/bw/genome_norm/ESCHGGlcNAc1_lane1sample12.bw", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_GlcPolIIGlucose_Sept2023/glc_glucoce_mouseESC/data/sequencing/bw/genome_norm/ESCHGGlcNAc2_lane1sample13.bw") # nolint

includerepeats <- FALSE

outputObjectPath <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/extendedGenomicCompartments/sept2023_mouse_HG/rep1-2/gclist.Rdat" #nolint

countstable <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_mRNASeq_oct2019/recompute_with_ensembl/data/sequencing/counts/featureCounts/featurecounts_MB1-WT-E14-T2i-replicate1.txt.gz.tabular" #nolint
countslength <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_mRNASeq_oct2019/recompute_with_ensembl/data/sequencing/counts/featureCounts/MB1-WT-E14-T2i-replicate1.txt.gz.tabular" # nolint
countsannotype <- "ensembl"

#"https://www.ensembl.org"
biomartstr <- "ENSEMBL_MART_ENSEMBL"
datasetstr <- "mmusculus_gene_ensembl"
hoststr <- "https://nov2020.archive.ensembl.org"
alternativemirroropt <- FALSE


