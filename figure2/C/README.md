# Genome Browser Screenshot

I. [Description](#description)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  


## Description

Genome browser view of O-GlcNac and RNA-seq signal in control and siOgt Knock-Down (KD) conditions. Both genes show downregulation following Ogt KD.

## Data

The bigwig files used in this figure can be obtained with:

```
#!/bin/bash

## The screenshots were performed with rep 2
wget https://zenodo.org/records/13444099/files/TPMbw_sictrl_rep2.bw
wget https://zenodo.org/records/13444099/files/TPMbw_siogt_rep2.bw

## Replicate 1 is also provided
wget https://zenodo.org/records/13444099/files/TPMbw_sictrl_rep1.bw
wget https://zenodo.org/records/13444099/files/TPMbw_siogt_rep1.bw
```

## Figure Generation

The bw files were uploaded to [IGV](https://igv.org/) v2.13.0 selecting the mouse mm10 genome.


## Pre-processing

The pre-processing was performed with the Galaxy workflows [Galaxy-Workflow-OGlcNac_RNASeqPE_mm10_STAR_bw](../B/galaxy-workflows/Galaxy-Workflow-OGlcNac_RNASeqPE_mm10_STAR_bw.ga). The .ga file can be imported in your own galaxy account.

The file to compute the count tables can be downloaded from `wget https://zenodo.org/records/13444099/files/Mus_musculus.GRCm38.102.chr.gtf.tar.gz`

FastQC 0.11.9 was used for quality control: `fastqc --outdir $outfolder --threads $nbcpu --quiet --extract --kmers 7 -f 'fastq' input.fastq.gz`.

Adapters and low quality reads were removed with trim-galore 0.4.3: `trim_galore --phred33 --quality 20  --stringency 1 -e 0.1 --length 20 --output_dir ./ $input.fastq.gz`.

Alignment was performed with STAR 2.6.0b: `STAR --runThreadN $nbcpu --genomeLoad NoSharedMemory --genomeDir 'mm10/rnastar_index2/mm10/files' --readFilesIn $input.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --outSAMstrandField None --outFilterIntronMotifs RemoveNoncanonical --outFilterIntronStrands RemoveInconsistentStrands --outSAMunmapped None --outSAMprimaryFlag OneBestScore --outSAMmapqUnique "255" --outFilterType Normal --outFilterMultimapScoreRange "1" --outFilterMultimapNmax "10" --outFilterMismatchNmax "10" --outFilterMismatchNoverLmax "0.3" --outFilterMismatchNoverReadLmax "1.0" --outFilterScoreMin "0" --outFilterScoreMinOverLread "0.66" --outFilterMatchNmin "0" --outFilterMatchNminOverLread "0.66" --outSAMmultNmax "-1" --outSAMtlen "1" --outBAMsortingBinsN "50"`.

Counts were obtained with subread v2.0.1: `featureCounts -a Mus_musculus.GRCm38.102.chr.gtf -F GTF -o ESCRNAseq_SRR11294181counts.txt -T $nbcpu -s 0 -Q 0 -t 'exon' -g 'gene_id' --minOverlap 1 --fracOverlap 0 --fracOverlapFeature 0 -C input.bam`.

The bigwig files were generated with deeptools v 3.0.2: `bamCoverage --numberOfProcessors $nbcpu --bam one.bam --outFileName $output.bigwig --outFileFormat 'bigwig' --binSize 50 --normalizeUsing BPM`
