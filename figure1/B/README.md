# Genome Browser Screenshot

I. [Description](#description)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  


## Description

Genome browser view of O-GlcNac, H3K4me3, and H3K27ac at the Zfp639 and Mfn1 genes. The peaks of H3K4me3 and H3K27ac indicates that both genes are actively transcribed. However, O-GlcNac peak can be seen at Zfp639 but not at Mfn1. This example indicates that O-GlcNac is seen at a subset of active promoters.

## Data

The bigwig files used in this figure can be obtained with:

```
#!/bin/bash

wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep1.bw
wget https://zenodo.org/records/12793186/files/H3K27ac_SRX19148013.bw
wget https://zenodo.org/records/12793186/files/H3K4me3_SRX5382140.bw
```

## Figure Generation

The bw files were uploaded to [IGV](https://igv.org/) v2.13.0 selecting the mouse mm10 genome.


## Pre-processing

The H3K4me3 and H3K27ac data was processed with [OGlcNac_ChIP-SeqPEmm10](../A/galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqPEmm10.ga) and ESCHGGlcNAc with [OGlcNac_ChIP-SeqSEmm10](../A/galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqSEmm10.ga).

Quality control was done with FastQC v0.11.9: `fastqc --outdir $outputfolder --threads $nbcpu --quiet --extract --kmers 7 -f 'fastq' $input.fastq.gz`.

Adapters and low quality reads were removed with trim-galore v0.4.3: `trim_galore --phred33 --quality 20 --stringency 1 -e 0.1 --length 20 --output_dir ./ $input.fastq.gz`.

Reads were aligned to mm10 with Bowtie 2.3.4.1 and the bam were sorted using samtools v1.9:
single: `bowtie2 -p $nbcpu -x m.musculus/mm10/mm10 -U $input.fastq.gz --sensitive --no-unal 2> $log |  samtools sort -@$nbcpu -O bam -o $output.bam`
paired: `bowtie2 -p $nbcpu -x m.musculus/mm10/mm10 -1 $input1.fastq.gz -2 $input2.fastq.gz -I 0 -X 500 --fr --dovetail --sensitive --no-unal 2> $log  | samtools sort -@$nbcpu -O bam -o $output.bam`.

Only primary alignments were kept using samtools v1.9: `samtools view -o $output.bam -h -b -q 20 -F 0x800 $input.bam`.

Reads not aligned to consensus chromosomes were excluded with samtools v1.9: `samtools view -o $output.bam -h -b $input.bam 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chrX' 'chrY'`.

The bigwig files were generated normalizing by the effective genome size using Deeptools 3.0.2: `bamCoverage --numberOfProcessors $nbcpu  --bam $input.bam --outFileName $output.bw --outFileFormat 'bigwig' --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 --scaleFactor 1.0 --extendReads 150 --minMappingQuality 1`
