# Genome Browser view of two downregulated genes following siOgt transfection

I. [Description](#description)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; IV.I. [Data](#data-1)  
&nbsp;&nbsp; IV.II. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; IV.II.I. [CutnRun](#cutnrun)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; IV.II.II. [RNA-seq](#rna-seq)  


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

### Data

| Target | ID | library layout | link 1| link 2 |
|--------|----|----------------|-------|--------|
| O-GlcNAc rep1 | E-MTAB-14308 | single | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/033/ERR13430733/ERR13430733.fastq.gz ||
| sictrl rep 2 | E-MTAB-14313 | paired | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/071/ERR13430771/ERR13430771_1.fastq.gz | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/071/ERR13430771/ERR13430771_2.fastq.gz |
| siogt rep 2 | E-MTAB-14313 | paired | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/073/ERR13430773/ERR13430773_1.fastq.gz | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/073/ERR13430773/ERR13430773_2.fastq.gz |

To download the data run:

```
#!/bin/bash

mkdir data/fastq

## Fastq files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/033/ERR13430733/ERR13430733.fastq.gz -P data/ESCHGGlcNAc_rep1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/071/ERR13430771/ERR13430771_1.fastq.gz -P data/sictrl_rep2_fwd.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/071/ERR13430771/ERR13430771_2.fastq.gz -P data/sictrl_rep2_rev.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/073/ERR13430773/ERR13430773_1.fastq.gz -P data/siogt_rep2_fwd.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/073/ERR13430773/ERR13430773_2.fastq.gz -P data/siogt_rep2_rev.fastq.gz
```


### Workflows

#### CutnRun

The pre-processing was performed with the Galaxy workflows [OGlcNac_ChIP-SeqSEmm10](../../figure1/A/galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqSEmm10.ga). The .ga file can be imported in your own galaxy account.

Quality control was done with FastQC v0.11.9: `fastqc --outdir $outputfolder --threads $nbcpu --quiet --extract --kmers 7 -f 'fastq' $input.fastq.gz`.

Adapters and low quality reads were removed with trim-galore v0.4.3: `trim_galore --phred33 --quality 20 --stringency 1 -e 0.1 --length 20 --output_dir ./ $input.fastq.gz`.

Reads were aligned to mm10 with Bowtie 2.3.4.1 and the bam were sorted using samtools v1.9: `bowtie2 -p $nbcpu -x m.musculus/mm10/mm10 -U $input.fastq.gz --sensitive --no-unal 2> $log |  samtools sort -@$nbcpu -O bam -o $output.bam`

Only primary alignments were kept using samtools v1.9: `samtools view -o $output.bam -h -b -q 20 -F 0x800 $input.bam`.

Reads not aligned to consensus chromosomes were excluded with samtools v1.9: `samtools view -o $output.bam -h -b $input.bam 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chrX' 'chrY'`.

Duplicates were removed with picard v2.18.2: `picard MarkDuplicates INPUT=$input.bam OUTPUT=$output.bam METRICS_FILE=$metrics.txt REMOVE_DUPLICATES='true' ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR`

Bigwig files normalized by the genome size were generated with deeptools v3.0.2:
`bamCoverage --numberOfProcessors $NBCPU --bam $input.bam --outFileName $output.bw --outFileFormat 'bigwig' --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 --scaleFactor 1.0  --extendReads 150 --minMappingQuality '1'`

#### RNA-seq

The pre-processing was performed with the Galaxy workflows [OGlcNac_RNASeqPE_mm10_STAR_bw](../B/galaxy-workflows/Galaxy-Workflow-OGlcNac_RNASeqPE_mm10_STAR_bw.ga). The .ga file can be imported in your own galaxy account.

The file to compute the count tables can be downloaded from `wget https://zenodo.org/records/12793186/files/Mus_musculus.GRCm38.102.chr.gtf.tar.gz`

FastQC v0.11.9 was used for quality control: `fastqc --outdir $outfolder --threads $nbcpu --quiet --extract --kmers 7 -f 'fastq' input.fastq.gz`.

Adapters and low quality reads were removed with trim-galore v0.4.3: `trim_galore --phred33 --quality 20  --stringency 1 -e 0.1 --length 20 --output_dir ./ --paired $input_1.fastq.gz $input_2.fastq.gz`.

Alignment was performed with STAR v2.6.0b: `STAR --runThreadN $nbcpu --genomeLoad NoSharedMemory --genomeDir 'mm10/rnastar_index2/mm10/files' --readFilesIn $input_1.fastq.gz $input_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --outSAMstrandField None --outFilterIntronMotifs RemoveNoncanonical --outFilterIntronStrands RemoveInconsistentStrands --outSAMunmapped None --outSAMprimaryFlag OneBestScore --outSAMmapqUnique "255" --outFilterType Normal --outFilterMultimapScoreRange "1" --outFilterMultimapNmax "10" --outFilterMismatchNmax "10" --outFilterMismatchNoverLmax "0.3" --outFilterMismatchNoverReadLmax "1.0" --outFilterScoreMin "0" --outFilterScoreMinOverLread "0.66" --outFilterMatchNmin "0" --outFilterMatchNminOverLread "0.66" --outSAMmultNmax "-1" --outSAMtlen "1" --outBAMsortingBinsN "50"`

Counts were obtained with subread v2.0.1: `featureCounts -a Mus_musculus.GRCm38.102.chr.gtf -F GTF -o $outputcounts.txt -T $nbcpu -s 0 -Q 12 -t 'exon' -g 'gene_id' --minOverlap 1 --fracOverlap 0 --fracOverlapFeature 0 -C input.bam`.

The bigwig files are normalized by TPM and were generated with deeptools v3.0.2:
`bamCoverage --numberOfProcessors $NBCPU  --bam $input.bam --outFileName $ouput.bw --outFileFormat 'bigwig'  --binSize 50  --normalizeUsing BPM`
