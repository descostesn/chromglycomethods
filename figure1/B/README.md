# O-GlcNac occupies promoters and transposable elements

I. [Description](#description)  
II. [Data](#data)  
III. [Figure Generation](#figure-generation)  
IV. [Pre-processing](#pre-processing)  


## Description

Genome browser view of O-GlcNac, H3K4me3, and H3K9me3 at the Setmar and Nphp3 genes. LINE1 retrotransposons are also indicated. O-GlcNac co-localizes at the gene promoters alongside H3K4me3 and also colocalizes at transposons with H3K9me3.

## Data

The bigwig files used in this figure can be obtained with:

```
#!/bin/bash

wget https://zenodo.org/records/13444099/files/ESCHGGlcNAc1_lane1sample12_multik10.bw
wget https://zenodo.org/records/13444099/files/H3K4me3_SRX5382140_multik10.bw
wget https://zenodo.org/records/13444099/files/H3K9me3_SRR925652_multik10.bw
```
 
## Figure Generation

The bw files were uploaded to [IGV](https://igv.org/) v2.13.0 selecting the mouse mm10 genome.


## Pre-processing

The ESCHGGlcNAc and H3K9me3 data was processed with [Galaxy-Workflow-ChIP-SeqSEmm10_multi10](galaxy-workflows/Galaxy-Workflow-ChIP-SeqSEmm10_multi10.ga) and H3K4me3 with [Galaxy-Workflow-ChIP-SeqPEmm10_multi10](galaxy-workflows/Galaxy-Workflow-ChIP-SeqPEmm10_multi10.ga) allowing a maximum of 10 matches for multireads. The multiread mode allowed retrieving signal on LINE1 elements.

Quality control was done with FastQC v0.11.9: `fastqc --outdir $outputfolder --threads $nbcpu --quiet --extract --kmers 7 -f 'fastq' $input.fastq.gz`.

Adapters and low quality reads were removed with trim-galore v0.4.3: `trim_galore --phred33 --quality 20 --stringency 1 -e 0.1 --length 20 --output_dir ./ $input.fastq.gz`.

Reads were aligned to mm10 with Bowtie 2.3.4.1 using the `-k 10` option and the bam were sorted using samtools v1.9:
single: `bowtie2  -p $nbcpu -x m.musculus/mm10/mm10' -U $input.fastq.gz --skip 0 --qupto 100000000 --trim5 0 --trim3 0 --phred33 -N 0 -L 22 -i 'S,1,1.15' --n-ceil 'L,0,0.15' --dpad 15 --gbar 4 --end-to-end --score-min 'L,-0.6,-0.6'  --mp '6,2' --np 1 --rdg 5,3 --rfg 5,3  -k 10 -D 15 -R 2 --seed 0 --no-unal 2> $log | samtools sort -@$nbcpu -O bam -o $output.bam`
paired: `bowtie2  -p $nbcpu -x m.musculus/mm10/mm10 -1 $input_f.fastq.gz -2 $input_r.fastq.gz -I 0 -X 500 --fr   --dovetail -N 0 -L 22 -i 'S,1,1.15' --n-ceil 'L,0,0.15' --dpad 15 --gbar 4 --end-to-end --score-min 'L,-0.6,-0.6' -k 10 --no-unal 2> $log | samtools sort -@$nbcpu -O bam -o $output.bam`

Only primary alignments were kept using samtools v1.9: `samtools view -o $output.bam -h -b -q 20 -F 0x800 $input.bam`.

Reads not aligned to consensus chromosomes were excluded with samtools v1.9: `samtools view -o $output.bam -h -b $input.bam 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chrX' 'chrY'`.

The bigwig files were generated normalizing by the effective genome size using Deeptools 3.0.2: `bamCoverage --numberOfProcessors $nbcpu  --bam $input.bam --outFileName $output.bw --outFileFormat 'bigwig' --binSize 50 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 --scaleFactor 1.0 --extendReads 150 --minMappingQuality 1`
