# Genome-wide O-GlcNac Cut&Run peaks distribution in functional genomic compartments in mouse ES cells

I. [Description](#description)  
II. [Details](#details)  
III. [Data](#data)  
IV. [Installation](#installation)  
V. [Figure Generation](#figure-generation)  
VI. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; VI.I. [Data](#data-1)  
&nbsp;&nbsp; VI.II. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.II.I. [ChIP-seq and CutnRun](#chip-seq-and-cutnrun)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.II.II. [ATAC-seq](#atac-seq)  
&nbsp;&nbsp; VI.III. [Peak Detection](#peak-detection)  



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
wget https://zenodo.org/records/12793186/files/gencode.vM25.annotation.gff -P annotations/
wget https://zenodo.org/records/12793186/files/refGeneUCSC-mm10-March2021.gff -P annotations/
wget https://zenodo.org/records/12793186/files/refseqNCBI-mm10-March2021.gff -P annotations/

# Histone marks
wget https://zenodo.org/records/12793186/files/H3K27ac_SRX19148013_peaks_broadPeak.gff -P data/
wget https://zenodo.org/records/12793186/files/H3K4me1_SRR5466745_peaks_broadPeak.gff -P data/
wget https://zenodo.org/records/12793186/files/H3K27me3_SRR10032683_peaks_hiddendomains.gff -P data/
wget https://zenodo.org/records/12793186/files/H3K4me3_SRX5382140_peaks_broadPeak.gff -P data/
wget https://zenodo.org/records/12793186/files/H3K9me3_SRR925652_peaks_broadPeak.gff -P data/

# Polycomb marks
wget https://zenodo.org/records/12793186/files/Suz12_SRR034190_peaks_broadPeak.gff -P data/
wget https://zenodo.org/records/12793186/files/Ring1B_SRR10095137_peaks_narrowPeak.gff -P data/

# CTD marks
wget https://zenodo.org/records/12793186/files/Ser5P_SRR391050_peaks_broadPeak.gff -P data/
wget https://zenodo.org/records/12793186/files/Ser2P_SRR391039_peaks_broadPeak.gff -P data/

# ATAC-seq
wget https://zenodo.org/records/12793186/files/ATAC_SRR5466767_peaks_narrow.gff -P data/

# O-GlcNac peak replicates
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep1_peaks.gff -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep2_peaks.gff -P data/

# O-GlcNac bigwig replicates
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep1.bw -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14308/ESCHGGlcNAc_rep2.bw -P data/
```


## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1A.yml](fig1A.yml), run:

```
conda env create -n fig1a --file ./fig1A.yml
conda activate fig1a
```

Open an R session and install our in-house package [genomecompR](https://github.com/descostesn/genomecompR/tree/main):

```
if (!require("devtools")) install.packages("devtools")
devtools::install_github("descostesn/genomecompR")
```


## Figure Generation

Run the script [upsetgenomecomp.R](upsetgenomecomp.R) from the current folder, it uses the files downloaded in the subfolders `data/` and `annotations/`. From the terminal:

```
Rscript upsetgenomecomp.R
```

It should output the following numbers:

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

<img src="complexUpset.png" alt="Upset plot" width="400"/>

## Pre-processing

### Data

| Target | ID | library layout | link |
|--------|----|----------------|------|
| O-GlcNAc rep1 | E-MTAB-14308 | single | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/033/ERR13430733/ERR13430733.fastq.gz |
| O-GlcNAc rep2 | E-MTAB-14308 | single | ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/034/ERR13430734/ERR13430734.fastq.gz |
| H3K27me3 | SRR10032683 | paired | see the provided [snakemake](snakemake/Snakefile) |
| H3K27ac | SRR23199099 | single | see the provided [snakemake](snakemake/Snakefile) |
| H3K4me1 | SRR5466745 | single | see the provided [snakemake](snakemake/Snakefile) |
| H3K4me3 | SRR8581420 | paired | see the provided [snakemake](snakemake/Snakefile) |
| Suz12   | SRR034190 | single | see the provided [snakemake](snakemake/Snakefile) |
| RING1B | SRR10095137 | single | see the provided [snakemake](snakemake/Snakefile) |
| H3K9me3 | SRR925652 | single | see the provided [snakemake](snakemake/Snakefile) |
| Ser5P | SRR391050 | single | see the provided [snakemake](snakemake/Snakefile) |
| Ser2P | SRR391039 | single | see the provided [snakemake](snakemake/Snakefile) |
| ATAC-seq | SRR5466767 | paired | see the provided [snakemake](snakemake/Snakefile) |


### Workflows

### ChIP-seq and CutnRun

The pre-processing was performed with the Galaxy workflows [OGlcNac_ChIP-SeqPEmm10](galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqPEmm10.ga) and [OGlcNac_ChIP-SeqSEmm10](galaxy-workflows/Galaxy-Workflow-OGlcNac_ChIP-SeqSEmm10.ga) for paired-end and single-end data respectively. The .ga files can be imported in one own galaxy account.

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


### ATAC-seq

The Galaxy workflow [OGlcNac_ATAC_Seq-PEmm10](galaxy-workflows/Galaxy-Workflow-OGlcNac_ATAC_Seq-PEmm10.ga) was kindly provided by Charles Girardot from EMBL [GBCS](https://www.embl.org/groups/genome-biology-computational-support/) and was modified for the needs of the study. The .ga file can be imported in your own galaxy account.

Adapters and low quality reads were removed with trim-galore v0.4.3: `trim_galore --phred33 --output_dir ./ --paired input_1.fastq.gz input_2.fastq.gz`

Reads were aligned to the mm10 genome with bowtie2 v2.3.4 and sorted with samtools v1.8: `bowtie2  -p $nbcpu -x m.musculus/mm10/mm10 -1 input_f.fastq.gz -2 input_r.fastq.gz --un-conc-gz $unaligned.fastq.gz -I 0 -X 2000 --fr --dovetail --sensitive 2> $ouput.log | samtools sort -@$nbcpu -O bam -T $TMPDIR -o $output.bam`

Only reads with valid primary alignments were kept with samtools v1.2: `samtools view -o $output.bam -h -b -q 20 -f 0x3 -F 0x800 input.bam`. Reads aligned to non-canonical chromosomes were filtered out with samtools v1.2: `samtools view -o $ouput.bam -h  -b $input.bam 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chrX' 'chrY'`. Duplicates were removed with picard v2.7.1: `picard MarkDuplicates  INPUT=$input.bam OUTPUT=$ouput.bam METRICS_FILE=$metrics.txt REMOVE_DUPLICATES='true' ASSUME_SORTED='true' DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'`.

To perform quality controls, from the output of Bowtie2, duplicated alignments were marked with picard v2.7.1: `picard MarkDuplicates INPUT=$input.bam OUTPUT=$output.bam METRICS_FILE=$metrics.txt REMOVE_DUPLICATES='false' ASSUME_SORTED='true' DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR`. Bam files were then sorted with picard v2.7.1: `picard ReorderSam INPUT=$input.bam OUTPUT=$ouput.bam REFERENCE="M.musculus/mm10/fasta/mm10.fa" ALLOW_INCOMPLETE_DICT_CONCORDANCE="true" ALLOW_CONTIG_LENGTH_DISCORDANCE="false" VALIDATION_STRINGENCY="LENIENT" QUIET=true VERBOSITY=ERROR`

FastQC v0.11.9 was run after the Bowtie2 alignment, marking, and removing duplicates ([rgFastQC.py](others/rgFastQC.py)): `rgFastQC.py -i $input.bam -d $outputdir -o $htmloutput -t $output.txt -f "bam" -j "Valid Uniquely Mapped Reads"`
Statistics on the reordered bam were obtained with samtools v1.9: `samtools flagstat $input.bam > $output.txt`
Alignment summary statistics were retrieved with picard v2.7.1: `picard CollectAlignmentSummaryMetrics INPUT=$input.bam OUTPUT=$outputstats.txt MAX_INSERT_SIZE=5000 METRIC_ACCUMULATION_LEVEL="ALL_READS" IS_BISULFITE_SEQUENCED="false" REFERENCE_SEQUENCE="M.musculus/mm10/fasta/mm10.fa" ASSUME_SORTED="true"  VALIDATION_STRINGENCY="LENIENT" QUIET=true VERBOSITY=ERROR`
Insert size metrics were obtained with picard v2.7.1: `picard CollectInsertSizeMetrics INPUT=$input.bam OUTPUT=$output.pdf HISTOGRAM_FILE=$hist.txt DEVIATIONS="10.0"   MINIMUM_PCT="0.05" REFERENCE_SEQUENCE="M.musculus/mm10/fasta/mm10.fa" ASSUME_SORTED="true" METRIC_ACCUMULATION_LEVEL="ALL_READS" VALIDATION_STRINGENCY="LENIENT" QUIET=true VERBOSITY=ERROR`

The log files of bowtie2, remove duplicates, FastQC, alignment summary, and insert size metrics were sent to MultiQC v1.7 to perform an overall quality assessment: `multiqc multiqc_WDir --filename "report"`

Alignments having an insert size lower or equal to 100 bp were removed with bamtools v2.4.0: `bamtools filter -in $input.bam -out $output.bam -length 100`



### Peak Detection

**MACS2**

| Target | Broad | q-value | Duplicates Thres. | Tag size |
|--------|-------|---------|-------------------|----------|
| H3K27ac | YES | 0.04 | 4 | 51 |
| H3K4me1 | YES | 0.04 | 6 | 76 |
| H3K4me3 | YES | 0.04 | 4 | 69 |
| Suz12   | YES | 0.04 | 4 | 36 |
| RING1B | NO | 0.001 | 1 | 75 |
| H3K9me3 | YES | 0.04 | 12 | 119 |
| Ser5P | YES | 1e-04 | 1 | 36 |
| Ser2P | YES | 0.03 | 1 | 51 |
| ESCHGGlcNAc_rep1 | NO | 0.04 | 7 | 82 |
| ESCHGGlcNAc_rep2 | NO | 0.04 | 7 | 82 |
| ATAC-seq (Galaxy) | NO | 0.1 | NR | NR |


* Macs2 v2.2.7.1 Broad: `macs2 callpeak -t $input.bam -c $control.bam -n $expname --outdir $outfold -f BAM -g 1.87e9 -s $tagsize --nomodel --extsize 150 --keep-dup $dupthresh --broad --broad-cutoff $qvalue`
* Macs2 v2.2.7.1 Narrow: `macs2 callpeak -t $input.bam -c $control.bam -n $expname --outdir $outfold -f BAM -g 1.87e9 -s $tagsize -q $qvalue --nomodel --extsize 150 --keep-dup $dupthres`
* Galaxy Macs2 v2.1.1.20160309: `macs2 callpeak -t $input.bam --name $expname --format BAMPE --gsize 1.87e9 --keep-dup '1' --qvalue '0.1' --nomodel --extsize '75' --shift '0'`


**hiddenDomains v3.1**

Used to detect H3K27me3 peaks: `hiddenDomains -g $chromInfoFile -b 300 -t $input.bam -c $control.bam -o $outfold`
