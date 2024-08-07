# Differential ATAC-seq in Mouse Embryonic Stem Cells


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


## Installation

Install conda following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Using the recipe [fig1E.yml](fig3D.yml), run:

```
conda env create -n fig3d --file ./fig3D.yml
conda activate fig3d
```


## Figure Generation

The MA plot was obtained with DESeq2 1.22.1. The workflow [OGlcNac_deseq2_atacseqPE_peaks](galaxy-workflow/Galaxy-Workflow-OGlcNac_deseq2_atacseqPE_peaks.ga). Two lists of bam files for each condition (NoDox/Dox) are given as input. The counts are retrieved on the union of peaks (ESC1bNoDox_vs_Dox.gtf) with the featureCounts function of subread v2.0.3: `featureCounts -a ESC1bNoDox_vs_Dox.gtf -F "GTF" -o $output -T $nbcpu -s 0 -Q 0 -t 'peak' -g 'peak_id' --minOverlap 1 --fracOverlap 0 --fracOverlapFeature 0 -p -B -C --countReadPairs $input.bam`. The lists of counts are used as input of DESeq2 1.22.1 that is launched with the script [deseq2.R](../../figure2/B/others/deseq2.R). The deseq2 table was then filtered on columns 3 (fold-change) and 7 (adj p-val) with `abs(c3)>0 and c7<0.05` to keep significant fold-changes that are not equal to zero. The lists of down- and up-regulated genes were separated to their respective files using the third column: `c3<0` for down and `c3>0` for up.


## Pre-processing

### ATAC-seq

The Galaxy workflow [OGlcNac_ATAC_Seq-PEmm10](../../figure1/A/galaxy-workflows/Galaxy-Workflow-OGlcNac_ATAC_Seq-PEmm10.ga) was kindly provided by Charles Girardot from EMBL [GBCS](https://www.embl.org/groups/genome-biology-computational-support/) and was modified for the needs of the study. The .ga file can be imported in your own galaxy account.

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
| ESC1bNoDoxrep1 | YES | 0.001 | 1 | 38 |
| ESC1bNoDoxrep2 | YES | 0.001 | 1 | 38 |
| ESC1bNoDoxrep3 | YES | 0.001 | 1 | 38 |
| ESC1bDoxrep1 | YES | 0.001 | 1 | 38 |
| ESC1bDoxrep2 | YES | 0.001 | 1 | 38 |
| ESC1bDoxrep3 | YES | 0.001 | 1 | 38 |

* Macs2 v2.2.7.1 Broad: `macs2 callpeak -t $input.bam -c NA -n $expname --outdir $outfold -f BAM -g 1.87e9 -s $tagsize --nomodel --extsize 1 --keep-dup $dupthresh --broad --broad-cutoff $qvalue`

### Peak union

Generate the gff file of the union of peaks by running:

```
Rscript union.R
```

The script should give the output:

```
Reading peak files
Reducing intervals
The union returned 106625 peaks
Writing results/ESC1bNoDox_vs_Dox.gtf
```

