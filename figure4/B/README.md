# Heatmap of the percentage of overlap of the ChIP-Atlas candidates with O-GlcNac peaks in mouse ESC

I. [Description](#description)  
II. [Details](#details)  
III. [Data](#data)  
IV. [Installation](#installation)  
V. [Figure Generation](#figure-generation)  
VI. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; VI.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [CutnRun](#cutnrun)  

## Description

The 3450 (nodox rep1), 4201 (nodox rep2), 2630 (dox rep1), and 1848 (dox rep2) O-GlcNac peaks were submitted to the [ChIP-Atlas](https://chip-atlas.org/) database to find overlap with ChIP-seq experiments using the [enrichment analysis](https://chip-atlas.org/enrichment_analysis) tool. The online tool is rapidly evolving and new experiments regularly added to [GEO](https://www.ncbi.nlm.nih.gov/geo/). We provide the results of the enrichment analysis in the folder [chipatlas_results](chipatlas_results/).

O-GlcNac peaks mainly overlap with RNA Polymerase II (RNAPol II) and associated factors: [Nelfcd](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NELFCD), [Supt5h](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SUPT5H), [Ncbp1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NCBP1), [Nelfe](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NELFE), and [Ints3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=INTS3).

The Nelf (negative elongation factor) complex interacts with DSIF to repress RNA Pol II from entering the elongation step of the transcriptional process. Nelfcd and Nelfe are subunits of NELF whereas Supt5h is a subunit of DSIF([[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5077205/), [2](https://www.biorxiv.org/content/10.1101/2020.01.23.917237v2), [3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10906531/)]). Ints3 is a subunit of the Integrator complex that was shown to regulate transcription elongation and initiation[[4](https://pubmed.ncbi.nlm.nih.gov/25201415/)]. Ncbp1 is also involved in initiation and promotes high-affinity mRNA-cap binding and associates with the CTD of RNA polymerase II[[5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10771035/)].

## Details

The results given by ChIP-Atlas are not manually curated. We manually identified the following problematic datasets:

See [problematic-nodoxrep1.txt](replacement-files/problematic-nodoxrep1.txt), [problematic-nodoxrep2.txt](replacement-files/problematic-nodoxrep2.txt), [problematic-doxrep1.txt](replacement-files/problematic-doxrep1.txt) and [problematic-doxrep2.txt](replacement-files/problematic-doxrep2.txt).

These datasets were replaced by the following valid one which were found in the sorted lists of the scripts. In other words, we looked manually for the next candidates that had a lower overlap (but still > 20%). NA means that no suitable candidate was found:

See [replacement-nodoxrep1.txt](replacement-files/replacement-nodoxrep1.txt), [replacement-nodoxrep2.txt](replacement-files/replacement-nodoxrep2.txt), [replacement-doxrep1.txt](replacement-files/replacement-doxrep1.txt) and [replacement-doxrep2.txt](replacement-files/replacement-doxrep2.txt).

## Data

As explained in the description section, The online tool is rapidly evolving and new experiments are regularly added [GEO](https://www.ncbi.nlm.nih.gov/geo/). We provide the results of the enrichment analysis in the folder [chipatlas_results](chipatlas_results/).

If one whishes to obtain results with the updated database, the O-GlcNac peaks can be obtained at:

```
#!/bin/bash

mkdir data

# O-GlcNac peaks
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcDoxAux_rep1_peaks.gff -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcDoxAux_rep2_peaks.gff -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcNoDoxAux_rep1_peaks.gff -P data/
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-14307/DLD1GlcNAcNoDoxAux_rep2_peaks.gff -P data/
```
