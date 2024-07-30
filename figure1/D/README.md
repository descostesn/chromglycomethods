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

The 527 O-GlcNac peaks were submitted to the [ChIP-Atlas](https://chip-atlas.org/) database to find overlap with ChIP-seq experiments using the [enrichment analysis](https://chip-atlas.org/enrichment_analysis) tool. The online tool is rapidly evolving and new experiments submitted to [GEO](https://www.ncbi.nlm.nih.gov/geo/) are regularly added. We provide the results of the enrichment analysis in the folder [chipatlas_results](chipatlas_results/).

O-GlcNac peaks mainly overlap with RNA Polymerase II (RNAPol II) factors: [TBP](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TBP), [Taf12](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TAF12&keywords=Taf12), [NelfA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NELFA&keywords=nelfa), [Med1](https://genecards.org/cgi-bin/carddisp.pl?gene=MED1&keywords=Med1), [Med12](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MED12&keywords=med12), [Med24](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MED24&keywords=med24), [Med26](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MED26&keywords=med26), and [Dr1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DR1&keywords=dr1). TBP (TATA-binding proteins) is part of the General Transcription Factor IID that coordinates the binding of RNAPol II to the core promoter. Taf12 directly interacts with TBP. NelfA (Negative Elongation Factor Complex Member A) is part of the NELF complex that regulates transcription elongation by RNAPol II. Med1, Med12, Med24, and Med26 are subunits of the Mediator Complex which, after its recruitment to the promoters, serves as a scaffold for the assembly of the RNA PolII pre-initiation complex. Dr1 (Down-Regulator Of Transcription 1) represses both basal and activated levels of transcription by RNAPol II via its interaction with TBP. Other proteins that are known to carry O-GlcNac are highlighted in green.

## Details

The results given by ChIP-Atlas are not manually curated. We manually identified the following problematic datasets:

See [problematic-rep1.txt](replacement-files/problematic-rep1.txt) and [problematic-rep2.txt](replacement-files/problematic-rep2.txt)

These datasets were replaced by the following valid one which were found in the sorted lists of the scripts. In other words, we looked manually for the next candidates that had a lower overlap (but still > 20%). NA means that no suitable candidate was found:

See [replacement-rep1.txt](replacement-files/replacement-rep1.txt) and [replacement-rep2.txt](replacement-files/replacement-rep2.txt)

## Data
## Installation
## Figure Generation

The script should output:

```
Reading ChIP-Atlas results
Reading mHG1hist
  Filtering lines having more overlap on random control
  Filtering on the Qvalue
  Filtering out the mock IP expriments
  Keeping 7300/7330(100%)
Reading mHG1pol
  Filtering lines having more overlap on random control
  Filtering on the Qvalue
  Filtering out the mock IP expriments
  Keeping 774/776(100%)
Reading mHG1TF
  Filtering lines having more overlap on random control
  Filtering on the Qvalue
  Filtering out the mock IP expriments
  Keeping 6772/6784(100%)
Reading mHG1hist
  Filtering lines having more overlap on random control
  Filtering on the Qvalue
  Filtering out the mock IP expriments
  Keeping 7301/7330(100%)
Reading mHG2pol
  Filtering lines having more overlap on random control
  Filtering on the Qvalue
  Filtering out the mock IP expriments
  Keeping 774/776(100%)
Reading mHG2TF
  Filtering lines having more overlap on random control
  Filtering on the Qvalue
  Filtering out the mock IP expriments
  Keeping 6779/6785(100%)
```
