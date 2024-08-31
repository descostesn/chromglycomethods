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

The Nelf (negative elongation factor) complex interacts with DSIF to repress RNA Pol II from entering the elongation step of the transcriptional process. Nelfcd and Nelfe are subunits of NELF whereas Supt5h is a subunit of DSIF([[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5077205/), [2](https://www.biorxiv.org/content/10.1101/2020.01.23.917237v2), [3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10906531/)]). ints3 is a subunit of the Integrator complex that was shown to regulate transcription elongation and initiation[4](https://pubmed.ncbi.nlm.nih.gov/25201415/). Ncbp1 is also involved in initiation and promotes high-affinity mRNA-cap binding and associates with the CTD of RNA polymerase II[5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10771035/).

