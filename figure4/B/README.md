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
