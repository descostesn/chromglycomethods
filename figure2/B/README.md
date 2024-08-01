# Differentially expressed genes upon degradation of OGT by siRNA

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  

## Description

    Volcano plot indicating the differentially expressed genes upon degradation of OGT by siRNA in mouse ES cells.
We found 867 down regulated (dark blue) and 594 upregulated genes (dark red). Among these, 44 down- and 12 up-regulated genes have a fold-change higher than two. 33 down- (light blue) and 2 up-regulated (orange) genes had an O-GlcNac peaks.

!~! check also previous numbers from script

## Data

```
#!/bin/bash

mkdir data

# The tabular file given by DEseq2
wget https://zenodo.org/records/12793186/files/resultDeseq2_siogt.txt  -P data/

## The down- and up-regulated genes having an O-GlcNac peak
wget https://zenodo.org/records/12793186/files/siogtdown_withOGlcNac.gff -P data/
wget https://zenodo.org/records/12793186/files/siogtup_withOGlcNac.gff  -P data/
```

