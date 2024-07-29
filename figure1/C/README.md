# Expression of O-GlcNac bound promoters

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  

## Description

Violin plot showing the mouse ESC RNA-seq expression levels of O-GlcNAc enriched gene promoters (236 genes). The median expression of genes having O-GlcNac at promoters is higher when compared to 236 randomly selected promoters without O-GlcNac, 236 randomly selected promoters, and 21,085 promoters. The three groups to which O-GlcNac promoters are compared to have at least 1 RNA-seq read.

The statistical test used is a two-sided Mann-Whitney test (mu = 0, paired = FALSE) giving the p-value:

glcprom vs noglcprom: 0.000692837371518593
glcprom vs randomprom: 0.0356250808138736
glcprom vs allprom: 0.000158394211256996

## Data

Download the following data:

ESCRNAseq_SRR11294181counts.txt
