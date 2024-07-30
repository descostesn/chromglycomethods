# Heatmap of O-GlcNac peaks along with ChIP-Atlas candidates

I. [Description](#description)  
II. [Data](#data)  
III. [Installation](#installation)  
IV. [Figure Generation](#figure-generation)  
V. [Pre-processing](#pre-processing)  
&nbsp;&nbsp; V.I. [Workflows](#workflows)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; VI.I.I. [ChIP-seq and CutnRun](#cutnrun)  


## Description

Heatmap of CuntRun O-GlcNac signal and ChIP-Seq ChIP-Atlas selected candidates signal 1 kb around the start of 1,200 O-GlcNac peaks. The list of peaks was obtained computing the union of the replicate peaks. Three groups were defined and highlight the RNA Polymerase II (RNAPol II) before, at, and after the O-GlcNac peaks. For the groups I and II, RNAPol II describes a double peak signal overlapping O-GlcNac independently of its position. Other selected candidates follow a similar pattern. It indicates that a sub-category of these proteins might carry O-GlcNac in a transient state. The description of the function of each protein can be found in [fig1D](../D/README.md).

## Data

	
mouseESC_fig1E_peakorder.bed
union_sept2023mouse_HG.gff