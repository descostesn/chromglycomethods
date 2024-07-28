#################
## This script aims at reprogramming the venn diagrams in a more sophisticated
## manner integrating the repeats' annotations. Some code was copied/pasted
## from ChIPPeakAnno (3.23.3) function assignChromosomeRegion and from ChIPSeeker
## (1.24.0) function upsetPlot.
## Run on R 4.2.0
# Descostes June 2020
#################


library("GenomicFeatures")
library("ChIPpeakAnno")
library("RColorBrewer")
library("ggplot2")
library("biomaRt")
library("reshape2")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")


################
# PARAMETERS
################


queryfilevec <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff"



#############
## FUNCTIONS
#############

source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/genomicRepartion/grouping_utils_enhancers.R")
#source("grouping_utils.R")


################
# MAIN
################

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## Filtering chromosomes
message("Filtering chromosomes")
seqlevels(txdb) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")

## Building the GRanges of annotations to which query is compared to
annotationsGRList <- buildRepeatsTarget(txdb, repeatsList, enhancerspath)
#save(annotationsGRList, file="/g/boulard/Projects/O-N-acetylglucosamine/analysis/tmp/annotationsGRList.Rdat")
#load("/g/boulard/Projects/O-N-acetylglucosamine/analysis/tmp/annotationsGRList.Rdat")

## Calculate number of annotations
cntRepeats <- lengths(annotationsGRList)
percentageRepVec <- 100 * cntRepeats / sum(cntRepeats)

## Building colors for piechart
pieColorVec <- c(brewer.pal(n = 12, name = "Paired"), "aliceblue", "azure4", 
        "darkgoldenrod1", "slategray")
names(pieColorVec) <- names(annotationsGRList)

message("Connecting to biomart")
if (isTRUE(all.equal(species, "mouse"))){
    ensembl <- tryUseMart(biomart = "ENSEMBL_MART_ENSEMBL",
    "mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org",
    alternativeMirror = TRUE)
} else {
    ensembl <- tryUseMart(biomart="ENSEMBL_MART_ENSEMBL",
    "hsapiens_gene_ensembl", host="https://nov2020.archive.ensembl.org",
    alternativeMirror = TRUE)
}


## Determining proportions on each target category for each query file
numbersPieList <- list()
percentagesList <- list()

for(i in seq_len(length(queryfilevec))){
    
    ## Processing query files
    queryFile <- queryfilevec[i]
    nameQuery <- pieTitleVec[i]
    outFold <- file.path(outputFolder, nameQuery)
    if(!file.exists(outFold))
        dir.create(outFold, recursive=TRUE)
    
    message("Processing ", nameQuery)
    
    message("\t Building GR with query file")
    queryGR <- unique(buildGR(queryFile))
    
    ## Performing overlap on the different categories
    res <- performOverlap(annotationsGRList, queryGR)
    overlapPriority <- res[[1]]
    overlap <- res[[2]]
    annoNamesVec <- res[[3]]
    res <- performPieChart(annoNamesVec, overlapPriority, pieColorVec, outFold, 
            nameQuery,percentageRepVec)
    numbersPieList <- c(numbersPieList, res[1])
    percentagesList <- c(percentagesList, res[2])
    subjectHitsNamesPriority <- res[[3]]
    
    ## Perform an upset diagram
    performUpset(annoNamesVec, overlap, nameQuery, outFold)
    
    ## Output the gff of queryGR per category defined by overlapPriority
    peaksIdxByCatPriorList <- savingPeaksPerCategory(overlapPriority, 
            subjectHitsNamesPriority, queryGR, outFold)
    
    ## Output the gff of the promoters
    outputGFFProm(annotationsGRList, queryGR, peaksIdxByCatPriorList, 
            symbolsTab, ensembl, outFold)
}

names(numbersPieList) <- names(percentagesList) <- pieTitleVec

## Perform barplot of all exp
message("Generating barplot for all experiments")
colVec <- c("antiquewhite1", "antiquewhite2", "antiquewhite3", 
        "aquamarine2", "aquamarine3", "aquamarine4", 
        "cadetblue2", "cadetblue3", "cadetblue4", 
        "chocolate1", "chocolate2", "chocolate3", "chocolate4",
        "black")
groupedBarplot(percentagesList, "overlap proportions", outputFolder, percentageRepVec, 
        colVec)
groupedBarplot(numbersPieList, "Nb_of_overlap", outputFolder, cntRepeats, 
        colVec)

