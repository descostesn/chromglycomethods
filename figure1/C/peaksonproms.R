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
library("Rargs")
library("biomaRt")
library("reshape2")



################
# PARAMETERS
################

paramsDefinition <- list()

paramsDefinition[["--queryFileVec"]] <- list(variableName="queryFileVec", numeric=F, mandatory=T, description="Path to the GFF file containing the peaks.")
paramsDefinition[["--repeatFilesVec"]] <- list(variableName="repeatFilesVec", numeric=F, mandatory=T, description="Vector containg the GFF paths of the repeats.")
paramsDefinition[["--repeatsNameVec"]] <- list(variableName="repeatsNameVec", numeric=F, mandatory=T, description="Vector containing the corresponding names of the repeats.")
paramsDefinition[["--pieTitleVec"]] <- list(variableName="pieTitleVec", numeric=F, mandatory=T, description="Title to display for the peaks on the piechart.")
paramsDefinition[["--outputFolder"]] <- list(variableName="outputFolder", numeric=F, mandatory=T, description="Path to the output folder.")
paramsDefinition[["--enhancerspath"]] <- list(variableName="enhancerspath", numeric=F, mandatory=T, description="Path to the file containing the enhancer coordinates. GFF format.")
## Optional
paramsDefinition[["--species"]] <- list(variableName="species", numeric=F, mandatory=F, description="Should be mouse or human.", default = "mouse")

queryFileVec <- 
c(paste("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG6hGlcNAc1_lane1sample15_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG6hGlcNAc2_lane1sample16_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG6hGlcNAc3_lane1sample17_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc2_lane1sample13_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc3_lane1sample14_peaks_narrowPeak.gff", sep = " "),

paste("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG48hGlcNAc1_lane1sample18_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG48hGlcNAc2_lane1sample19_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG48hGlcNAc3_lane1sample20_peaks_narrowPeak.gff", sep = " "),

paste("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG1wGlcNAc1_lane1sample21_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG1wGlcNAc2_lane1sample22_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCLG1wGlcNAc3_lane1sample23_peaks_narrowPeak.gff", sep = " "))

repeatFilesVec <- paste("/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/LINE.gff", 
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/LTR.gff",
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/SINE.gff", 
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/Satellite.gff",
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/DNA.gff",
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/Simple_repeat.gff",
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/RNA.gff",
		"/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/Low_complexity.gff", sep=" ")

repeatsNameVec <- paste("LINE", "LTR", "SINE", "Satellite", "DNA", "Simple_repeat", "RNA", "Low_complexity", sep=" ")

pieTitleVec <- c(paste("LG6hGlcNAc1s15",
"LG6hGlcNAc2s16",
"LG6hGlcNAc3s17",
"HGGlcNAc1s12",
"HGGlcNAc2s13",
"HGGlcNAc3s14", sep = " "),

paste("LG48hGlcNAc1s18",
"LG48hGlcNAc2s19",
"LG48hGlcNAc3s20", sep = " "),

paste("LG1wGlcNAc1s21",
"LG1wGlcNAc2s22",
"LG1wGlcNAc3s23", sep = " "))


outputFolder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_mouseGlucose_enhancers/",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_mouseGlucose_enhancers/",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_mouseGlucose_enhancers/")

enhancerspath <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/enhancerAtlas2/E14.gff"




# Retreives the parameters
getParams(paramsDefinition)



#############
## FUNCTIONS
#############

source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/genomicRepartion/grouping_utils_enhancers.R")
#source("grouping_utils.R")


################
# MAIN
################

if (!isTRUE(all.equal(length(queryFileVec), length(pieTitleVec))))
    stop("One title per experiment should be given.")

if (!isTRUE(all.equal(length(enhancerspath), 1)))
    stop("Only one list of enhancers is supported.")

if (!file.exists(outputFolder))
        dir.create(outputFolder, recursive=TRUE)

if (!isTRUE(all.equal(species, "mouse")) &&
    !isTRUE(all.equal(species, "human")))
    stop("species should be mouse or human")

if (isTRUE(all.equal(species, "mouse"))) {
    library("TxDb.Mmusculus.UCSC.mm10.knownGene")
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
} else {
    library("TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
}

## Building GR with repeats
message("Building list of repeats")
repeatsList <- lapply(repeatFilesVec, buildGR)
#save(repeatsList, file="/g/boulard/Projects/O-N-acetylglucosamine/analysis/tmp/repeatsList.Rdat")
#load("/g/boulard/Projects/O-N-acetylglucosamine/analysis/tmp/repeatsList.Rdat")

## Filtering chromosomes
message("Filtering chromosomes")
if (isTRUE(all.equal(species, "mouse"))) {
    seqlevels(txdb) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
		"chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
		"chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
} else {
    seqlevels(txdb) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
		"chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
		"chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chrX", "chrY")
}

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

for(i in seq_len(length(queryFileVec))){
    
    ## Processing query files
    queryFile <- queryFileVec[i]
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

