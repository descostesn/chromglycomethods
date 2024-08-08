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
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("GenomicRanges")
library("IRanges")


################
# PARAMETERS
################

queryfilevec <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/makeunion/sept2023/human/glcPolII_samples1-2-3-4/union_glcPolII_sept2023.gff" # nolint
repeatfilesvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/LINE.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/LTR.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/SINE.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/Satellite.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/DNA.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/Simple_repeat.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/RNA.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/repeatmasker/classes/Low_complexity.gff") # nolint
pietitlevec <- "unionPeaksPolIIGlc" # nolint
outputfolder <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/test" # nolint
enhancerspath <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/enhancerAtlas2/DLD1.gff" # nolint
species <- "human"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
chromvec <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chrX", "chrY")
biomartconnection <- "hsapiens_gene_ensembl"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#############
## FUNCTIONS
#############


checkchromosomes <- function(fi, chromvec) {

    idxna <- which(is.na(match(fi$V1, chromvec)))

    if (!isTRUE(all.equal(length(idxna), 0)))
        return(fi[-idxna,])
    else
        return(fi)
}

buildgr <- function(currentpath, chromvec) {

    message("\t Processing ", currentpath)
    fi <- read.table(currentpath, stringsAsFactors = FALSE)
    fi <- checkchromosomes(fi, chromvec) # nolint
    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = fi$V9),
            strand = fi$V7)
    return(gr)
}

.retrieveElementsListOfGR <- function(lEx, f, f2){
    return(f(f2(lEx)))
}


################
# MAIN
################

if (!isTRUE(all.equal(length(queryfilevec), length(pietitlevec))))
    stop("One title per experiment should be given.")

if (!isTRUE(all.equal(length(enhancerspath), 1)))
    stop("Only one list of enhancers is supported.")

if (!file.exists(outputfolder))
        dir.create(outputfolder, recursive = TRUE)

message("Filtering database chromosomes")
seqlevels(txdb) <- chromvec

## Building GR with repeats
message("Building list of repeats")
repeatslist <- lapply(repeatfilesvec, buildgr, chromvec)


## Building the GRanges of annotations to which query is compared to
annotationsGRList <- buildRepeatsTarget(txdb, repeatslist, enhancerspath)
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
ensembl <- tryUseMart(biomart = "ENSEMBL_MART_ENSEMBL",
    biomartconnection, host = "https://nov2020.archive.ensembl.org",
    alternativeMirror = TRUE)

## Determining proportions on each target category for each query file
numbersPieList <- list()
percentagesList <- list()

for(i in seq_len(length(queryfilevec))){
    
    ## Processing query files
    queryFile <- queryfilevec[i]
    nameQuery <- pietitlevec[i]
    outFold <- file.path(outputfolder, nameQuery)
    if(!file.exists(outFold))
        dir.create(outFold, recursive=TRUE)
    
    message("Processing ", nameQuery)
    
    message("\t Building GR with query file")
    queryGR <- unique(buildgr(queryFile))
    
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

names(numbersPieList) <- names(percentagesList) <- pietitlevec

## Perform barplot of all exp
message("Generating barplot for all experiments")
colVec <- c("antiquewhite1", "antiquewhite2", "antiquewhite3", 
        "aquamarine2", "aquamarine3", "aquamarine4", 
        "cadetblue2", "cadetblue3", "cadetblue4", 
        "chocolate1", "chocolate2", "chocolate3", "chocolate4",
        "black")
groupedBarplot(percentagesList, "overlap proportions", outputfolder, percentageRepVec, 
        colVec)
groupedBarplot(numbersPieList, "Nb_of_overlap", outputfolder, cntRepeats, 
        colVec)

