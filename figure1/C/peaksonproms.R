#################
## This script aims at reprogramming the venn diagrams in a more sophisticated
## manner integrating the repeats' annotations. Some code was copied/pasted
## from ChIPPeakAnno (3.23.3) function assignChromosomeRegion and from ChIPSeeker
## (1.24.0) function upsetPlot.
## Run on R 4.2.0
# Descostes June 2020
#################


# library("ChIPpeakAnno")
# library("RColorBrewer")
# library("ggplot2")
# library("reshape2")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("GenomicFeatures")
library("GenomicRanges")
library("S4Vectors")
library("biomaRt")

################
# PARAMETERS
################


queryfile <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff" # nolint
usealtmirror <- TRUE

#############
## FUNCTIONS
#############

buildgr <- function(currentpath) {
    message("\t Processing ", currentpath)
    fi <- read.table(currentpath, stringsAsFactors = FALSE)
    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = fi$V9), strand = fi$V7)
    return(gr)
}

## Function by Ilyess Rachedi
tryusemart <- function(biomart = "ensembl", dataset, host, alternativemirror) {
    c <- 1
    repeat {
        message("# Attempt ", c, "/5 # Connection to Ensembl ... ")
        if (!alternativemirror)
            ensembl <- try(biomaRt::useMart(biomart, dataset = dataset,
                host = host), silent = TRUE)
        else
            ensembl <- try(biomaRt::useEnsembl(biomart, dataset = dataset,
                host = host, mirror = "useast"), silent = TRUE)

        if (isTRUE(is(ensembl, "try-error"))) {
            c <- c + 1
            error_type <- attr(ensembl, "condition")
            message(error_type$message)

            if (c > 5)
                stop("There is a problem of connexion to Ensembl for now. ",
                "Please retry later or set alternativemirror=TRUE.")
        }else {
            message("Connected with success.")
            return(ensembl)
        }
    }
}

################
# MAIN
################

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## Filtering chromosomes
message("Filtering chromosomes")
seqlevels(txdb) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")

## Retrieve promoters coordinates
promotersgr <- unique(GenomicFeatures::promoters(txdb, upstream = 1000,
                    downstream = 1000))
querygr <- unique(buildgr(queryfile))

## Perform overlap between O-GlcNac peaks and promoters
resultoverlap <- GenomicRanges::findOverlaps(querygr, promotersgr,
    ignore.strand = FALSE)
idxkeep <- which(!duplicated(S4Vectors::queryHits(resultoverlap)))
resultoverlap <- resultOverlap[idxkeep, ]
promotersgr <- promotersgr[S4Vectors::subjectHits(result), ]

## Retrieving information on genes
ensembl <- tryusemart(biomart = "ENSEMBL_MART_ENSEMBL",
    "mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org",
    alternativemirror = usealtmirror)





!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    alternativemirror = TRUE)
} else {
    ensembl <- tryUseMart(biomart="ENSEMBL_MART_ENSEMBL",
    "hsapiens_gene_ensembl", host="https://nov2020.archive.ensembl.org",
    alternativemirror = TRUE)
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

