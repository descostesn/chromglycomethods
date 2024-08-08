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
repeatsnamevec <- c("LINE", "LTR", "SINE", "Satellite", "DNA", "Simple_repeat", "RNA", "Low_complexity") # nolint
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



.retrieveelementslistofgr <- function(lex, f, f2) {
    return(f(f2(lex)))
}

.retrievefirstexongr <- function(firstexonlist){

    chrvec <- unlist(lapply(firstexonlist, function (x) {
                        .retrieveelementslistofgr(x, as.character, seqnames)}))
    startvec <- unlist(lapply(firstexonlist, function(x) {
                        .retrieveelementslistofgr(x, as.numeric, start)}))
    endvec <- unlist(lapply(firstexonlist, function(x) {
                        .retrieveelementslistofgr(x, as.numeric, end)}))
    strandvec <- unlist(lapply(firstexonlist, function(x) {
                        .retrieveelementslistofgr(x, as.character, strand)}))
    exonidvec <- unlist(lapply(firstexonlist, function(x) mcols(x)$exon_id))
    firstexongr <- GRanges(seqnames = chrvec,
        ranges = IRanges::IRanges(start = startvec, end = endvec),
        strand = strandvec, exonID = exonidvec) # nolint

    return(firstexongr)
}

buildrepeatstarget <- function(txdb, repeatslist, enhancerspath,
    repeatsnamevec) {

    ## Retrieving Promoter, 5' UTR, 3' UTR, Exon, Intron, Downstream
    message("\t Retrieving genomic features: Promoter, 5' UTR, 3' UTR, Exon,",
        " Intron, Downstream")
    promotersgr <- unique(GenomicFeatures::promoters(txdb, upstream = 1000,
                    downstream = 1000))
    intronsgr <- unique(unlist(GenomicFeatures::intronsByTranscript(txdb)))
    exonsgr <- GenomicFeatures::exons(txdb, columns = NULL)
    fiveutrsgr <- unique(unlist(GenomicFeatures::fiveUTRsByTranscript(txdb)))
    threeutrsgr <- unique(unlist(GenomicFeatures::threeUTRsByTranscript(txdb)))

    ## Retrieving the first exons and subtracting it to exons
    message("\t Retrieving first exons")
    allexonsbygenegr <- GenomicFeatures::exonsBy(txdb, by = "gene")
    firstexonlist <- lapply(allexonsbygenegr, function(x) x[1, ])
    firstexongr <- .retrievefirstexongr(firstexonlist)

    ## Removing firstexongr from exonsgr
    exonsgr <- GenomicRanges::setdiff(exonsgr, firstexongr)

    ## Retrieving enhancers
    message("\t Retrieving enhancers")
    enhancersgr <- buildgr(enhancerspath)

    ## Building list of annotations. The order of elements define the priorities
    message("\t Building list of annotations")
    gffnamesvec <- c("promoters", "introns", "firstExons", "otherExons",
            "fiveUTR", "threeUTR", "enhancers")
    annotationslist <- c(repeatslist, promotersgr, intronsgr, firstexongr,
            exonsgr, fiveutrsgr, threeutrsgr, enhancersgr)
    annotationslist <- lapply(annotationslist,
            function(.anno) { mcols(.anno) <- NULL; .anno} )
    names(annotationslist) <- c(repeatsNameVec, gffnamesvec)
    annotationsGRList <- GRangesList(annotationslist)
    
    ## Computing the other locations
    newAnno <- c(unlist(annotationsGRList))
    newAnno.rd <- reduce(trim(newAnno))
    otherLocations <- gaps(newAnno.rd, end=seqlengths(txdb))
    otherLocations <-  otherLocations[strand(otherLocations)!="*"]
    names(otherLocations) <- NULL
    annotationsGRList$otherLocations <- otherLocations
    
    return(annotationsGRList)
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

## Building the GRanges of annotations to which query is compared to
message("Building list of repeats")
repeatslist <- lapply(repeatfilesvec, buildgr, chromvec)
annotationsgrlist <- buildrepeatstarget(txdb, repeatslist, enhancerspath,
    repeatsnamevec)

## Calculate number of annotations
cntRepeats <- lengths(annotationsgrlist)
percentageRepVec <- 100 * cntRepeats / sum(cntRepeats)

## Building colors for piechart
pieColorVec <- c(brewer.pal(n = 12, name = "Paired"), "aliceblue", "azure4", 
        "darkgoldenrod1", "slategray")
names(pieColorVec) <- names(annotationsgrlist)

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
    res <- performOverlap(annotationsgrlist, queryGR)
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
    outputGFFProm(annotationsgrlist, queryGR, peaksIdxByCatPriorList, 
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

