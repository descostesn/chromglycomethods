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
library("GenomeInfoDb")
library("S4Vectors")
library("tibble")
library("ggupset")

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
    names(annotationslist) <- c(repeatsnamevec, gffnamesvec)
    annotationsgrlist <- GenomicRanges::GRangesList(annotationslist)

    ## Computing the other locations
    newanno <- c(unlist(annotationsgrlist))
    newannord <- reduce(trim(newanno))
    otherlocations <- gaps(newannord, end = GenomeInfoDb::seqlengths(txdb))
    otherlocations <-  otherlocations[strand(otherlocations) != "*"]
    names(otherlocations) <- NULL
    annotationsgrlist$otherlocations <- otherlocations

    return(annotationsgrlist)
}


performoverlap <- function(annotationsgrlist, querygr) {

    #############
    ## Performing the overlap and determining preferencce order
    ## The order of repeats matters since a peak could overlap two different
    ## repeats if they are close enough. Note the repeats in mouse do not
    ## overlap between them.
    #############

    annonamesvec <- names(annotationsgrlist)
    mcols(querygr) <- NULL
    resultoverlap <- GenomicRanges::findOverlaps(querygr, annotationsgrlist,
        ignore.strand = FALSE)
    idxkeep <- which(!duplicated(S4Vectors::queryHits(resultoverlap)))
    resultoverlappriority <- resultoverlap[idxkeep, ]

    if (isTRUE(all.equal(length(resultoverlappriority), 0)) ||
            length(resultoverlappriority) > length(resultoverlap))
        stop("Pb in the script.")

    if (!isTRUE(all.equal(length(resultoverlappriority), length(querygr))))
        stop("All peaks of query should have a unique mapping location, pb in ",
                "the script")

    return(list(resultoverlappriority, resultoverlap, annonamesvec))
}


performpiechart <- function(annonamesvec, overlappriority, piecolorvec, outfold,
    namequery, percentagerepvec, querygr) {

    ## Calculate number of overlaps
    subjecthitsnamespriority <- annonamesvec[subjectHits(overlappriority)]
    counts <- table(subjecthitsnamespriority)
    percentagevec <- 100 * counts / length(querygr)

    if (!isTRUE(all.equal(sum(percentagevec), 100)))
        stop("For the priority piechart, the percentages sum is not equal to ",
                "100: ", sum(percentagevec))

    piecolorvechits <-  piecolorvec[names(percentagevec)]

    ## Plotting the priority piechart
    png(file = file.path(outfold, paste0(namequery, "priorityPiechart.png")),
            width = 10, height = 10)
    par(mfrow = c(1, 2))
    pie(percentagevec, labels = names(percentagevec), col = piecolorvechits,
            main = namequery)
    pie(percentagerepvec, labels = names(percentagerepvec), col = piecolorvec,
            main = "Proportions")
    dev.off()

    return(list(counts, percentagevec, subjecthitsnamespriority))
}


performupset <- function(annonamesvec, overlap, namequery, outfold){

    message("Plotting upset")
    subjecthitsnames <- annonamesvec[subjectHits(overlap)]
    regionsperpeaklist <- split(subjecthitsnames, queryHits(overlap))
    regionslist <- lapply(regionsperpeaklist, function(currentpeakregions,
        namescategories) {
            isvec <- rep(FALSE, length(namescategories))
            names(isvec) <- namescategories
            currentpeakregions <- unique(currentpeakregions)
            isvec[currentpeakregions] <- TRUE
            return(isvec)
        }, annonamesvec)
    regionsmatrix <- do.call(rbind, regionslist)

    if (!isTRUE(all.equal(nrow(regionsmatrix), length(regionsperpeaklist))) ||
            !isTRUE(all.equal(ncol(regionsmatrix), length(annonamesvec))))
        stop("regionsmatrix does not have the correction dimensions.")

    namescolregions <- colnames(regionsmatrix)
    res <- tibble::tibble(anno = lapply(seq_len(nrow(regionsmatrix)),
                    function(i) namescolregions[regionsmatrix[i, ]]))

    g <- ggplot2::ggplot(res, aes_(x = ~anno)) + ggplot2::geom_bar() +
            ggplot2::xlab(NULL) + ggplot2::ylab(NULL) +
            ggplot2::theme_minimal() +
            ggupset::scale_x_upset(n_intersections = 20, order_by = "freq")

    ggsave(filename = paste0(namequery, "-priorityUpSet.png"), plot = g,
            device = png(), path = outfold)
}



.restrictannooverlap <- function(querygr, peaksidxbycatpriorlist,
    currentannogr) {

    result <- GenomicRanges::findOverlaps(
        querygr[peaksidxbycatpriorlist$promoters, ],
            currentannogr, ignore.strand = FALSE)
    result <- result[which(!duplicated(queryHits(result))), ]
    currentannogr <- currentannogr[subjectHits(result), ]
    return(currentannogr)
}


## Function by Ilyess Rachedi
.trygetbm <- function(attributes, ensembl, values = NULL, filters = NULL) {
    c <- 1
    repeat {
        message("# Attempt ", c, "/5 # ",
                "Retrieving information about genes from biomaRt ...")
        if (is.null(values) && is.null(filters))
            res <- try(biomaRt::getBM(attributes = attributes, mart = ensembl),
                silent = TRUE)
        else
            res <- try(biomaRt::getBM(attributes = attributes, mart = ensembl,
                values = values, filters = filters), silent = TRUE)

        if( isTRUE(is(res, "try-error"))) {
            c <- c + 1
            error_type <- attr(res, "condition")
            message(error_type$message)
            if (c > 5)
                stop("There is a problem of connexion to Ensembl for ",
                        "now. Please retry later.")
        } else {
            message("Information retrieved with success.")
            return(res)
        }
    }
}

.retrievegeneinfo  <- function(ensembl, currentannogr) {

    attributes <- c("chromosome_name", "ensembl_gene_id",
            "ensembl_transcript_id_version", "external_gene_name",
            "start_position", "end_position", "strand", "transcript_start",
            "transcript_end", "transcription_start_site")
    symbolstab <- .trygetbm(attributes, ensembl, values = names(currentannogr),
            filters = "ensembl_transcript_id_version")
    symbolstab$strand[which(symbolstab$strand == 1)] <- "+"
    symbolstab$chromosome_name <- paste0("chr",symbolstab$chromosome_name)
    if (!isTRUE(all.equal(length(which(symbolstab$strand == -1)), 0)))
        symbolstab$strand[which(symbolstab$strand == -1)] <- "-"
    return(symbolstab)
}

outputgffprom <- function(annotationsgrlist, querygr, peaksidxbycatpriorlist,
        symbolstab, ensembl, outfold) {

    currentannogr <- annotationsgrlist$promoters
    currentannogr <- .restrictannooverlap(querygr, peaksidxbycatpriorlist,
            currentannogr)
    symbolstab <- .retrievegeneinfo(ensembl, currentannogr)
    idxTable <- match(names(currentannogr), 
            symbolstab$ensembl_transcript_id_version)
    idxNA <- which(is.na(idxTable))
    currentannogr <- currentannogr[-idxNA,]
    idxTable <- idxTable[-idxNA]
    ## Transcripts
    writePromWithUnique(IDVec = symbolstab$ensembl_transcript_id_version[idxTable],
            startvec = symbolstab$transcript_start[idxTable],
            endvec = symbolstab$transcript_end[idxTable], 
            symbolstab = symbolstab, idxTable = idxTable, outfold = outfold, 
            filename = "transcripts_fromProm")
    ## Genes
    writePromWithUnique(IDVec = symbolstab$ensembl_gene_id[idxTable],
            startvec = symbolstab$start_position[idxTable],
            endvec = symbolstab$end_position[idxTable], symbolstab = symbolstab,
            idxTable =  idxTable, outfold = outfold, filename = "genes_fromProm")
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
cntrepeats <- lengths(annotationsgrlist)
percentagerepvec <- 100 * cntrepeats / sum(cntrepeats)

## Building colors for piechart
piecolorvec <- c(brewer.pal(n = 12, name = "Paired"), "aliceblue", "azure4",
        "darkgoldenrod1", "slategray")
names(piecolorvec) <- names(annotationsgrlist)

message("Connecting to biomart")
ensembl <- tryusemart(biomart = "ENSEMBL_MART_ENSEMBL",
    biomartconnection, host = "https://nov2020.archive.ensembl.org",
    alternativeMirror = TRUE)

## Determining proportions on each target category for each query file
message("Determining proportions on each target category for each query file")
numberspielist <- list()
percentageslist <- list()

for (i in seq_len(length(queryfilevec))) {

    ## Processing query files
    queryfile <- queryfilevec[i]
    namequery <- pietitlevec[i]
    outfold <- file.path(outputfolder, namequery)
    if (!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)

    message("\t Processing ", namequery)

    message("\t\t Building GR with query file")
    querygr <- unique(buildgr(queryfile))

    ## Performing overlap on the different categories
    message("\t\t Performing overlap on the different categories")
    res <- performoverlap(annotationsgrlist, querygr)
    overlappriority <- res[[1]]
    overlap <- res[[2]]
    annonamesvec <- res[[3]]
    message("\t\t Plotting piechart")
    res <- performpiechart(annonamesvec, overlappriority, piecolorvec, outfold,
            namequery, percentagerepvec, querygr)
    numberspielist <- c(numberspielist, res[1])
    percentageslist <- c(percentageslist, res[2])
    subjecthitsnamespriority <- res[[3]]

    ## Perform an upset diagram
    performupset(annonamesvec, overlap, namequery, outfold)

    ## Output the gff of querygr per category defined by overlappriority
    peaksidxbycatpriorlist <- savingPeaksPerCategory(overlappriority,
            subjecthitsnamespriority, querygr, outfold)

    ## Output the gff of the promoters
    outputgffprom(annotationsgrlist, querygr, peakidxbycatpriorlist,
            symbolstab, ensembl, outfold)
}

names(numberspielist) <- names(percentageslist) <- pietitlevec

## Perform barplot of all exp
message("Generating barplot for all experiments")
colVec <- c("antiquewhite1", "antiquewhite2", "antiquewhite3", 
        "aquamarine2", "aquamarine3", "aquamarine4", 
        "cadetblue2", "cadetblue3", "cadetblue4", 
        "chocolate1", "chocolate2", "chocolate3", "chocolate4",
        "black")
groupedBarplot(percentageslist, "overlap proportions", outputfolder, percentagerepvec, 
        colVec)
groupedBarplot(numberspielist, "Nb_of_overlap", outputfolder, cntrepeats, 
        colVec)

