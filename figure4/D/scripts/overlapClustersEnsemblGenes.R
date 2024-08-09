#################
# This script performs the overlap of two exp given as gff and output the
# individual overlapping peaks.
#
# Descostes May 2016 - update Nov 2023 - R-4.2.0
#################


library("GenomicRanges")
library("ChIPpeakAnno")


################
# PARAMETERS
################


clusterpathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/5groups/test/cluster_1-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/5groups/test/cluster_2-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/5groups/test/cluster_3-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/5groups/test/cluster_4-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/5groups/test/cluster_5-compartmentsgff/promoters.gff") # nolint
expnamevec <- paste0("cluster", seq_len(5))
ensemblpath <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/ensembl/Homo_sapiens_GRCh38110_chrfiltered_genes.gff" # nolint
ensemblname <- "Homo_sapiens_GRCh38110genes"
outputfoldervec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_1-compartmentsgff/promoters_vs_ensemblgenes/test", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_2-compartmentsgff/promoters_vs_ensemblgenes/test", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_3-compartmentsgff/promoters_vs_ensemblgenes/test", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_4-compartmentsgff/promoters_vs_ensemblgenes/test", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_5-compartmentsgff/promoters_vs_ensemblgenes/test") # nolint
expname <- "promspeaks_ensemblgenes"



#############
## FUNCTIONS
#############


checkingoutputfolder <- function(output_path) {
    if (!file.exists(output_path))
        dir.create(output_path, recursive = TRUE)
}

checkparams <- function(clusterpathvec, expnamevec, outputfoldervec) {

    if (!isTRUE(all.equal(length(clusterpathvec), length(expnamevec))))
        stop("One name should be given to each cluster") # nolint

    extsuffix <- unique(sapply(basename(clusterpathvec),
        function(x) strsplit(x, "\\.")[[1]][2]))

    if (!isTRUE(all.equal(length(extsuffix), 1)) &&
        !isTRUE(all.equal(extsuffix, "gff")))
            stop("The files should be in gff format.") # nolint

    if (!isTRUE(all.equal(length(clusterpathvec), length(outputfoldervec))))
        stop("One output folder should be given by cluster file") # nolint

    invisible(lapply(outputfoldervec, checkingoutputfolder))
}


.checkduplicatedranges <- function(currentgff) {

    ## Check for duplicated ranges
    tmp <- paste(currentgff[, 1], currentgff[, 4], currentgff[, 5], sep="-")
    idxdup <- which(duplicated(tmp))
    if (!isTRUE(all.equal(length(idxdup), 0))) {
        message("\t\t Removing ", length(idxdup), "/", nrow(currentgff), # nolint
            " duplicated ranges")
        currentgff <- currentgff[-idxdup, ]
    }
    return(currentgff)
}

buildgffrangeslist <- function(clusterpathvec) {
    grlist <- list()

    for (i in seq_len(length(clusterpathvec))) {

        currentgff <- read.delim(clusterpathvec[i], stringsAsFactors = FALSE,
            header = FALSE)
        currentgff <- .checkduplicatedranges(currentgff)

        ## Check for duplicated names
        genenamevec <- currentgff[, 3]
        if (length(genenamevec) != length(unique(genenamevec)))
            genenamevec <- make.unique(genenamevec, sep = "-")

        ## Building genomicranges object
        grlist[[i]] <- GenomicRanges::GRanges(
            seqnames = currentgff[, 1],
            ranges = IRanges::IRanges(start = currentgff[, 4],
                                  end = currentgff[, 5],
                                  names = genenamevec),
            strand = currentgff[, 7])
    }
    return(grlist)
}


buildgrensembl <- function(currentpath) {

    fi <- read.delim(currentpath, stringsAsFactors = FALSE, header = FALSE)
    fi <- .checkduplicatedranges(fi)
    genenamevec <- unlist(lapply(strsplit(unlist(lapply(
                strsplit(fi$V9, ";"), "[", 2)), "="), "[", 2))
    if (length(genenamevec) != length(unique(genenamevec)))
            genenamevec <- make.unique(genenamevec, sep = "-")
    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = genenamevec),
            strand = fi$V7)
    return(gr)
}

removeduplicatedens <- function(gff2) {
    featuregff2 <- gff2$feature
    if (!isTRUE(all.equal(length(featuregff2), length(unique(featuregff2))))) {
        idxdup <- which(duplicated(featuregff2))
        message("\t\t Removing ", length(idxdup), "/", nrow(gff2), # nolint
            " duplicated ensembl genes.") # nolint
        gff2 <- gff2[-idxdup, ]
    }
    return(gff2)
}


##############
# MAIN
##############

# Retreives the parameters
checkparams(clusterpathvec, expnamevec, outputfoldervec)

message("Reading gff input and converting to genomicranges Data") # nolint
grlist <- buildgffrangeslist(clusterpathvec)
names(grlist) <- expnamevec
grensembl <- buildgrensembl(ensemblpath)

message("Performing overlap of each cluster with ensembl annotations")  # nolint
# currentgr=grlist[[1]]
# currentname=expnamevec[1]
# outfold=outputfoldervec[1]
# grens=grensembl
reslist <- mapply(function(currentgr, currentname, outfold, grens,
    ensemblname) {

    message("\t Performing overlap for ", currentname)
    resoverlap <- ChIPpeakAnno::findOverlapsOfPeaks(currentgr, grens)

    message("\t Converting result to gff format")
    gff1 <- data.frame(
    seqname = as.character(resoverlap$overlappingPeaks[[1]][, 2]),
    source = "vennDiagram_overlapGFF",
    feature = resoverlap$overlappingPeaks[[1]][, 1],
    start = resoverlap$overlappingPeaks[[1]][, 3],
    end = resoverlap$overlappingPeaks[[1]][, 4],
    score = resoverlap$overlappingPeaks[[1]][, 5],
    strand = as.character(resoverlap$overlappingPeaks[[1]][, 6]),
    frame = ".", group = ".")
    gff1$feature <-  gsub("-[0-9]+", "",gff1$feature)

    gff2 <- data.frame(
    seqname = as.character(resoverlap$overlappingPeaks[[1]][, 8]),
    source = "vennDiagram_overlapGFF",
    feature = resoverlap$overlappingPeaks[[1]][, 7],
    start = resoverlap$overlappingPeaks[[1]][, 9],
    end = resoverlap$overlappingPeaks[[1]][, 10],
    score = resoverlap$overlappingPeaks[[1]][, 11],
    strand = resoverlap$overlappingPeaks[[1]][, 12],
    frame = ".", group = ".")
    gff2 <- removeduplicatedens(gff2)

    message("\t The number of genes for ", currentname, " is ", nrow(gff2))
    message("\t Writing coordinates to ", outfold)
    write.table(gff1,
    file = file.path(outfold, paste0(currentname, ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(gff2,
    file = file.path(outfold, paste0(ensemblname, "-", currentname,
        ".gff")), sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = FALSE)

}, grlist, expnamevec, outputfoldervec,
    MoreArgs = list(grensembl, ensemblname), SIMPLIFY = FALSE)

