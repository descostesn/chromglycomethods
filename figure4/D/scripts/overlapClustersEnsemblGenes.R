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

!!!!!!!!!!!!!!!!!!! NEED TO ADAPT THE SCRIPT



#############
## FUNCTIONS
#############


checkingoutputfolder <- function(output_path) {
    if (!file.exists(output_path))
        dir.create(output_path, recursive = TRUE)
}

checkparams <- function(clusterpathvec, expnamevec, outputfoldervec) {

    if (length(clusterpathvec) != 2 || length(expnamevec) != 2)
        stop("\n This script takes only two exp as input\n")

    extsuffix <- unique(sapply(basename(clusterpathvec),
        function(x) strsplit(x, "\\.")[[1]][2]))

    if (!isTRUE(all.equal(length(extsuffix), 1)) &&
        !isTRUE(all.equal(extsuffix, "gff")))
            stop("The files should be in gff format.")

    if (!isTRUE(all.equal(length(clusterpathvec), length(outputfoldervec))))
        stop("One output folder should be given by cluster file")

    checkingoutputfolder(outputfolder)
}


buildgffrangeslist <- function(clusterpathvec) {
    grlist <- list()

    for (i in seq_len(length(clusterpathvec))) {

        current_gff <- read.delim(clusterpathvec[i], stringsAsFactors = FALSE,
            header = FALSE)

        ## Check for duplicated names
        genenamevec <- current_gff[, 3]
        if (length(genenamevec) != length(unique(genenamevec)))
            genenamevec <- make.unique(genenamevec, sep = "-")

        ## Check for duplicated ranges
        tmp <- paste(current_gff[, 1], current_gff[, 4], current_gff[, 5])
        idxdup <- which(duplicated(tmp))
        if (!isTRUE(all.equal(length(idxdup), 0))) {
            message("Removing ", length(idxdup), "/", nrow(current_gff))
            current_gff <- current_gff[-idxdup, ]
            genenamevec <- genenamevec[-idxdup]
        }

        grlist[[i]] <- GenomicRanges::GRanges(
            seqnames = current_gff[, 1],
            ranges = IRanges::IRanges(start = current_gff[, 4],
                                  end = current_gff[, 5],
                                  names = genenamevec),
            strand = current_gff[, 7])
    }
    return(grlist)
}


buildgrensembl <- function(currentpath) {

    message("\t Processing ", currentpath)

    !!!!!!!!!
    genenamevec <- unlist(lapply(strsplit(unlist(lapply(
                strsplit(current_gff$V9, ";"), "[", 2)), "="), "[", 2))
    if (length(genenamevec) != length(unique(genenamevec)))
            genenamevec <- make.unique(genenamevec, sep = "-")
    !!!!!!!!!!
    fi <- read.table(currentpath, stringsAsFactors = FALSE)
    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = fi$V9),
            strand = fi$V7)
    return(gr)
}



##############
# MAIN
##############

# Retreives the parameters
checkparams(clusterpathvec, expnamevec, outputfolder)

message("Reading gff input and converting to rangedData")
grlist <- buildgffrangeslist(clusterpathvec)
grensemble <- buildgrensembl()
mapply(function(currentgr, currentname, enstab))




message("Performing the overlap")
resoverlap <- ChIPpeakAnno::findOverlapsOfPeaks(grlist[[1]],
                        grlist[[2]])

message("Writting the overlapping peaks")
gff1 <- data.frame(
    seqname = as.character(resoverlap$overlappingPeaks[[1]][, 2]),
    source = "vennDiagram_overlapGFF",
    feature = resoverlap$overlappingPeaks[[1]][, 1],
    start = resoverlap$overlappingPeaks[[1]][, 3],
    end = resoverlap$overlappingPeaks[[1]][, 4],
    score = resoverlap$overlappingPeaks[[1]][, 5],
    strand = as.character(resoverlap$overlappingPeaks[[1]][, 6]),
    frame = ".", group = ".")

gff2 <- data.frame(
    seqname = as.character(resoverlap$overlappingPeaks[[1]][, 8]),
    source = "vennDiagram_overlapGFF",
    feature = resoverlap$overlappingPeaks[[1]][, 7],
    start = resoverlap$overlappingPeaks[[1]][, 9],
    end = resoverlap$overlappingPeaks[[1]][, 10],
    score = resoverlap$overlappingPeaks[[1]][, 11],
    strand = resoverlap$overlappingPeaks[[1]][, 12],
    frame = ".", group = ".")

message("Writing output files")
write.table(gff1,
    file = file.path(outputfolder, paste0(expnamevec[1], ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gff2,
    file = file.path(outputfolder, paste0(ensemblname, ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Done.")
