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


gfffilevec <- c("ESCHGGlcNAc_rep1.gff",
"ESCHGGlcNAc_rep2.gff")
outputfolder <- "results"
expnamevec <- c("Glc1", "Glc2")


#############
## FUNCTIONS
#############


checkingoutputfolder <- function(output_path) {
    if (!file.exists(output_path))
        dir.create(output_path, recursive = TRUE)
}

checkparams <- function(gfffilevec, expnamevec, outputfolder) {

    if (length(gfffilevec) != 2 || length(expnamevec) != 2)
        stop("\n This script takes only two exp as input\n")

    extsuffix <- unique(sapply(basename(gfffilevec),
        function(x) strsplit(x, "\\.")[[1]][2]))

    if (!isTRUE(all.equal(length(extsuffix), 1)) &&
        !isTRUE(all.equal(extsuffix, "gff")))
            stop("The files should be in gff format.")

    checkingoutputfolder(outputfolder)
}


buildgffrangeslist <- function(gfffilevec, extractensemblgenename = FALSE) {
    gffgrangeslist <- list()

    for (i in seq_len(length(gfffilevec))) {

        currentgff <- read.delim(gfffilevec[i], stringsAsFactors = FALSE,
            header = FALSE)

        ## In case the overlap is done with ensembl genes, extract gene name
        if (isTRUE(all.equal(i, 2)) && extractensemblgenename)
            genenamevec <- unlist(lapply(strsplit(unlist(lapply(
                strsplit(currentgff$V9, ";"), "[", 2)), "="), "[", 2))
        else
            genenamevec <- currentgff[, 3]

        ## Check for duplicated names
        if (length(genenamevec) != length(unique(genenamevec)))
            genenamevec <- make.unique(genenamevec, sep = "-")

        ## Check for duplicated ranges
        tmp <- paste(currentgff[, 1], currentgff[, 4], currentgff[, 5])
        idxdup <- which(duplicated(tmp))
        if (!isTRUE(all.equal(length(idxdup), 0))) {
            message("Removing ", length(idxdup), "/", nrow(currentgff))
            currentgff <- currentgff[-idxdup, ]
            genenamevec <- genenamevec[-idxdup]
        }

        gffgrangeslist[[i]] <- GenomicRanges::GRanges(
            seqnames = currentgff[, 1],
            ranges = IRanges::IRanges(start = currentgff[, 4],
                                  end = currentgff[, 5],
                                  names = genenamevec),
            strand = currentgff[, 7])
    }
    return(gffgrangeslist)
}



##############
# MAIN
##############

checkparams(gfffilevec, expnamevec, outputfolder)

message("Reading gff input and converting to rangedData")
gffgrangeslist <- buildgffrangeslist(gfffilevec)

message(expnamevec[1], " has ", length(gffgrangeslist[[1]]), " peaks and ",
expnamevec[2], " has ", length(gffgrangeslist[[2]]), " peaks")

message("Performing the overlap")
resultoverlap <- ChIPpeakAnno::findOverlapsOfPeaks(gffgrangeslist[[1]],
                        gffgrangeslist[[2]])

res <-  resultoverlap$overlappingPeaks[[1]]
message("The intersection gives ", nrow(res), " peaks.")

gffpeak1 <- data.frame(
    seqname = as.character(res[, 2]),
    source = "vennDiagram_overlapGFF",
    feature = res[, 1],
    start = res[, 3],
    end = res[, 4],
    score = res[, 5],
    strand = as.character(res[, 6]),
    frame = ".", group = ".")

gffpeak2 <- data.frame(
    seqname = as.character(res[, 8]),
    source = "vennDiagram_overlapGFF",
    feature = res[, 7],
    start = res[, 9],
    end = res[, 10],
    score = res[, 11],
    strand = res[, 12],
    frame = ".", group = ".")

message("Writing the coordinates of intersecting peaks for each replicate to ",
outputfolder, "/")
write.table(gffpeak1,
    file = file.path(outputfolder, paste0(expnamevec[1], ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gffpeak2,
    file = file.path(outputfolder, paste0(expnamevec[2], ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Done.")
