##############
## This script computes the union of intervals from several bed/gff files and
## output a bed file.
## Descostes
##############


library("IRanges")
library("GenomicRanges")


################
# PARAMETERS
################

filepathvec <- c("../data/DLD1GlcNAcNoDoxAux_rep1.gff",
"../data/DLD1GlcNAcNoDoxAux_rep2.gff",
"../data/DLD1GlcNAcDoxAux_rep1_peaks.gff",
"../data/DLD1GlcNAcDoxAux_rep2_peaks.gff")

outputfile <- "../results/union_OGlcNac_noauxaux-fig4C.bed"
inputformat <- "gff"
featurename <- "union-glc-noauxaux"


##############
# MAIN
##############

if (!isTRUE(all.equal(inputformat, "bed")) &&
    !isTRUE(all.equal(inputformat, "gff")))
    stop("input format should be bed or gff")

if (!file.exists(dirname(outputfile)))
    dir.create(dirname(outputfile), recursive = TRUE)

message("Reading peak files")
peaklist <- lapply(filepathvec, read.table, stringsAsFactors = FALSE,
    header = FALSE)

peaksdf <- do.call(rbind, peaklist)

if (isTRUE(all.equal(inputformat, "bed"))) {
    peaksgr <- GenomicRanges::GRanges(seqnames = peaksdf[, 1],
        ranges = IRanges::IRanges(start = peaksdf[, 2], end = peaksdf[, 3]))
} else {
    peaksgr <- GenomicRanges::GRanges(seqnames = peaksdf[,1],
        ranges = IRanges::IRanges(start = peaksdf[, 4], end = peaksdf[, 5]))
}

message("Reducing intervals")
peaksgr <- GenomicRanges::reduce(peaksgr)

res <- data.frame(chrom = as.character(seqnames(peaksgr)),
    start = start(peaksgr),
    end = end(peaksgr),
    name = featurename,  score = 0, strand = '+',
    thickStart = start(peaksgr), thickEnd = end(peaksgr),
    itemrgb = "255,0,0", blockCount = 1,
    blockSizes = end(peaksgr) - start(peaksgr), blockStarts = start(peaksgr))

message("The union returned ", nrow(res), " peaks")
message("Writing ", outputfile)
write.table(res, file = outputfile, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE)
