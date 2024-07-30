##############
## This script computes the union of intervals from several bed/gff files and
## output a gff file.
## Descostes
##############


library("IRanges")
library("GenomicRanges")


################
# PARAMETERS
################


filepathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc2_lane1sample13_peaks_narrowPeak.gff")
outputfile <- "/g/romebioinfo/tmp/union_sept2023mouse_HG1-2.gff"
inputformat <- "gff"
featurename <- "union-glc-rep1_2"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

filepathvec <- c("data/ESCHGGlcNAc_rep1_peaks.gff",
    "data/ESCHGGlcNAc_rep2_peaks.gff")
outputfile <- "result/union_sept2023mouse_HG1-2.gff"
inputformat <- "gff"
featurename <- "union-glc-rep1_2"

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

peaksdf <- do.call(rbind, files_vec)

if (isTRUE(all.equal(inputformat, "bed"))) {
    peaksgr <- GenomicRanges::GRanges(seqnames = peaksdf[, 1],
        ranges = IRanges::IRanges(start = peaksdf[, 2], end = peaksdf[, 3]))
} else {
    peaksgr <- GenomicRanges::GRanges(seqnames = peaksdf[,1],
        ranges = IRanges::IRanges(start = peaksdf[, 4], end = peaksdf[, 5]))
}

message("Reducing intervals")
peaksgr <- GenomicRanges::reduce(peaksgr)

res <- data.frame(seqname = as.character(seqnames(peaksgr)),
    source = "union of peaks", feature = featurename, start = start(peaksgr),
    end = end(peaksgr), score = 0, strand = '+', frame = ".",
    group = make.unique(rep("group", length(start(peaksgr))), sep = "-"))

write.table(res, file = outputfile, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE)

