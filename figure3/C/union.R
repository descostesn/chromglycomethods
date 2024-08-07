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

filepathvec <- c("data/ESC_1b_atac_NoDox_rep1_peaks.bed",
"data/ESC_1b_atac_NoDox_rep2_peaks.bed",
"data/ESC_1b_atac_Dox_rep3_peaks.bed",
"data/ESC_1b_atac_Dox_rep1_peaks.bed",
"data/ESC_1b_atac_Dox_rep2_peaks.bed",
"data/ESC_1b_atac_NoDox_rep3_peaks.bed")
outputfile <- "results/ESC1bNoDox_vs_Dox.gtf"
inputformat <- "bed"
analysisname <- "unionESCatac"
featurename <- "peak"

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

groupvec <- paste0("group-", seq_len(length(peaksgr)))
res <- data.frame(chrom = as.character(seqnames(peaksgr)),
    name = analysisname, feature = featurename,
    start = start(peaksgr), end = end(peaksgr),
    score = 0, strand = '+', frame = ".",
    group = paste0("peak_id \"group-", groupvec, "\"; peak_name \"",
    groupvec, "\";"))

message("The union returned ", nrow(res), " peaks")
message("Writing ", outputfile)
write.table(res, file = outputfile, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE)
