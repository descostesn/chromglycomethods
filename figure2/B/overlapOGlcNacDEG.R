###############
# This script aimns at representing differing and common intervals for
# 3 gff files.
#
# Descostes - R-4.2.0
###############


library("ChIPpeakAnno")
library("VennDiagram")


#################
# PARAMETERS
################


gff_file_vec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHGPeaksRep1_vs_Rep2/Glc1.gff", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_down-ensembl.gff", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_up-ensembl.gff")  # nolint
outfolder <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHG_vs_DEGsiogt/mergedrep1-2/test" # nolint
comparisonname <- "mergedrep1-2_vs_DEGsiogt"
expnamevec <- c("mergedrep1-2", "Down", "Up")
genome_version <- "mm10"
colvec <- c("#E69F00", "#56B4E9", "#E95680")
outformat <- "png"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gff_file_vec <- c("results/Glc1.gff", "results/log0_siogtdown-ensembl.gff",
"results/log0_siogtup-ensembl.gff")
outfolder <- "results"
comparisonname <- "glcrep1-2_vs_degsiogt"
expnamevec <- c("glcrep1-2", "Down", "Up")
colvec <- c("#E69F00", "#56B4E9", "#E95680")
outformat <- "png"


################


################
# FUNCTION
################

retrieve_peaks_number <- function(peak_list, label_name) {
    peak_number <- lapply(peak_list, length)[[label_name]]
    if (is.null(peak_number)) return(0) else return(peak_number)
}

checkparams <- function(gff_file_vec, expnamevec, outformat) {

    if (length(gff_file_vec) > 5)
        stop("This script takes at most 5 gff files as input")

    if (length(gff_file_vec) < 2)
        stop("This script takes a minimum of 2 gff files as input")

    if (!isTRUE(all.equal(length(gff_file_vec), length(expnamevec))))
        stop("The number of exp names should be equal to the number of gff ",
            "files")

    if (!isTRUE(all.equal(outformat, "ps")) &&
        !isTRUE(all.equal(outformat, "png")) &&
        !isTRUE(all.equal(outformat, "pdf")))
        stop("outformat should be png, pdf or ps")
}

convert2gr <- function(gfflist) {
    return(lapply(gfflist, function(currentgff) {
    return(GenomicRanges::GRanges(seqnames = currentgff[, 1],
        ranges = IRanges::IRanges(start = currentgff[, 4],
            end = currentgff[, 5]), names = currentgff[, 2],
            strand = currentgff[, 7]))
    }))
}

readgff <- function(gff_file_vec) {
    return(lapply(gff_file_vec, function(currentpath) {
    currentgff <- read.delim(currentpath, header = FALSE)

    ## Removing duplicates
    semicolon_feature <- paste(currentgff[, 1], currentgff[, 4],
        currentgff[, 5], sep = ";")
    idxdup <- which(duplicated(semicolon_feature))
    if (!isTRUE(all.equal(length(idxdup), 0)))
        currentgff <- currentgff[-idxdup, ]

    ## Removing any chrNA or chrMT
    idxchrna <-  which(currentgff$V1 == "chrNA" | currentgff$V1 == "chrMT")
    if (!isTRUE(all.equal(length(idxchrna), 0)))
        currentgff <- currentgff[-idxchrna, ]

    return(currentgff)}))
}

.defineformat <- function(outformat, outfolder, comparisonname) {
    if (outformat == "png") {
        png(filename = file.path(outfolder, paste0(comparisonname, ".png")),
            width = 1000, height = 1000, bg = "transparent")
    }else {
        pdf(file = file.path(outfolder, paste0(comparisonname, ".pdf")),
            width = 10, height = 10)
    }
}

eulerthree <- function(ol, outformat, outfolder, comparisonname, expnamevec,
    colvec) {
    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1")
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2")
    area3 <- retrieve_peaks_number(ol$peaklist, "peaks3")
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2")
    area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3")
    area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3")
    area1_area2_area3 <- retrieve_peaks_number(ol$peaklist,
        "peaks1///peaks2///peaks3")

    .defineformat(outformat, outfolder, comparisonname)
    VennDiagram::draw.triple.venn(
        area1 = area1 + area1_area2 + area1_area3 + area1_area2_area3,
        area2 = area2 + area1_area2 + area2_area3 + area1_area2_area3,
        area3 = area3 + area1_area3 + area2_area3 + area1_area2_area3,
        n12 = area1_area2 + area1_area2_area3,
        n23 = area2_area3 + area1_area2_area3,
        n13 = area1_area3 + area1_area2_area3,
        n123 = area1_area2_area3, category = expnamevec, euler.d = TRUE,
        scaled = TRUE, ext.text = FALSE, col = colvec, fill = colvec,
        cex = rep(2, 7), cat.cex = rep(2.5, 3))
    dev.off()
}


################



##############
# MAIN
##############

checkparams(gff_file_vec, expnamevec, outformat)
if (!file.exists(outfolder))
    dir.create(outfolder, recursive = TRUE)

message("Reading GFF files and converting to GRanges")
gfflist <- readgff(gff_file_vec)
gffgrlist <- convert2gr(gfflist)

message("Performing the overlap")
peaks1 <- gffgrlist[[1]]
peaks2 <- gffgrlist[[2]]
peaks3 <- gffgrlist[[3]]
ol <- ChIPpeakAnno::findOverlapsOfPeaks(peaks1, peaks2, peaks3, maxgap = 0,
        connectedPeaks = "keepAll")

# Making the venn diagram venneuler
message("Generating the venn diagram")
eulerthree(ol, outformat, outfolder, comparisonname, expnamevec, colvec)

## Writing list of element per overlap
outfoldpeaks <- file.path(outfolder, "peaks_per_circle/")
if (!file.exists(outfoldpeaks))
    dir.create(outfoldpeaks, recursive = TRUE)

message("Writing list of element per overlap to ", outfoldpeaks, "/")
invisible(mapply(function(peaktab, peakname, expnamevec, outfoldpeaks) {

    if (nchar(peakname) > 6) {
        idx <- as.numeric(unlist(lapply(strsplit(unlist(
            strsplit(peakname, "///")), "peaks"), "[", 2)))
        outfile <- paste(expnamevec[idx], collapse = "_")
    } else {
        idx <- as.numeric(unlist(lapply(strsplit(peakname, "peaks"), "[", 2)))
        outfile <- expnamevec[idx]
    }
    message("\t", outfile)
    outfile <- file.path(outfoldpeaks, paste0(outfile, ".gff"))
    gffres <- data.frame(seqname = as.character(seqnames(peaktab)),
            source = "vennDiagram_overlapGFF",
            feature = as.character(unlist(lapply(
                elementMetadata(peaktab)$peakNames, paste, collapse = "-"))),
            start = start(ranges(peaktab)), end = end(ranges(peaktab)),
            score = 0, strand = as.character(strand(peaktab)),
            frame = ".", group = ".")
    write.table(gffres, file = outfile, sep = "\t", quote = FALSE,
        col.names = FALSE, row.names = FALSE)
}, ol$peaklist, names(ol$peaklist), MoreArgs = list(expnamevec, outfoldpeaks)))
