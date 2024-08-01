###########################
# The deseq2 results are used to generate a volcano plot using a padj < 0.05.
# The venn diagrams list between glc peaks and down-reg genes is used to
# highlight those genes.
#
# Descostes - R 4.2.0 (R/4.2.0-foss-2021b) - June 2024
###########################

library("ggplot2")
library("clusterProfiler")
library("ggrepel")

##################
# PARAMETERS
##################


deseq2paths <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/deseq2/DESeq2_result_file_on_data_9_data_8_and_others.tabular") # nolint
overdownpeakspath <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHGmergedRep1-2_vs_downsiogt/siogtdown.gff" # nolint
overuppeakspath <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHGmergedRep1-2_vs_upsiogt/siogtup.gff" # nolint
outfold <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/volcanoCustom/siogt" # nolint

!!!!!!!!!!!!!!!!!


deseq2paths <- c("data/resultDeseq2_siogt.txt")
overdownpeakspath <- "data/siogtdown_withOGlcNac.gff"
overuppeakspath <- "data/siogtup_withOGlcNac.gff"
outfold <- "results/"

##################
#FUNCTIONS
##################

verifyna <- function(idx) {
    idxna <- which(is.na(idx))
    lna <- length(idxna)
    if (!isTRUE(all.equal(lna, 0)))
        stop("Some ensembl ids were not retrieved")
}

convertsymbols <- function(gfftab, colids) {
    res <- clusterProfiler::bitr(gfftab[, colids], fromType = "ENSEMBL",
                    toType = "SYMBOL", OrgDb = "org.Mm.eg.db", drop = FALSE)
    return(res)
}

readfilteraddsymbols <- function(filepath, colnamevec, colids,
    checkna = FALSE) {
    filetab <- read.delim(filepath, header = FALSE)
    idxdup <- which(duplicated(filetab[, colids]))
    if (!isTRUE(all.equal(length(idxdup), 0)))
        filetab <- filetab[-idxdup, ]
    symbolvec <- convertsymbols(filetab, colids)
    idxdupsym <- which(duplicated(symbolvec[,1]))
    if (!isTRUE(all.equal(length(idxdupsym), 0)))
        symbolvec <- symbolvec[-idxdupsym, ]
    if (checkna)
        if (!isTRUE(all.equal(length(which(is.na(symbolvec[, 2]))), 0)))
            stop("Symbols were not found")
    filetab <- cbind(filetab, symbolvec[, 2])
    colnames(filetab) <- colnamevec
    return(filetab)
}



##################
# MAIN
##################

if (!file.exists(outfold))
    dir.create(outfold, recursive = TRUE)

## Reading files
message("Reading files")
colnamedeseq2 <- c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", # nolint
            "pvalue", "padj", "symbol")
deseq2 <- suppressWarnings(readfilteraddsymbols(deseq2paths, colnamedeseq2, 1))

colnamespeaks <- c("chromosome_name", "source", "ensembl_gene_id",
    "start_position", "end_position", "length", "strand", "group", "frame",
    "symbol")
downpeaks <- readfilteraddsymbols(overdownpeakspath, colnamespeaks, 3,
    checkna = TRUE)
uppeaks <- readfilteraddsymbols(overuppeakspath, colnamespeaks, 3,
    checkna = TRUE)

## Add a status column to deseq2tab: up, down, glcpeaks, nodiff
message("Defining status")
statusvec <- rep("nodiff", nrow(deseq2))
idxdown <- which(deseq2$log2FoldChange < 0 & deseq2$padj < 0.05)
idxup <- which(deseq2$log2FoldChange > 0 & deseq2$padj < 0.05)
idxglcdown <- match(downpeaks$ensembl_gene_id, deseq2$geneID)
idxglcdown <- intersect(idxdown, idxglcdown)
idxglcup <- match(uppeaks$ensembl_gene_id, deseq2$geneID)
idxglcup <- intersect(idxup, idxglcup)
verifyna(idxglcdown)
verifyna(idxglcup)
statusvec[idxdown] <- "down"
statusvec[idxup] <- "up"
statusvec[idxglcdown] <- "downwithglc"
statusvec[idxglcup] <- "upwithglc"
deseq2 <- cbind(deseq2, statusvec)

message("The number of genes per category is: ")
print(table(deseq2$statusvec))


## Plotting volcano
message("Plotting volcano")
colvec <- c("deepskyblue4", "darkturquoise", "darkgrey", "orange4",
    "orange")
idxnoglc <- which(deseq2$statusvec != "downwithglc" &
    deseq2$statusvec != "upwithglc")
p <- ggplot(data = deseq2[idxnoglc, ],
    aes(x = log2FoldChange, y = -log10(padj), col = statusvec)) + geom_point() +
    theme_minimal() + scale_color_manual(values = colvec)
p1 <- p + geom_point(data = deseq2[which(deseq2$statusvec == "downwithglc"), ],
        aes(x = log2FoldChange, y = -log10(padj), col = statusvec)) +
        geom_point(data = deseq2[which(deseq2$statusvec == "upwithglc"), ],
        aes(x = log2FoldChange, y = -log10(padj), col = statusvec)) +
        scale_color_manual(values = colvec)

## Adding labels for down-reg genes with glc
deseq2downglc <- deseq2[which(deseq2$statusvec == "downwithglc"), ]
deseq2upglc <- deseq2[which(deseq2$statusvec == "upwithglc"), ]
p2 <- p1 + geom_label_repel(data = deseq2downglc, aes(label = symbol),
                   force = 2, nudge_y = 1,  max.overlaps = Inf,
                   show.legend = FALSE, min.segment.length = 0,
                   direction = "y") +
            geom_label_repel(data = deseq2upglc, aes(label = symbol),
                   force = 2, nudge_y = 1,  max.overlaps = Inf,
                   show.legend = FALSE, min.segment.length = 0,
                   direction = "y")
p2bis <- p1 + geom_label_repel(data = deseq2downglc, aes(label = symbol),
                   force = 2, nudge_y = 1,  max.overlaps = Inf,
                   show.legend = FALSE, min.segment.length = 0) +
              geom_label_repel(data = deseq2upglc, aes(label = symbol),
                   force = 2, nudge_y = 1,  max.overlaps = Inf,
                   show.legend = FALSE, min.segment.length = 0)


ggsave(filename = "volcano-siogt.svg", plot = p1, device = svg(),
     path = outfold)
ggsave(filename = "volcano-siogt-glcdown-vertical.svg", plot = p2,
    device = svg(), path = outfold)
ggsave(filename = "volcano-siogt-glcdown.svg", plot = p2bis, device = svg(),
     path = outfold)