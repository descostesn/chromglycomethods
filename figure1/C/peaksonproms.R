#################
## This script aims at reprogramming the venn diagrams in a more sophisticated
## manner integrating the repeats' annotations. Some code was copied/pasted
## from ChIPPeakAnno (3.23.3) function assignChromosomeRegion and from ChIPSeeker
## (1.24.0) function upsetPlot.
## Run on R 4.2.0
# Descostes June 2020
#################


# library("ChIPpeakAnno")
# library("RColorBrewer")
# library("ggplot2")
# library("reshape2")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("GenomicFeatures")
library("GenomicRanges")
library("S4Vectors")
library("biomaRt")

################
# PARAMETERS
################


queryfile <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff" # nolint
usealtmirror <- TRUE

#############
## FUNCTIONS
#############

buildgr <- function(currentpath) {
    message("\t Processing ", currentpath)
    fi <- read.table(currentpath, stringsAsFactors = FALSE)
    gr <- GenomicRanges::GRanges(seqnames = fi$V1,
            ranges = IRanges::IRanges(start = fi$V4, end = fi$V5,
                names = fi$V9), strand = fi$V7)
    return(gr)
}

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

## Function by Ilyess Rachedi
trygetbm <- function(attributes, ensembl, values = NULL, filters = NULL) {
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

retrievegeneinfo  <- function(ensembl, annogr) {
    attributes <- c('chromosome_name', 'ensembl_gene_id', 'external_gene_name', # nolint
            'start_position', 'end_position', 'strand', # nolint
            'ensembl_transcript_id_version') # nolint
    symbolstab <- trygetbm(attributes, ensembl, values = names(annogr),
        filters='ensembl_transcript_id_version') # nolint
    symbolstab$strand[which(symbolstab$strand == 1)] <- '+' # nolint
    symbolstab$chromosome_name <- paste0("chr", symbolstab$chromosome_name)
    if (!isTRUE(all.equal(length(which(symbolstab$strand == -1)), 0)))
        symbolstab$strand[which(symbolstab$strand == -1)] <- '-' # nolint
    return(symbolstab)
}



################
# MAIN
################

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## Filtering chromosomes
message("Filtering chromosomes")
seqlevels(txdb) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")

## Retrieve promoters coordinates
promotersgr <- unique(GenomicFeatures::promoters(txdb, upstream = 1000,
                    downstream = 1000))
querygr <- unique(buildgr(queryfile))

## Perform overlap between O-GlcNac peaks and promoters
resultoverlap <- GenomicRanges::findOverlaps(querygr, promotersgr,
    ignore.strand = FALSE)
idxkeep <- which(!duplicated(S4Vectors::queryHits(resultoverlap)))
resultoverlap <- resultOverlap[idxkeep, ]
promotersgr <- promotersgr[S4Vectors::subjectHits(result), ]

## Retrieving information on genes
ensembl <- tryusemart(biomart = "ENSEMBL_MART_ENSEMBL",
    "mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org",
    alternativemirror = usealtmirror)
symbolstab <- retrievegeneinfo(ensembl, promotersgr)
idxtable <- match(names(promotersgr), symbolstab$ensembl_transcript_id_version)
idxna <- which(is.na(idxtable))
if (!isTRUE(all.equal(length(idxna), 0))) {
    promotersgr <- promotersgr[-idxna, ]
    idxtable <- idxtable[-idxna]
}
promgff <- data.frame(seqname = symbolstab$chromosome_name[idxtable],
            source = symbolstab$external_gene_name[idxtable],
            feature = symbolstab$ensembl_gene_id[idxtable],
            start = symbolstab$start_position[idxtable],
            end = symbolstab$end_position[idxtable],
            score = 0, strand = symbolstab$strand[idxtable], frame = ".",
            group = ".")
promgff <- promgff[-which(duplicated(promgff$feature)), ]

