#################
# This script measure the expression of genes for several categories:
# 1) Genes having a glc peak at promoter
# 2) Genes not having a glc peak at promoter (same nb of prom than previous)
# 3) Randomly selected genes
# 4) All genes
#
# The comparison is done in the form of a violin plot.
#
# > sessionInfo()
# MatrixGenerics_1.10.0, httr_1.4.2, bit64_4.0.5, assertthat_0.2.1,
# BiocFileCache_2.6.1, blob_1.2.3, GenomeInfoDbData_1.2.9, Rsamtools_2.14.0,
# yaml_2.3.5, progress_1.2.2, pillar_1.7.0, RSQLite_2.2.12, lattice_0.20-45,
# glue_1.6.2, digest_0.6.29, XVector_0.38.0, colorspace_2.0-3, Matrix_1.4-1,
# XML_3.99-0.9, pkgconfig_2.0.3, zlibbioc_1.44.0, purrr_0.3.4, scales_1.2.0,
# BiocParallel_1.32.6, tibble_3.1.6, KEGGREST_1.38.0, farver_2.1.0,
# generics_0.1.2, ellipsis_0.3.2, cachem_1.0.6, withr_2.5.0,
# SummarizedExperiment_1.28.0, cli_3.6.1, magrittr_2.0.3, crayon_1.5.1,
# memoise_2.0.1, fansi_1.0.3, xml2_1.3.3, tools_4.2.0, prettyunits_1.1.1,
# hms_1.1.1, BiocIO_1.8.0, lifecycle_1.0.3, matrixStats_0.62.0, stringr_1.4.0,
# munsell_0.5.0, DelayedArray_0.24.0, Biostrings_2.66.0, compiler_4.2.0,
# rlang_1.1.0, grid_4.2.0, RCurl_1.98-1.6, rjson_0.2.21, rappdirs_0.3.3,
# labeling_0.4.2, bitops_1.0-7, restfulr_0.0.15, gtable_0.3.0, codetools_0.2-18,
# DBI_1.1.2, curl_4.3.2, R6_2.5.1, GenomicAlignments_1.34.1, dplyr_1.0.8,
# rtracklayer_1.58.0, fastmap_1.1.0, bit_4.0.4, utf8_1.2.2, filelock_1.0.2,
# stringi_1.7.6, parallel_4.2.0, Rcpp_1.0.8.3, vctrs_0.6.1, png_0.1-7,
# dbplyr_2.1.1, tidyselect_1.1.2
#
# Descostes - April 2024 with R-4.3.3
#################


library("ggplot2")
library("RColorBrewer")

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("GenomicFeatures")
library("GenomicRanges")
library("S4Vectors")
library("biomaRt")


##################
# PARAMETERS
##################


queryfile <- "data/ESCHGGlcNAc_rep1.gff"
countspath <- "data/ESCRNAseq_SRR11294181counts.txt"
namescountsvec <- "rep1"
countslength <- "data/ESCRNAseq_SRR11294181countslength.txt"
outputfolder <- "./"

## Would an alternative mirror be used to connect to ensembl
usealtmirror <- TRUE
## Number of bp upstream and downstream TSS for the promoter coordinates
upstreambp <- 1000
downstreambp <- 1000



##################
#FUNCTIONS
##################

nonZeroIndex <- function(idx, step) { # nolint: object_name_linter.
        if (isTRUE(all.equal(length(idx), 0)))
            stop("Problem in retrieving expressed genes", step)
}

zeroIndex <- function(idx, idxna, step) { # nolint: object_name_linter.
    if (!isTRUE(all.equal(length(idxna), 0))) {
        msg <- paste0("Problem in retrieving idx, there should not be NA: ",
            step, " - ", length(idxna), " ids not found")
        warning(msg)
        message(msg)
        idx <- idx[-idxna]
    }
    return(idx)
}

defineNoGlcProm <- function(promglc, allcounts) { # nolint
    idx <- match(promglc[, 3], allcounts[, 1])
    idxna <- which(is.na(idx))
    idx <- zeroIndex(idx, idxna, "defineNoGlcProm")

    nonglcidx <- seq_len(nrow(allcounts))[-idx]
    allcountsnoglc <- allcounts[nonglcidx, ]
    allcountsnoglc <- allcountsnoglc[which(allcountsnoglc[, 2] != 0), ]

    if (nrow(allcountsnoglc) < nrow(promglc))
        stop("Problem in the script. Non glc lower than glc.")

    promnoglc <- allcountsnoglc[sample(seq_len(nrow(allcountsnoglc)),
    nrow(promglc)), ]
    return(promnoglc)
}

defineGlcProm <- function(promglc, allcounts) { # nolint
    idx <- match(promglc[, 3], allcounts[, 1])
    idxna <- which(is.na(idx))
    idx <- zeroIndex(idx, idxna, "defineGlcProm")
    return(allcounts[idx, ])
}

prepareMatrix <- function(counts_glcprom, counts_noglcprom, counts_random, # nolint
    allcountsnozero) {

        df_matrix <- data.frame(TPM = c(counts_glcprom$TPM,
        counts_noglcprom$TPM,
        counts_random$TPM,
        allcountsnozero$TPM),
        category <- c(rep("glcprom", nrow(counts_glcprom)),
        rep("noglcprom", nrow(counts_noglcprom)),
        rep("randomprom", nrow(counts_random)),
        rep("allprom", nrow(allcountsnozero))))

        labelvec <- unique(df_matrix$category)
        df_matrix$category <- factor(df_matrix$category, levels = labelvec)

        return(list(df_matrix, category, labelvec))
}

performPlot <- function(labelvec, df_matrix, currentname, outfold,
    showplot = FALSE) {

    colvec <- brewer.pal(n = length(labelvec), name = "Paired")
    maintitle <- "gene expression"

    g <- ggplot(df_matrix, aes(category, asinh(TPM), fill = factor(category),
                            colour = factor(category))) +
            geom_violin(trim = FALSE) + scale_colour_manual(values = colvec) +
            scale_fill_manual(values = rep("white", length(colvec))) +
            ggtitle(maintitle) + theme_classic()
    g1 <- g + theme(legend.position = "none",
            panel.background = element_rect(fill = 'white',
            colour = 'black'),
            plot.title = element_text(size = 20, face = "bold",
            margin = margin(10, 0, 10, 0)),
            axis.text.x = element_text(angle = 90, size = 10))
    g2 <- g1 + geom_boxplot(width = 0.1, outlier.colour = "NA")

    if (showplot) {
        print(g2)
    } else {
        message("Writing file")
        filename <- paste0("violingplotExpression-", currentname, ".pdf")
        ggsave(filename, plot = g2, device = "pdf", path = outfold)
    }
}

wilcoxonTest <- function(counts_glcprom, counts_noglcprom, counts_random, # nolint
    allcountsnozero, labelvec, currentname, outfold){

    message("Computing mann-whitney-wilcoxson on each violin")
    tpmlist <- list(asinh(counts_glcprom$TPM), asinh(counts_noglcprom$TPM),
     asinh(counts_random$TPM), asinh(allcountsnozero$TPM))
    all_combination_vec <- combn(seq_len(length(tpmlist)), 2)
    result_wilcox_list <- apply(all_combination_vec, MARGIN=2, function(x) {
                        result <- wilcox.test(tpmlist[[x[1]]], tpmlist[[x[2]]],
                        mu = 0, paired = FALSE, conf.int = TRUE)
                        return(result)
                    })
    namescompvec <- apply(all_combination_vec, MARGIN = 2, function(x) {
        paste(labelvec[x[[1]]], labelvec[x[[2]]], sep = "-")})
    names(result_wilcox_list) <- namescompvec

    message("Writing file")
    invisible(mapply(function(test, compname, repname, outfold) {
        filename <- paste("stats", repname, compname, ".txt", sep = "_")
        write(unlist(test), file=file.path(outfold, filename), ncolumns = 1)
    }, result_wilcox_list, names(result_wilcox_list),
    MoreArgs = list(currentname, outfold)))
}

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


retrieveglcproms <- function(upstreambp, downstreambp, queryfile,
    usealtmirror) {

    message("Retrieve the promoters coordinates")
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene # nolint
    ## Filtering chromosomes
    message("Filtering chromosomes")
    seqlevels(txdb) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
    ## Retrieve promoters coordinates
    promotersgr <- unique(GenomicFeatures::promoters(txdb,
        upstream = upstreambp, downstream = downstreambp))

    ## Reading the OGlcNac peaks
    message("Reading the OGlcNac peaks")
    querygr <- unique(buildgr(queryfile))

    ## Perform overlap between O-GlcNac peaks and promoters
    message("Perform overlap between O-GlcNac peaks and promoters")
    resultoverlap <- GenomicRanges::findOverlaps(querygr, promotersgr,
        ignore.strand = FALSE)
    idxkeep <- which(!duplicated(S4Vectors::queryHits(resultoverlap)))
    resultoverlap <- resultoverlap[idxkeep, ]
    promotersgr <- promotersgr[S4Vectors::subjectHits(resultoverlap), ]

    ## Retrieving information on genes
    ensembl <- tryusemart(biomart = "ENSEMBL_MART_ENSEMBL",
        "mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org",
        alternativemirror = usealtmirror)
    symbolstab <- retrievegeneinfo(ensembl, promotersgr)
    idxtable <- match(names(promotersgr),
        symbolstab$ensembl_transcript_id_version)
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

    return(promgff)
}



##################
# MAIN
##################

# Retreives the parameters
getParams(paramsDefinition)

if (!file.exists(outputfolder))
    dir.create(outputfolder, recursive = TRUE)

###
# Reading inputs
###

## Reading gene lengths
genelengthvec <- read.delim(countslength, header = TRUE)

## Reading counts and computing TPM
counts_list <- lapply(countspath, function(x, genelengthvec) {
    fi <- read.delim(x, header = TRUE)
    if (!isTRUE(all.equal(nrow(fi), nrow(genelengthvec))))
        stop("Problem with gene lengths")

    ## Calculate TPM
    countsvec <- fi[, 2] / genelengthvec$Length
    countsvectpm <- (countsvec*1e6)/sum(countsvec)
    fi <- cbind(fi, countsvectpm)
    colnames(fi) <- c("Name", "NumReads", "TPM")
    return(fi)
}, genelengthvec)

## Reading prom glc
promglc <- retrieveglcproms(upstreambp, downstreambp, queryfile,
    usealtmirror)

## Making counts list with only expressed one
counts_nozero_list <- lapply(counts_list, function(x) {
    idxzero <- which(x$NumReads == 0)
    nonZeroIndex(idxzero, "step 1")
    return(x[-idxzero, ])
})


###
# Defining each replicate expression table, for the categories non glc prom,
# random prom, and glc prom and generating the boxplot of TPM values
###

invisible(mapply(function(allcounts, allcountsnozero, currentname, promglc,
outfold) {

    message("Processing ", currentname)

    ## Defining expression of non glc prom
    counts_noglcprom <- defineNoGlcProm(promglc, allcounts)

    ## Defining expression of random prom
    counts_random <- allcountsnozero[sample(seq_len(nrow(allcountsnozero)),
    nrow(promglc)), ]

    ## Retrieving expression of glc prom
    counts_glcprom <- defineGlcProm(promglc, allcounts)

    ## Creating data frame for ggplot
    res <- prepareMatrix(counts_glcprom, counts_noglcprom, counts_random, 
        allcountsnozero)
    df_matrix <- res[[1]]
    category <- res[[2]]
    labelvec <- res[[3]]

    ## Generating the violin plot
    performPlot(labelvec, df_matrix, currentname, outfold)

    ## Computing mann-whitney-wilcoxson on each violin
    wilcoxonTest(counts_glcprom, counts_noglcprom, counts_random,
        allcountsnozero, labelvec, currentname, outfold)
}, counts_list, counts_nozero_list, namescountsvec,
MoreArgs = list(promglc, outputfolder)))
