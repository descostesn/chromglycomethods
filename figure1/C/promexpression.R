#################
# This script measure the expression of genes for several categories:
# 1) Genes having a glc peak at promoter
# 2) Genes not having a glc peak at promoter (same nb of prom than previous)
# 3) Randomly selected genes
# 4) All genes
#
# The comparison is done in the form of a violin plot.
#
# Descostes - April 2024 with R-4.3.3
#################

library("ggplot2")
library("RColorBrewer")



##################
# PARAMETERS
##################


countspath
namescountsvec
countslength
promsglcpath
outputfolder


countspath <- "data/ESCRNAseq_SRR11294181counts.txt"
namescountsvec <- "rep1"
countslength <- "data/ESCRNAseq_SRR11294181countslength.txt" # nolint
promsglcpath <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_mouseGlucose_enhancers/HGGlcNAc1s12/renamed/genes_fromProm-HGGlcNAc1s12.gff", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_mouseGlucose_enhancers/HGGlcNAc2s13/renamed/genes_fromProm-HGGlcNAc2s13.gff") #nolint

outputfolder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/expression_vs_glc/violinplot_bycategory_promexpression/hackett_rnaseq/ESCHGGlcNAc1_lane1sample12/", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/expression_vs_glc/violinplot_bycategory_promexpression/hackett_rnaseq/ESCHGGlcNAc2_lane1sample13/") # nolint




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
promglc <- read.delim(promsglcpath, header = FALSE)

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

# allcounts <- counts_list[[1]]
# allcountsnozero <- counts_nozero_list[[1]]
# currentname <- namescountsvec[1]
# outfold <- outputfolder
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
