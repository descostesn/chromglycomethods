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


clusterpathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_1-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_2-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_3-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_4-compartmentsgff/promoters.gff", # nolint
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_5-compartmentsgff/promoters.gff") # nolint
ensemblpath <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/human/hg38/ensembl/Homo_sapiens_GRCh38110_chrfiltered_genes.gff" # nolint
output_folder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_1-compartmentsgff/promoters_vs_ensemblgenes/test",
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_2-compartmentsgff/promoters_vs_ensemblgenes/test",
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_3-compartmentsgff/promoters_vs_ensemblgenes/test",
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_4-compartmentsgff/promoters_vs_ensemblgenes/test",
        "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/cluster_5-compartmentsgff/promoters_vs_ensemblgenes/test")
expname <- "promspeaks_ensemblgenes"
extractensemblgenename <- "TRUE"

!!!!!!!!!!!!!!!!!!! NEED TO ADAPT THE SCRIPT



#############
## FUNCTIONS
#############


checkingoutputfolder <- function(output_path) {
    if (!file.exists(output_path))
        dir.create(output_path, recursive = TRUE)
}

checkparams <- function(gff_file_vec, expname_vec, output_folder) {

    if (length(gff_file_vec) != 2 || length(expname_vec) != 2)
        stop("\n This script takes only two exp as input\n")

    extsuffix <- unique(sapply(basename(gff_file_vec),
        function(x) strsplit(x, "\\.")[[1]][2]))

    if (!isTRUE(all.equal(length(extsuffix), 1)) &&
        !isTRUE(all.equal(extsuffix, "gff")))
            stop("The files should be in gff format.")

    checkingoutputfolder(output_folder)
}


buildgffrangeslist <- function(gff_file_vec, extractensemblgenename) {
    gff_granges_list <- list()

    for (i in seq_len(length(gff_file_vec))) {

        current_gff <- read.delim(gff_file_vec[i], stringsAsFactors = FALSE,
            header = FALSE)

        ## In case the overlap is done with ensembl genes, extract gene name
        if (isTRUE(all.equal(i, 2)) && extractensemblgenename)
            genenamevec <- unlist(lapply(strsplit(unlist(lapply(
                strsplit(current_gff$V9, ";"), "[", 2)), "="), "[", 2))
        else
            genenamevec <- current_gff[, 3]

        ## Check for duplicated names
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

        gff_granges_list[[i]] <- GenomicRanges::GRanges(
            seqnames = current_gff[, 1],
            ranges = IRanges::IRanges(start = current_gff[, 4],
                                  end = current_gff[, 5],
                                  names = genenamevec),
            strand = current_gff[, 7])
    }
    return(gff_granges_list)
}



##############
# MAIN
##############

# Retreives the parameters
checkparams(gff_file_vec, expname_vec, output_folder)

message("Reading gff input and converting to rangedData")
gff_granges_list <- buildgffrangeslist(gff_file_vec, extractensemblgenename)

message("Performing the overlap")
result_overlap <- ChIPpeakAnno::findOverlapsOfPeaks(gff_granges_list[[1]],
                        gff_granges_list[[2]])

message("Writting the overlapping peaks")
gff_table_peak1 <- data.frame(
    seqname = as.character(result_overlap$overlappingPeaks[[1]][, 2]),
    source = "vennDiagram_overlapGFF",
    feature = result_overlap$overlappingPeaks[[1]][, 1],
    start = result_overlap$overlappingPeaks[[1]][, 3],
    end = result_overlap$overlappingPeaks[[1]][, 4],
    score = result_overlap$overlappingPeaks[[1]][, 5],
    strand = as.character(result_overlap$overlappingPeaks[[1]][, 6]),
    frame = ".", group = ".")

gff_table_peak2 <- data.frame(
    seqname = as.character(result_overlap$overlappingPeaks[[1]][, 8]),
    source = "vennDiagram_overlapGFF",
    feature = result_overlap$overlappingPeaks[[1]][, 7],
    start = result_overlap$overlappingPeaks[[1]][, 9],
    end = result_overlap$overlappingPeaks[[1]][, 10],
    score = result_overlap$overlappingPeaks[[1]][, 11],
    strand = result_overlap$overlappingPeaks[[1]][, 12],
    frame = ".", group = ".")

message("Writing output files")
write.table(gff_table_peak1,
    file = file.path(output_folder, paste0(expname_vec[1], ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gff_table_peak2,
    file = file.path(output_folder, paste0(expname_vec[2], ".gff")),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Done.")
