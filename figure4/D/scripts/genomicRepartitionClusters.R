###############################################################################
# This script aims at determining the type of peaks in the different groups
# of the heatmap produced by using the global genomic repartition of the
# union of peaks.
#
# Descostes - R-4.3.2
###############################################################################

library(dplyr)
library(RColorBrewer)


##################
# PARAMETERS
##################


coordpath <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/heatmapsandprofiles/sept2023Glc/human/polIIGlc/unionpeaks/coordheatmapunion5groups.bed" # nolint

peakspathvec <- c(
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-DNA.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-enhancers.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-firstExons.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-introns.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-LINE.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-Low_complexity.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-LTR.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-otherExons.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-otherLocations.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-promoters.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-Satellite.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-Simple_repeat.gff",
    "/g/boulard/Projects/O-N-acetylglucosamine/analysis/genomicDistribution/Sept2023_glcPolII_enhancers_unionofpeaks/unionPeaksPolIIGlc/peaks-SINE.gff"
) # nolint

unionfilepath <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/makeunion/sept2023/human/glcPolII_samples1-2-3-4/test/union_OGlcNac_noauxaux-fig4C.bed"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

coordpath <- "data/peakscoord-fig4C.bed"

peakspathvec <- c("results/peaks-DNA.gff", "results/peaks-enhancers.gff",
    "results/peaks-firstExons.gff", "results/peaks-introns.gff",
    "results/peaks-LINE.gff", "results/peaks-Low_complexity.gff",
    "results/peaks-LTR.gff", "results/peaks-otherExons.gff",
    "results/peaks-otherLocations.gff", "results/peaks-promoters.gff",
    "results/peaks-Satellite.gff", "results/peaks-Simple_repeat.gff",
    "results/peaks-SINE.gff")

unionfilepath <- "data/union_OGlcNac_noauxaux-fig4C.bed"


##################
# FUNCTIONS
##################

computecompfreq <- function(clustcoord, unioncoord, compcoord, outfold) {
    clustername <- unique(clustcoord$cluster)
    message("Processing ", clustername)

    ## Retrieving the rows in the union file that was used to build the
    ## heatmap. The row number is indicated in the string of the source
    ## column of clustcoord.
    idxunion <- as.numeric(gsub("union_files_r", "", clustcoord$source))
    idxunionna <- which(is.na(idxunion))
    if (!isTRUE(all.equal(length(idxunionna), 0))) {
        message(length(idxunionna), " peaks do not have a ref in the union")
        idxunion <- idxunion[-idxunionna]
    }
    unionclust <- unioncoord[idxunion, ]

    ## Makes the row of the union file corresponding with the compartment file
    res <- unionclust %>% left_join(compcoord, by = c("chrom", "start", "end"))

    if (!isTRUE(all.equal(nrow(res), nrow(clustcoord) - length(idxunionna)))) {
        stop("Missing rows in the results")
    }

    idxna <- which(is.na(res$compartment.y))
    nbna <- length(idxna)

    if (!isTRUE(all.equal(nbna, 0))) {
        message("There are ", nbna, "/", nrow(res), " peaks without comp")
        res <- res[-idxna, ]
    }

    ## Splitting coord by compartment
    resbycomplist <- split(res, as.factor(res$compartment.y))

    ## Formatting each element to gff
    resbycompgfflist <- lapply(resbycomplist, function(currentref) {
        return(data.frame(seqname = currentref$chrom,
        source = "detailsGroupsHeatmapHumanSept2023",
        feature = currentref$compartment.y,
        start = currentref$start,
        end = currentref$end,
        score = currentref$score.y, strand = currentref$strand.y,
        frame = currentref$frame.x, group = currentref$group.y))
    })

    ## Write the peak coord associated to its compartment in the heatmap folder
    outfoldtmp <- file.path(outfold, paste0(clustername, "-compartmentsgff"))
    if (!file.exists(outfoldtmp))
        dir.create(outfoldtmp, recursive = TRUE)
    invisible(mapply(function(clustergff, namecomp, clustername, outfoldtmp) {
        write.table(clustergff, file = file.path(outfoldtmp,
            paste0(namecomp, ".gff")), sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE)
    }, resbycompgfflist, names(resbycompgfflist),
        MoreArgs = list(clustername, outfoldtmp)))

    ## Computing the frequency of each compartment
    frqcomp <- table(res$compartment.y)
    return(frqcomp)
}

##################
# MAIN
##################


## Reading the peaks coordinates of the heatmap and spliting by clusters
coord <- read.table(coordpath, header = FALSE, stringsAsFactors = FALSE)
colnames(coord) <- c(
    "chrom", "start", "end", "source", "score", "strand",
    "start2", "end2", "score2", "strand2", "length", "id", "cluster"
)
coordgrouplist <- split(coord, as.factor(coord$cluster))
output_folder <- dirname(coordpath)

## Reading peaks coord for each compartment of the genomic repartition
complist <- lapply(peakspathvec, read.table,
    header = FALSE,
    stringsAsFactors = FALSE
)
comptab <- do.call("rbind", complist)
colnames(comptab) <- c(
    "chrom", "source", "compartment", "start", "end", "score",
    "strand", "group", "frame"
)

## Reading the union of peaks file
uniontab <- read.table(unionfilepath, header = FALSE, stringsAsFactors = FALSE)
colnames(uniontab) <- c("chrom", "start", "end", "source", "score", "strand",
    "start2", "end2", "itemrgb", "blockcount", "blocksize")

## For each cluster group in coordgrouplist, retrieve the genomic compartment
frqcomplist <- lapply(coordgrouplist, function(
    clustcoord, compcoord,
    unioncoord, outfold) {
    return(computecompfreq(clustcoord, unioncoord, compcoord, outfold))
}, comptab, uniontab, output_folder)


##  Plotting piecharts of compartments
piecolorvec <- c(brewer.pal(n = 12, name = "Paired"), "aliceblue")
nbelementvec <- lengths(frqcomplist)
idmax <- which.max(nbelementvec)
names(piecolorvec) <- names(frqcomplist[[idmax]])

invisible(mapply(function(frqcomp, clustname, outfold, colvec) {
    pdf(file = file.path(outfold, paste0(clustname, "-piechart.pdf")),
            width = 10, height = 10)
    pie(frqcomp, labels = paste(names(frqcomp), frqcomp, sep = "-"),
        col = piecolorvec[names(frqcomp)], main = clustname)
    dev.off()
}, frqcomplist, names(frqcomplist), MoreArgs = list(output_folder,
    piecolorvec)))
