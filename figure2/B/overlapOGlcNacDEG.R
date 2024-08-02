###############
# This script aimns at representing differing and common intervals for a maximum of 4 gff files.
# Descostes November 2015
# update: feb 2018
# Format: May 2024 (except second part)
# R-4.2.0
###############


library("GenomicFeatures")
library("NGSprofiling")
library("ChIPpeakAnno")
library("rtracklayer")
library("biomaRt")
library("annotate")
library("gplots")
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
col_vec <- c("#E69F00", "#56B4E9", "#E95680")
outformat <- "png"
max_gap <- 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gff_file_vec <- c("results/Glc1.gff", "results/log0_siogtdown-ensembl.gff",
"results/log0_siogtup-ensembl.gff")
outfolder <- "results"
comparisonname <- "glcrep1-2_vs_degsiogt"
expnamevec <- c("glcrep1-2", "Down", "Up")
col_vec <- c("#E69F00", "#56B4E9", "#E95680")
outformat <- "png"
max_gap <- 0

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

################



##############
# MAIN
##############

checkparams(gff_file_vec, expnamevec, outformat)
checkingOutputFolder(outfolder)

message("Reading GFF files and converting to GRanges")
gfflist <- readgff(gff_file_vec)
gffgrlist <- convert2gr(gfflist)

message("Performing the overlap")
peaks1 <- gffgrlist[[1]]
peaks2 <- gffgrlist[[2]]
peaks3 <- gffgrlist[[3]]
ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, maxgap = max_gap,
        connectedPeaks = "keepAll")

# Making the venn diagram venneuler

eulerthree <- function(ol) {
    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1")
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2")
    area3 <- retrieve_peaks_number(ol$peaklist, "peaks3")
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2")
    area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3")
    area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3")
    area1_area2_area3 <- retrieve_peaks_number(ol$peaklist,
        "peaks1///peaks2///peaks3")
}


    if (outformat == "png") {
        png(filename=paste(outfolder, comparisonname, "overlap-chippeakanno-formatted.png",sep=""), width = 1000, height = 1000, bg = "transparent")
    }else if(outformat == "ps"){
        cairo_ps(filename=paste(outfolder, comparisonname, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste0(outfolder, comparisonname,"overlap-chippeakanno-formatted.pdf"), width=10, height=10)
    }
    draw.triple.venn(area1 = area1 + area1_area2 + area1_area3 + area1_area2_area3, 
            area2 = area2 + area1_area2 + area2_area3 + area1_area2_area3, 
            area3 = area3 + area1_area3 + area2_area3 + area1_area2_area3, 
            n12 = area1_area2 + area1_area2_area3, 
            n23 = area2_area3 + area1_area2_area3, 
            n13 = area1_area3 + area1_area2_area3, 
            n123 = area1_area2_area3,
            category = expnamevec,
            euler.d = TRUE, 
            scaled = TRUE, 
            ext.text = FALSE,
            col = col_vec, 
            fill = col_vec,
            cex = rep(2, 7), 
            cat.cex = rep(2.5, 3));
    dev.off();


cat("Retrieving list of element per overlap\n");

outfolder_peaks <- paste(outfolder, "peaks_per_circle/", sep="");
checkingOutputFolder(outfolder_peaks);

for(i in 1:length(ol$peaklist)) 
{
    cat("\t", i, "/", length(ol$peaklist), "\n");
    
    if(nchar(names(ol$peaklist)[i]) > 6)
    {
        index_exp_vec <- as.numeric(unlist(lapply(strsplit(unlist(strsplit(names(ol$peaklist)[i],"///")),"peaks"),"[",2)));
        output_file <- paste(expnamevec[index_exp_vec], collapse="_");
    }
    else
    {
        index_exp <- as.numeric(unlist(lapply(strsplit(names(ol$peaklist)[i],"peaks"),"[",2)));
        output_file <- expnamevec[index_exp];
    }
    
    gff_table <- data.frame(seqname=as.character(seqnames(ol$peaklist[[i]])), 
            source="vennDiagram_overlapGFF", 
            feature = as.character(unlist(lapply(elementMetadata(ol$peaklist[[i]])$peakNames, paste,collapse="-"))), 
            start=start(ranges(ol$peaklist[[i]])),
            end=end(ranges(ol$peaklist[[i]])),
            score=0,
            strand=as.character(strand(ol$peaklist[[i]])),
            frame=".",
            group=".")
    
    write.table(gff_table, file=paste(outfolder_peaks, output_file, ".gff", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE);
}



