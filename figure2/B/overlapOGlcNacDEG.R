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
# library("org.Mm.eg.db")


#################
# PARAMETERS
################


gff_file_vec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHGPeaksRep1_vs_Rep2/Glc1.gff", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_down-ensembl.gff", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_up-ensembl.gff")  # nolint
output_folder <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHG_vs_DEGsiogt/mergedrep1-2/test" # nolint
comparison_title <- "mergedrep1-2_vs_DEGsiogt"
expnames_tab <- c("mergedrep1-2", "Down", "Up")
genome_version <- "mm10"
organism_name <- "mouse"
center_criterion <- "max"
col_vec <- c("#E69F00", "#56B4E9", "#E95680")
output_format <- "png"





paramsDefinition[["--bigwigNameVec"]] <- list(variableName="bigwig_name_vec", numeric=F, mandatory=F, description="Vector of big wig file names.", default=NA) # nolint
paramsDefinition[["--maxgap"]] <- list(variableName="max_gap", numeric=T, mandatory=F, description="Non-negative integer. Peak intervals with a separation of maxgap or less are considered to be overlapped.", default=0) # nolint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gff_file_vec <- c("results/Glc1.gff", "results/log0_siogtdown-ensembl.gff",
"results/log0_siogtup-ensembl.gff")
output_folder <- "results"
comparison_title <- "Glcrep1-2_vs_DEGsiogt"
expnames_tab <- "Glcrep1-2 Down Up"
genome_version <- "mm10"
organism_name <- "mouse"
center_criterion <- "max"
col_vec <- c("#E69F00", "#56B4E9", "#E95680")
output_format <- "png"


################


################
# FUNCTION
################

retrieve_peaks_number <- function(peak_list, label_name) {
    peak_number <- lapply(peak_list, length)[[label_name]]
    if (is.null(peak_number)) return(0) else return(peak_number)
}

checkparams2 <- function(output_format) {

    if (!isTRUE(all.equal(output_format, "ps")) &&
        !isTRUE(all.equal(output_format, "png")) &&
        !isTRUE(all.equal(output_format, "pdf")))
        stop("output_format should be png, pdf or ps")
}

checkparams <- function(gff_file_vec, expnames_tab, center_criterion,
    output_format) {

    if (length(gff_file_vec) > 5)
        stop("This script takes at most 5 gff files as input")

    if (length(gff_file_vec) < 2)
        stop("This script takes a minimum of 2 gff files as input")

    if (!isTRUE(all.equal(length(gff_file_vec), length(expnames_tab))))
        stop("The number of exp names should be equal to the number of gff ",
            "files")

    if (!isTRUE(all.equal(organism_name, "mouse")) &&
        !isTRUE(all.equal(organism_name, "human")))
            stop("The organism name should be mouse or human")

    if (!isTRUE(all.equal(center_criterion, "coordinates")) &&
        !isTRUE(all.equal(center_criterion, "max")))
        stop("Center criterion should be 'coordinates' or 'max'")

    checkparams2(output_format)
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

loadorg <- function(organism_name) {
    if (isTRUE(all.equal(organism_name, "mouse"))) {
        library("org.Mm.eg.db")
    }else {
        library("org.Hs.eg.db")
    }
}

################



##############
# MAIN
##############

checkparams(gff_file_vec, expnames_tab, organism_name, center_criterion,
    output_format)
checkingOutputFolder(output_folder)
loadorg(organism_name)

message("Reading GFF files and converting to GRanges")
gfflist <- readgff(gff_file_vec)
gffgrlist <- convert2gr(gfflist)

message("Performing the overlap")

if (length(gffgrlist) == 2) {
    peaks1 <- gffgrlist[[1]]
    peaks2 <- gffgrlist[[2]]
    ol <- findOverlapsOfPeaks(peaks1, peaks2, maxgap = max_gap,
        connectedPeaks = "keepAll")
}else if (length(gffgrlist) == 3) {
    peaks1 <- gffgrlist[[1]]
    peaks2 <- gffgrlist[[2]]
    peaks3 <- gffgrlist[[3]]

    ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, maxgap = max_gap,
        connectedPeaks = "keepAll")
}else if (length(gffgrlist) == 4) {
    peaks1 <- gffgrlist[[1]]
    peaks2 <- gffgrlist[[2]]
    peaks3 <- gffgrlist[[3]]
    peaks4 <- gffgrlist[[4]]

    ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, peaks4, maxgap = max_gap,
        connectedPeaks = "keepAll")
} else if (length(gffgrlist) == 5) {
    peaks1 <- gffgrlist[[1]]
    peaks2 <- gffgrlist[[2]]
    peaks3 <- gffgrlist[[3]]
    peaks4 <- gffgrlist[[4]]
    peaks5 <- gffgrlist[[5]]

    ol <- findOverlapsOfPeaks(peaks1, peaks2, peaks3, peaks4, peaks5,
        maxgap = max_gap, connectedPeaks = "keepAll")
}

message("Making the venn diagram")
fileout <- paste0(comparison_title, "overlap-chippeakanno")

if (isTRUE(all.equal(output_format, "png"))) {
    png(filename = file.path(output_folder, paste0(fileout, ".png")),
      width = 1000, height = 1000, bg = "transparent")
}else if (isTRUE(all.equal(output_format, "ps"))){
    cairo_ps(filename = file.path(output_folder, paste0(fileout, ".ps")),
      width = 7, height = 7, bg = "transparent")
}else if(isTRUE(all.equal(output_format, "pdf"))) {
    pdf(file = file.path(output_folder, paste0(fileout, ".pdf")),
      width = 10, height = 10)
}else {
    stop("Format ", output_format, " is not supported")
}
result <- makeVennDiagram(ol, NameOfPeaks= expnames_tab,
    totalTest = max(sapply(gfflist, nrow)) * 2, ignore.strand = TRUE,
    connectedPeaks = "keepAll")
dev.off()

result$p.value[, "pval"] <- p.adjust(result$p.value[, "pval"], method = "BH")
write.table(result$p.value, file = file.path(output_folder, "pvalue_stats.txt"),
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


# Making the venn diagram venneuler

if (length(gffgrlist) == 2) {

    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1")
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2")
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2")

    if (output_format == "png") {
        png(filename=paste(output_folder, comparison_title,
          "overlap-chippeakanno-formatted.png", sep = ""), width = 1000,
          height = 1000, bg = "transparent")
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, comparison_title,
          "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7,
          bg = "transparent")
    } else {
        pdf(file=paste0(output_folder, comparison_title,
          "overlap-chippeakanno-formatted.pdf"), width = 10, height = 10)
    }
    draw.pairwise.venn(area1 = area1 + area1_area2,
            area2 = area2 + area1_area2,
            cross.area = area1_area2,
            category = expnames_tab,
            euler.d = TRUE,
            scaled = TRUE,
            ext.text = FALSE,
            col = col_vec,
            fill = col_vec,
            cex = rep(2, 3),
            cat.cex = rep(2.5, 2))
    dev.off()

}else if (length(gffgrlist) == 3) {

    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1")
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2")
    area3 <- retrieve_peaks_number(ol$peaklist, "peaks3")
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2")
    area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3")
    area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3")
    area1_area2_area3 <- retrieve_peaks_number(ol$peaklist,
        "peaks1///peaks2///peaks3")

    if (output_format == "png") {
        png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 1000, height = 1000, bg = "transparent")
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste0(output_folder, comparison_title,"overlap-chippeakanno-formatted.pdf"), width=10, height=10)
    }
    draw.triple.venn(area1 = area1 + area1_area2 + area1_area3 + area1_area2_area3, 
            area2 = area2 + area1_area2 + area2_area3 + area1_area2_area3, 
            area3 = area3 + area1_area3 + area2_area3 + area1_area2_area3, 
            n12 = area1_area2 + area1_area2_area3, 
            n23 = area2_area3 + area1_area2_area3, 
            n13 = area1_area3 + area1_area2_area3, 
            n123 = area1_area2_area3,
            category = expnames_tab,
            euler.d = TRUE, 
            scaled = TRUE, 
            ext.text = FALSE,
            col = col_vec, 
            fill = col_vec,
            cex = rep(2, 7), 
            cat.cex = rep(2.5, 3));
    dev.off();
    
}else if(length(gffgrlist) == 4){
    
    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1");  
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2");
    area3 <- retrieve_peaks_number(ol$peaklist, "peaks3");
    area4 <- retrieve_peaks_number(ol$peaklist, "peaks4");
    
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2");
    area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3");
    area1_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks4");
    area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3");
    area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks4");
    area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks4");
    
    area1_area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3");
    area1_area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks4");
    area1_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks4");
    area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks4");
    
    area1_area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks4");
    
    if(output_format == "png"){
        png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 1000, height = 1000, bg = "transparent")
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste0(output_folder, comparison_title,".pdf"), width=10, height=10)
    }
    draw.quad.venn(area1 = area1 + area1_area2 + area1_area3 + area1_area4 + area1_area2_area3 + area1_area2_area4 + area1_area3_area4 + area1_area2_area3_area4, 
            area2 = area2 + area1_area2 + area2_area3 + area2_area4 + area1_area2_area3 + area1_area2_area4 + area2_area3_area4 + area1_area2_area3_area4, 
            area3 = area3 + area1_area3 + area2_area3 + area3_area4 + area1_area2_area3 + area1_area3_area4 + area2_area3_area4 + area1_area2_area3_area4, 
            area4 = area4 + area1_area4 + area2_area4 + area3_area4 + area1_area2_area4 + area1_area3_area4 + area2_area3_area4 + area1_area2_area3_area4, 
            n12 = area1_area2 + area1_area2_area3 + area1_area2_area4 + area1_area2_area3_area4, 
            n13 = area1_area3 + area1_area3_area4 + area1_area2_area3 + area1_area2_area3_area4, 
            n14 = area1_area4 + area1_area2_area4 + area1_area3_area4 + area1_area2_area3_area4 , 
            n23 = area2_area3 + area1_area2_area3 + area2_area3_area4 + area1_area2_area3_area4, 
            n24 = area2_area4 + area1_area2_area4 + area2_area3_area4 + area1_area2_area3_area4,
            n34 = area3_area4 + area1_area3_area4 + area2_area3_area4 + area1_area2_area3_area4, 
            n123 = area1_area2_area3 + area1_area2_area3_area4, 
            n124 = area1_area2_area4 + area1_area2_area3_area4, 
            n134 = area1_area3_area4 + area1_area2_area3_area4, 
            n234 = area2_area3_area4 + area1_area2_area3_area4, 
            n1234 = area1_area2_area3_area4, 
            category = expnames_tab, 
            euler.d = TRUE, 
            scaled = TRUE, 
            ext.text = FALSE,
            col = col_vec, 
            fill = col_vec,
            cex = rep(2, 15), 
            cat.cex = rep(2.5, 4));
    dev.off();
}else{
    
    
    
    area1 <- retrieve_peaks_number(ol$peaklist, "peaks1");  
    area2 <- retrieve_peaks_number(ol$peaklist, "peaks2");
    area3 <- retrieve_peaks_number(ol$peaklist, "peaks3");
    area4 <- retrieve_peaks_number(ol$peaklist, "peaks4");
    area5 <- retrieve_peaks_number(ol$peaklist, "peaks5");
    
    area1_area2 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2");
    area1_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3");
    area1_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks4");
    area1_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks5");
    area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3");
    area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks4");
    area2_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks5");
    area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks4");
    area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks5");
    area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks4///peaks5");
    
    area1_area2_area3 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3");
    area1_area2_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks4");
    area1_area2_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks5");
    area1_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks4");
    area1_area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks5");
    area1_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks4///peaks5");
    area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks4");
    area2_area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks5");
    area2_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks4///peaks5");
    area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks3///peaks4///peaks5");
    
    area1_area2_area3_area4 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks4");
    area1_area2_area3_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks5");
    area1_area2_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks4///peaks5");
    area1_area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks3///peaks4///peaks5");
    area2_area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks2///peaks3///peaks4///peaks5");
    
    area1_area2_area3_area4_area5 <- retrieve_peaks_number(ol$peaklist, "peaks1///peaks2///peaks3///peaks4///peaks5");
    
    if(output_format == "png"){
        png(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.png",sep=""), width = 1000, height = 1000, bg = "transparent")
    }else if(output_format == "ps"){
        cairo_ps(filename=paste(output_folder, comparison_title, "overlap-chippeakanno-formatted.ps",sep=""), width = 7, height = 7, bg = "transparent");
    }else{
        pdf(file=paste0(output_folder, comparison_title,".pdf"), width=10, height=10)
    }
    
    draw.quintuple.venn(
            area1 = area1 + area1_area2 + area1_area3 + area1_area4 + area1_area5 + area1_area2_area3 + area1_area2_area4 + area1_area2_area5 + area1_area3_area4 + area1_area3_area5 + area1_area4_area5 + area1_area2_area3_area4 + 
                    area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5, 
            area2 = area2 + area1_area2 + area2_area3 + area2_area4 + area2_area5 + area1_area2_area3 + area1_area2_area4 + area1_area2_area5 + area2_area3_area4 + area2_area3_area5 + area2_area4_area5 + area1_area2_area3_area4 +
                    area1_area2_area3_area5 + area1_area2_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5, 
            area3 = area3 + area1_area3 + area2_area3 + area3_area4 + area3_area5 + area1_area2_area3 + area1_area3_area4 + area1_area3_area5 + area2_area3_area4 + area2_area3_area5 + area3_area4_area5 + area1_area2_area3_area4 +
                    area1_area2_area3_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5, 
            area4 = area4 + area1_area4 + area2_area4 + area3_area4 + area4_area5 + area1_area2_area4 + area1_area3_area4 + area1_area4_area5 + area2_area3_area4 + area2_area4_area5 + area3_area4_area5 + area1_area2_area3_area4 +
                    area1_area2_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            area5 = area5 + area1_area5 + area2_area5 + area3_area5 + area4_area5 + area1_area2_area5 + area1_area3_area5 + area1_area4_area5 + area2_area3_area5 + area2_area4_area5 + area3_area4_area5 + area1_area2_area3_area5 +
                    area1_area2_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n12 = area1_area2 + area1_area2_area3 + area1_area2_area4 + area1_area2_area5 + area1_area2_area3_area4 + area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area2_area3_area4_area5, 
            n13 = area1_area3 + area1_area2_area3 + area1_area3_area4 + area1_area3_area5 + area1_area2_area3_area4 + area1_area2_area3_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5, 
            n14 = area1_area4 + area1_area2_area4 + area1_area3_area4 + area1_area4_area5 + area1_area2_area3_area4 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n15 = area1_area5 + area1_area2_area5 + area1_area3_area5 + area1_area4_area5 + area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n23 = area2_area3 + area1_area2_area3 + area2_area3_area4 + area2_area3_area5 + area1_area2_area3_area4 + area1_area2_area3_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5, 
            n24 = area2_area4 + area1_area2_area4 + area2_area3_area4 + area2_area4_area5 + area1_area2_area3_area4 + area1_area2_area4_area5  + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n25 = area2_area5 + area1_area2_area5 + area2_area3_area5 + area2_area4_area5 + area1_area2_area3_area5 + area1_area2_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n34 = area3_area4 + area1_area3_area4 + area2_area3_area4 + area3_area4_area5 + area1_area2_area3_area4 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n35 = area3_area5 + area1_area3_area5 + area2_area3_area5 + area3_area4_area5 + area1_area2_area3_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n45 = area4_area5 + area1_area4_area5 + area2_area4_area5 + area3_area4_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n123 = area1_area2_area3 + area1_area2_area3_area4 + area1_area2_area3_area5 + area1_area2_area3_area4_area5, 
            n124 = area1_area2_area4 + area1_area2_area3_area4 + area1_area2_area4_area5 + area1_area2_area3_area4_area5, 
            n125 = area1_area2_area5 + area1_area2_area3_area5 + area1_area2_area4_area5 + area1_area2_area3_area4_area5,
            n134 = area1_area3_area4 + area1_area2_area3_area4 + area1_area3_area4_area5 + area1_area2_area3_area4_area5, 
            n135 = area1_area3_area5 + area1_area2_area3_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n145 =  area1_area4_area5 + area1_area2_area4_area5 + area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n234 = area2_area3_area4 + area1_area2_area3_area4 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,,
            n235 =  area2_area3_area5 + area1_area2_area3_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n245 = area2_area4_area5 + area2_area3_area4_area5 + area1_area2_area4_area5 + area1_area2_area3_area4_area5,
            n345 = area3_area4_area5 + area1_area3_area4_area5 + area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n1234 = area1_area2_area3_area4 + area1_area2_area3_area4_area5,
            n1235 = area1_area2_area3_area5 + area1_area2_area3_area4_area5,
            n1245 =  area1_area2_area4_area5 + area1_area2_area3_area4_area5,
            n1345 = area1_area3_area4_area5 + area1_area2_area3_area4_area5,
            n2345 =  area2_area3_area4_area5 + area1_area2_area3_area4_area5,
            n12345 =area1_area2_area3_area4_area5,
            category = expnames_tab, 
            euler.d = TRUE, 
            scaled = TRUE, 
            ext.text = FALSE,
            col = col_vec, 
            fill = col_vec,
            cex = rep(2, 31), 
            cat.cex = rep(2.5, 5));
    dev.off();
    
}


cat("Retrieving list of element per overlap\n");

output_folder_peaks <- paste(output_folder, "peaks_per_circle/", sep="");
checkingOutputFolder(output_folder_peaks);

for(i in 1:length(ol$peaklist)) 
{
    cat("\t", i, "/", length(ol$peaklist), "\n");
    
    if(nchar(names(ol$peaklist)[i]) > 6)
    {
        index_exp_vec <- as.numeric(unlist(lapply(strsplit(unlist(strsplit(names(ol$peaklist)[i],"///")),"peaks"),"[",2)));
        output_file <- paste(expnames_tab[index_exp_vec], collapse="_");
    }
    else
    {
        index_exp <- as.numeric(unlist(lapply(strsplit(names(ol$peaklist)[i],"peaks"),"[",2)));
        output_file <- expnames_tab[index_exp];
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
    
    write.table(gff_table, file=paste(output_folder_peaks, output_file, ".gff", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE);
}



