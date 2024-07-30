##############
# This script computes the union of intervals from several bed files and output a gff file.
# Descostes, September 2015, update nov 2018
##############


library("IRanges")
library("GenomicRanges")


################
# PARAMETERS
################


filepathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc2_lane1sample13_peaks_narrowPeak.gff")
outputfile <- "/g/romebioinfo/tmp"
inputformat <- "gff"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

filepathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc2_lane1sample13_peaks_narrowPeak.gff")
outputfile <- "/g/romebioinfo/tmp"
inputformat <- "gff"


##############
# MAIN
##############


# Retreives the parameters
getParams(paramsDefinition);

if(inputformat != "bed" && inputformat != "gff")
{
	stop("input format should be bed or gff\n");
}


if(!file.exists(output_path))
{
	dir.create(output_path, recursive = TRUE)
}


# Reading bed files

cat("Reading files\n");

files_vec <- list();

for(i in 1:length(filepathvec)) 
{
	files_vec[[i]] <- read.table(filepathvec[i], stringsAsFactors=FALSE);
}

union_files <- do.call(rbind,files_vec);

if(inputformat == "bed"){
	
	##intervals_rangedData <- RangedData(IRanges(start = union_files[,2], end = union_files[,3]), space = union_files[,1]);
	intervals_gRanges <- GRanges(seqnames= union_files[,1], ranges= IRanges(start = union_files[,2], end = union_files[,3]))
	
}else{
	##intervals_rangedData <- RangedData(IRanges(start = union_files[,4], end = union_files[,5]), space = union_files[,1]);
	intervals_gRanges <- GRanges(seqnames= union_files[,1], ranges= IRanges(start = union_files[,4], end = union_files[,5]))
}

cat("Reducing intervals\n");
intervals_gRanges <- reduce(intervals_gRanges);

gff_file_to_write <- cbind(seqName = as.character(seqnames(intervals_gRanges)), source = "union_files", feature = "union_files", start = start(intervals_gRanges), end = end(intervals_gRanges), 
		score = 0, strand = '+', frame = ".", group = make.unique(rep("group", length(start(intervals_gRanges))), sep="-"));

write.table(gff_file_to_write, file=paste(output_path, outputfile, sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE);
