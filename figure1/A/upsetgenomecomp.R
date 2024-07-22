#################
# This script aims at reproducing fig 1A
# The different categories are defined as follows:
# 1) Active promoter: K27ac peaks in TSS-/+1Kb.
# 2) Transcription initiation: TSS-1/+1kb overlapping with Ser5P peaks.
# 3) Heterochromatin: H3K9me3 peak intervals.
# 4) Bivalent promoters: H3K4me3/H3K27me3 peaks overlapping TSS-/+1Kb. 
# 5) Polycomb domain: Suz12 and RING1B peaks overlapping each other.
# 6) Transcription elongation: TSS+1kb to TES overlapping Ser2P peaks.
# 7) Transcription termination: TES+50bp intervals.
#
# sessionInfo()
# loaded via a namespace (and not attached):
#  TH.data_1.1-0, colorspace_2.0-3, rjson_0.2.21, ellipsis_0.3.2,
# htmlTable_2.4.0, XVector_0.38.0, GenomicRanges_1.50.2, base64enc_0.1-3,
# rstudioapi_0.13,   MatrixModels_0.5-0, bit64_4.0.5, mvtnorm_1.1-3,
# AnnotationDbi_1.60.2, fansi_1.0.3, xml2_1.3.3, org.Dr.eg.db_3.16.0,
# codetools_0.2-18, splines_4.2.0, cachem_1.0.6, chipenrich.data_2.22.0,
# knitr_1.38, Formula_1.2-4, Rsamtools_2.14.0, dbplyr_2.1.1, cluster_2.1.3,
# png_0.1-7, readr_2.1.2, compiler_4.2.0, httr_1.4.2, backports_1.4.1,
# assertthat_0.2.1, Matrix_1.4-1, fastmap_1.1.0, cli_3.6.1, org.Rn.eg.db_3.16.0,
# org.Mm.eg.db_3.16.0, prettyunits_1.1.1, htmltools_0.5.2, quantreg_5.88,
# tools_4.2.0, gtable_0.3.0, glue_1.6.2, GenomeInfoDbData_1.2.9, dplyr_1.0.8, 
# rappdirs_0.3.3, Rcpp_1.0.8.3, Biobase_2.58.0, vctrs_0.6.1, Biostrings_2.66.0,
# nlme_3.1-157, rtracklayer_1.58.0, xfun_0.30, stringr_1.4.0, lifecycle_1.0.3, 
# restfulr_0.0.15, XML_3.99-0.9, polspline_1.1.20, org.Hs.eg.db_3.16.0 
# zoo_1.8-10, zlibbioc_1.44.0, MASS_7.3-57, scales_1.2.0, hms_1.1.1,
# MatrixGenerics_1.10.0, sandwich_3.0-1, parallel_4.2.0,
# SummarizedExperiment_1.28.0 ggupset_0.3.0, SparseM_1.81, RColorBrewer_1.1-3,
# curl_4.3.2, yaml_2.3.5, memoise_2.0.1, gridExtra_2.3, ggplot2_3.4.1,
# rms_6.3-0, biomaRt_2.54.1, rpart_4.1.16, latticeExtra_0.6-29, stringi_1.7.6,
# RSQLite_2.2.12, S4Vectors_0.36.2, BiocIO_1.8.0, checkmate_2.1.0,
# filelock_1.0.2, BiocGenerics_0.44.0, BiocParallel_1.32.6, GenomeInfoDb_1.34.9,
# rlang_1.1.0, pkgconfig_2.0.3, matrixStats_0.62.0, bitops_1.0-7,
# lattice_0.20-45, purrr_0.3.4, GenomicAlignments_1.34.1 patchwork_1.1.1,
# htmlwidgets_1.5.4, bit_4.0.4, tidyselect_1.1.2, plyr_1.8.7, magrittr_2.0.3,
# R6_2.5.1, IRanges_2.32.0, generics_0.1.2, Hmisc_4.7-0, multcomp_1.4-18,
# DelayedArray_0.24.0, DBI_1.1.2, pillar_1.7.0, foreign_0.8-82, mgcv_1.8-40,
# survival_3.3-1, KEGGREST_1.38.0, RCurl_1.98-1.6, nnet_7.3-17, tibble_3.1.6,
# org.Dm.eg.db_3.16.0, crayon_1.5.1, utf8_1.2.2, BiocFileCache_2.6.1,
# tzdb_0.3.0, progress_1.2.2, jpeg_0.1-9, grid_4.2.0, data.table_1.14.2,
# blob_1.2.3, chipenrich_2.22.0, digest_0.6.29, stats4_4.2.0,
# ComplexUpset_1.3.3, munsell_0.5.0
#
# Descostes - R 4.2.0
#################



library("chipgc")



#############
## PARAMS
#############

geneannovec <- c(gencode = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/genesAnnoForEnhancers/gencode.vM25.annotation.gff", # nolint
        refgene = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/genesAnnoForEnhancers/refGeneUCSC-mm10-March2021.gff", # nolint
        refseq = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/genesAnnoForEnhancers/refseqNCBI-mm10-March2021.gff") # nolint

peakspathcategoriesvec <- c(H3K27ac = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_chipatlas_candidates/mouse/HG/0.04/no_model_broad/H3K27ac_SRX19148013_peaks_broadPeak.gff", # nolint
        H3K4me1 = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/public_data/histones_marks/0.04/no_model_broad/Wysocka-H3K4me1-rep1-input1_peaks_broadPeak.gff", # nolint
        H3K27me3 = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/hiddenDomains/public_data/histones_marks/300/Bell-H3K27me3rep1-input1_vis.gff", # nolint
        H3K4me3 = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_chipatlas_candidates/mouse/HG/0.04/no_model_broad/H3K4me3_SRX5382140_peaks_broadPeak.gff", # nolint
        Suz12 = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/public_data/histones_marks/0.04/no_model_broad/Wysocka-Suz12-input2-single_SRR034190_peaks_broadPeak.gff", # nolint
        RING1B = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/public_data/Christensen_Ring1B_SRR10095137/0.001/no_model/RING1B_peaks_narrowPeak.gff", # nolint
        H3K9me3 = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/public_data/histones_marks/0.04/no_model_broad/H3K9me3_SRR925652_peaks_broadPeak.gff", # nolint
        Ser5P = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/public_data/histones_marks/1e-04/no_model_broad/Pombo-RNAPIISer5P-rep2-input1-single_SRR391050_peaks_broadPeak.gff", # nolint
        Ser2P = "/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/public_data/histones_marks/0.03/no_model_broad/Pombo-RNAPIISer2P-input1tech3-single_SRR391039_peaks_broadPeak.gff", # nolint
        ATACSeq = "/g/boulard/Projects/O-N-acetylglucosamine/data/public_data/mouse/mm10/ATAC-Seq_ESC/data/sequencing/macs2/peak-calls/UniqueNoDupeShiftedNFR/narrow/ATACRep2_SRR466767_02.gff") # nolint

peakspathqueryvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc1_lane1sample12_peaks_narrowPeak.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/peak_detection/macs2/Sept2023_glcPolII/mouse/glcGlucose/0.04/no_model/ESCHGGlcNAc2_lane1sample13_peaks_narrowPeak.gff") # nolint

querynamevec <- c("ESCHGGlcNAc1", "ESCHGGlcNAc2")

repeatsannovec <- c(LINE = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/LINE.gff",  # nolint
        LTR = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/LTR.gff",  # nolint
        SINE = "/g/boulard/Projects/O-N-acetylglucosamine/data/Annotations/mouse/mm10/repeatmasker/classes/gff/SINE.gff") # nolint

outputfolder <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/extendedGenomicCompartments/sept2023_mouse_HG/rep1-2" # nolint

glcnacbwvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_GlcPolIIGlucose_Sept2023/glc_glucoce_mouseESC/data/sequencing/bw/genome_norm/ESCHGGlcNAc1_lane1sample12.bw", # nolint
"/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_GlcPolIIGlucose_Sept2023/glc_glucoce_mouseESC/data/sequencing/bw/genome_norm/ESCHGGlcNAc2_lane1sample13.bw") # nolint

includerepeats <- FALSE

outputObjectPath <- "/g/boulard/Projects/O-N-acetylglucosamine/analysis/extendedGenomicCompartments/sept2023_mouse_HG/rep1-2/gclist.Rdat" #nolint

countstable <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_mRNASeq_oct2019/recompute_with_ensembl/data/sequencing/counts/featureCounts/featurecounts_MB1-WT-E14-T2i-replicate1.txt.gz.tabular" #nolint
countslength <- "/g/boulard/Projects/O-N-acetylglucosamine/data/Sofia_mRNASeq_oct2019/recompute_with_ensembl/data/sequencing/counts/featureCounts/MB1-WT-E14-T2i-replicate1.txt.gz.tabular" # nolint
countsannotype <- "ensembl"

#"https://www.ensembl.org"
biomartstr <- "ENSEMBL_MART_ENSEMBL"
datasetstr <- "mmusculus_gene_ensembl"
hoststr <- "https://nov2020.archive.ensembl.org"
alternativemirroropt <- FALSE


#############
## MAIN
#############

## Checking parameters
checkParams(peakspathqueryvec, glcnacbwvec, querynamevec, geneannovec,
        peakspathcategoriesvec, repeatsannovec, countstable, includerepeats,
        countsannotype)

## Define vector names
names(peakspathqueryvec) <- querynamevec

## Connecting to ensembl
geneidtab <- retrieveconversiontab(biomartstr, datasetstr, hoststr,
  alternativemirroropt)

## Start processing for each peak file
# querypath <- peakspathqueryvec[1]; queryname <- querynamevec[1]
# glcnacbw <- glcnacbwvec[1]; peakscat <- peakspathcategoriesvec;
# geneannos <- geneannovec; outputfold <- outputfolder; 
# includerep <- includerepeats; rnacounts <- countstable;
# refseqpath <- geneannovec["refseq"]
gclist <- mapply(function(querypath, queryname, glcnacbw, peakscat, geneannos,
                outputfold, includerep, rnacounts, refseqpath,
                geneidtab, countsannotype) {

        message("Processing ", queryname)

        peakspathvec <- c(querypath, peakscat)
        outfold <- file.path(outputfold, queryname)
        createfolder(outfold)

        ## PART 1: Build initial information
        # Define intervals of the different genomic compartments
        #includerepeats <- includerep
        gc <- buildIntervalsObject(peakspathvec, geneannos, includerep)
        # Retrieve query chip-seq signal on peaks
        # theobject <- gc
        # includerepeats <- includerep
        # bwpath <- glcnacbw
        gc <- retrieveGlcPeakVal(gc, includerep, glcnacbw)

        ## PART 2: Upset diagram of the glucnac peaks overlap with each
        ## compartment
        upsetDiagram(gc, outfold, includerep)

        ## PART 3: Boxplot of glucnac levels accross the different
        ## compartments
        boxplotGlcnacLevels(gc, glcnacbw, outfold, includerep)

        ## PART 4: Extract glucnac coordinates in each compartment
        outputGlcPeaksCoordPerCompartment(gc, outfold,
                peakspathvec[1], includerep)

        ## PART 5: Extract coordinates for compartments having a peak
        gc <- extractCompCoordWithPeak(gc, outfold, includerep)

        ## PART 6: Retrieve gene expression in each compartment. The gene
        ## must be associated with a glc peak
        #theobject <- gc
        #plotviolin <- FALSE
        violinplotExpression(gc, outfold, includerep, rnacounts,
                refseqpath, geneidtab, countsannotype)

        ##!! PART 7: Retrieve expression for ALL genes in each compartment
        ##!! PART 8: Retrieve levels of each mark in each compartment


            return(gc)
        }, peakspathqueryvec, querynamevec, glcnacbwvec,
        MoreArgs = list(peakspathcategoriesvec, geneannovec, outputfolder,
                includerepeats, countstable, geneannovec["refseq"],
                geneidtab, countsannotype), SIMPLIFY = FALSE)

save(gclist, file = outputObjectPath)
#load(outputObjectPath)

## Generating the complex upset
complexUpsetDiagram(gclist, includerepeats, outputfolder)
