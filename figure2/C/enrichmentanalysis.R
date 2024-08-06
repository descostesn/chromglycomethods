##############
# This script performs gene ontologies from gene lists with cluster profiler.
# It enables also the comparison of ontologies between different samples.
# Performed with R 4.3.2
# Descostes Feb 2020 - April 2023 - update Dec 2023
##############



pacman::p_load(clusterProfiler, ReactomePA, qusage, topGO, Rgraphviz, ggplot2,
    stringr, biomaRt)


################
# PARAMETERS
################


gff_vec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_down-ensembl.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_up-ensembl.gff", 
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHGmergedRep1-2_vs_downsiogt/siogtdown.gff")
expnames_vec <- c("down", "up", "downglc")
species_name <- "mouse"
output_folder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/clusterProfiler/rnaseq_siogt_formichetti/down_up_downglc/test")
output_format <- "png"

!!!!!!!!!!!!

gff_vec <- c("data/log0_down-ensembl.gff", "data/log0_up-ensembl.gff",
    "data/siogtdown.gff")
expnames_vec <- c("down", "up", "downglc")
species_name <- "mouse"
output_folder <- c("results")
output_format <- "png"



################




################
# FUNCTION
################


source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/tools/clusterProfiler/utils/biomart.R") #nolint
source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/tools/clusterProfiler/utils/comparisons.R") #nolint
source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/tools/clusterProfiler/utils/go.R") #nolint
source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/tools/clusterProfiler/utils/io.R") #nolint
source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/tools/clusterProfiler/utils/plots.R") #nolint
source("/g/boulard/Projects/O-N-acetylglucosamine/src/R/tools/clusterProfiler/utils/utils.R") #nolint

################



##############
# MAIN
##############


## Building gff and expnames lists
names(gffvec) <- expnamesvec

## Checking and retrieving parameters
checkParams(species_name, gffvec, expnamesvec, output_folder)
database_name <- returnDB(species_name)
kegg_name <- returnDB(species_name, kegg = TRUE)


## Reading gff files
message("Reading gff files")
different_id_list <- lapply(gffvec, function(currentgff, dbname) {
    fi <- read.delim(currentgff, header = FALSE)
    return(clusterProfiler::bitr(fi[, 3], fromType = "ENSEMBL",
                                    toType = c("SYMBOL", "ENTREZID",
                                            "ENSEMBL", "UNIPROT",
                                            "ENSEMBLPROT"), OrgDb = dbname))
}, database_name)

## Defining background
message("Defining background")
background_id_vec <- backgroundDef(different_id_list, backgroundpath,
        species_name)

## Retrieving info from biomart
message("Retrieving info from biomart")

if (isTRUE(all.equal(species_name, "mouse"))) {
        datasetname <- "mmusculus_gene_ensembl"
} else {
        datasetname <- "hsapiens_gene_ensembl"
}

res <- retrieveInfo(biomartname = "ENSEMBL_MART_ENSEMBL",
        datasetname =  datasetname,
        versionname = "Ensembl Genes 109", alternativemirrorchoice = TRUE)
ensembl <- res[[1]]
genesinfo <- res[[2]]

## Filtering out non-canonical chromosomes from genesinfo
message("Filtering out non-canonical chromosomes from genesinfo")
idxremove <- unique(c(grep("PATCH", genesinfo[, 1]),
                      grep("HS", genesinfo[, 1]),
                      grep("GL", genesinfo[, 1]),
                      grep("KI", genesinfo[, 1]),
                      grep("MT", genesinfo[, 1]),
                      grep("HG", genesinfo[, 1])))

if (!isTRUE(all.equal(length(idxremove), 0))) {
    message("\t Removing ", length(idxremove), "/", nrow(genesinfo),
        " annotations with non canonical chromosomes")
    genesinfo <- genesinfo[-idxremove, ]
}

## Computing GO, bar, dot and cnet plots
message("Computing GO, bar, dot and cnet plots")

invisible(mapply(function(diff_ids, currentname, dbname, levelnum, background,
                        kegg, species, outfold, outformat, ens, infos) {

            message("\t Processing ", currentname)
            results <- caclulate_initial_GO(diff_ids, dbname, levelnum,
                    background, kegg, species)
            names(results) <- c("ggoCC", "ggoBP", "ggoMF", "egoCC", "egoBP",
                    "egoMF", "kk", "mkk", "react")
            loadplots(results, outfold, levelnum, currentname, outformat,
                                ens, infos, diff_ids)
        }, different_id_list, names(different_id_list),
        MoreArgs = list(database_name, levelnum, background_id_vec, kegg_name,
                species_name, output_folder, output_format, ensembl,
                genesinfo), SIMPLIFY = FALSE))


if (length(gffvec) > 1) {

    message("Creating the input list of entrezID")
    idcomplist <- lapply(different_id_list,
        function(x) return(unique(x$ENTREZID)))

    message("Creating the entrezID-symbol table")
    idtable <- do.call("rbind", different_id_list)
    idx <- which(duplicated(paste(idtable$SYMBOL, idtable$ENTREZID, sep = "-")))
    idtable <- idtable[-idx, ]

    message("Performing clusters comparison")
    results <- clustersComparison(idcomplist, database_name, levelnum,
        kegg_name, species_name, background_id_vec)

    message("Output dotplots of the comparisons")
    outfoldcomp <- file.path(output_folder, levelnum, "compare_cluster")
    if (!file.exists(outfoldcomp))
        dir.create(outfoldcomp, recursive = TRUE)

    ## Removing empty categories
    results <- results[sapply(results, function(x) !is.null(nrow(x)))]

    if (!isTRUE(all.equal(length(results), 0)))
        dotplotComparisons(results, output_format, outfoldcomp, idtable,
            ensembl, genesinfo, different_id_list)
    else
        message("No enrichment found for any of the categories.")
}



message("Done")
