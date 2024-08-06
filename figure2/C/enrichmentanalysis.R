##############
# This script performs gene ontologies from gene lists with cluster profiler.
# It enables also the comparison of ontologies between different samples.
# Performed with R 4.3.2
# Descostes Feb 2020 - April 2023 - update Dec 2023
##############


library("org.Mm.eg.db")
library("clusterProfiler")

# pacman::p_load(clusterProfiler, ReactomePA, qusage, topGO, Rgraphviz, ggplot2,
#     stringr, biomaRt)


################
# PARAMETERS
################


gffvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_down-ensembl.gff",
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/differential_analysis/deseq2/RNASeq_siogt_formichetti/lfc0/upanddown/log0_up-ensembl.gff", 
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/venndiagrams/Sofia_polIIGlc_sept2023/mouse/glucose-levels/glcHGmergedRep1-2_vs_downsiogt/siogtdown.gff")
expnamesvec <- c("down", "up", "downglc")
species_name <- "mouse"
output_folder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/clusterProfiler/rnaseq_siogt_formichetti/down_up_downglc/test")
output_format <- "png"
database_name <- "org.Mm.eg.db"
kegg_name <- "mmu"
backgroundpath <- NA
datasetname <- "mmusculus_gene_ensembl"
ensemblversion <- "Ensembl Genes 109"
levelnum <- 3
!!!!!!!!!!!!

gffvec <- c("data/log0_down-ensembl.gff", "data/log0_up-ensembl.gff",
    "data/siogtdown_withOGlcNac.gff")
expnamesvec <- c("down", "up", "downglc")
species_name <- "mouse"
output_folder <- c("results")
output_format <- "png"
database_name <- "org.Mm.eg.db"
kegg_name <- "mmu"
backgroundpath <- NA
datasetname <- "mmusculus_gene_ensembl"
ensemblversion <- "Ensembl Genes 109"
levelnum <- 3

################




################
# FUNCTION
################

checkparams <- function(species_name, gff_list, expnames_list, output_folder) {

    if (!isTRUE(all.equal(species_name, "human")) &&
            !isTRUE(all.equal(species_name, "mouse")))
        stop("\n The only supported species are mouse and human\n")

    if (!isTRUE(all.equal(length(gff_list), length(expnames_list))))
        stop("\n One name has to be given per gff file\n")

    invisible(mapply(function(currentgff, currentnames) {

        if (!isTRUE(all.equal(length(currentgff), length(currentnames))))
            stop("The nb of names is different than nb of files in the lists.")
    }, gff_list, expnames_list))

    if (!isTRUE(all.equal(length(output_folder), 1)))
        stop("\n outputfolder should be unique\n")

    if (!file.exists(output_folder))
        dir.create(output_folder, recursive = TRUE)
}

backgrounddef <- function(activeidslist, backgroundpath, species_name) { #nolint

    if (is.na(backgroundpath)) {
        message("\n\n ### No background provided\n\n")
        background_id_vec <- unique(unlist(lapply(activeidslist,
                                function(x) return(unique(x$ENTREZID)))))
    }else {
        message("\n\n ### Reading background file\n\n")
        background_id_vec <- readLines(backgroundpath)
        ## Convert the ensembl IDs to entrez IDs
        if (isTRUE(all.equal(species_name, "human")))
            background_id_vec <- as.character(
                    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                    background_id_vec, "ENTREZID", "ENSEMBL"))
        else if (isTRUE(all.equal(species_name, "mouse")))
            background_id_vec <- as.character(
                    AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                        background_id_vec, "ENTREZID", "ENSEMBL"))
        else
            stop("Problem with the species in backgrounddef")
    }
    return(background_id_vec)
}

tryGetBM <- function(attributes, ensembl, values = NULL, filters = NULL) { #nolint

    c <- 1
    repeat {
        message("# Attempt ", c, "/5 # ", "Retrieving information about genes",
            " from biomaRt ...")

        if (is.null(values) && is.null(filters))
            res <- try(biomaRt::getBM(attributes = attributes, mart = ensembl),
                silent = TRUE)
        else
            res <- try(biomaRt::getBM(attributes = attributes, mart = ensembl,
                values = values, filters = filters), silent = TRUE)

        if (isTRUE(is(res, "try-error"))) {
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


## Code taken from sharedInternals.R of the CONCLUS package
tryUseMart <- function(biomart = "ensembl", dataset, version, #nolint
                       alternativeMirror) { #nolint
    c <- 1

    repeat {
        message("# Attempt ", c, "/5 # ", "Connection to Ensembl ... ")

        if (!alternativeMirror)
            ensembl <- try(biomaRt::useMart(biomart, dataset = dataset,
                version = version), silent = TRUE)
        else
            ensembl <- try(biomaRt::useEnsembl(biomart, dataset = dataset,
                mirror = "useast"), silent = TRUE)

        if (isTRUE(is(ensembl, "try-error"))) {
            c <- c + 1
            error_type <- attr(ensembl, "condition")
            message(error_type$message)

            if (c > 5) {
                stop(
                    "There is a problem of connexion to Ensembl for ",
                    "now. Please retry later or set ",
                    "alternativeMirror=TRUE. ATTENTION: If you use the ",
                    "alternative mirror at a time a more recent version ",
                    "of FVB or mm39 is availablde, the alternative mirror ",
                    "will pick the more recent version."
                )
            }
        } else {
            message("Connected with success.")
            return(ensembl)
        }
    }
}


retrieveInfo <- function(biomartname = "ENSEMBL_MART_ENSEMBL", #nolint
        datasetname = "mmusculus_gene_ensembl",
        versionname = "Ensembl Genes 105", alternativemirrorchoice = TRUE) {

    message("\t Connecting to biomart")
    ensembl <- tryUseMart(biomart = biomartname, dataset = datasetname,
        version = versionname, alternativeMirror = alternativemirrorchoice)
    message("\t Retrieving gene info")
    infovec <- c("chromosome_name", "start_position", "end_position",
        "external_gene_name", "strand")
    genesinfo <- tryGetBM(infovec, ensembl)
    return(list(ensembl, genesinfo))
}

################



##############
# MAIN
##############


#######
## First part that prepares data
#######

## Building gff and expnames lists
names(gffvec) <- expnamesvec

## Checking and retrieving parameters
checkparams(species_name, gffvec, expnamesvec, output_folder)

## Reading gff files
message("Reading gff files and return conversion table")
different_id_list <- lapply(gffvec, function(currentgff, dbname) {
    fi <- read.delim(currentgff, header = FALSE)
    return(clusterProfiler::bitr(fi[, 3], fromType = "ENSEMBL",
                                    toType = c("SYMBOL", "ENTREZID",
                                            "ENSEMBL", "UNIPROT",
                                            "ENSEMBLPROT"), OrgDb = dbname))
}, database_name)

## Defining background
message("Defining background")
background_id_vec <- backgrounddef(different_id_list, backgroundpath,
        species_name)

## Retrieving info from biomart
message("Retrieving info from biomart")
res <- retrieveInfo(biomartname = "ENSEMBL_MART_ENSEMBL",
        datasetname =  datasetname,
        versionname = ensemblversion, alternativemirrorchoice = TRUE)
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



#######
## Second part that compares the three categories
#######

message("Creating the input list of entrezID")
idcomplist <- lapply(different_id_list, function(x) return(unique(x$ENTREZID)))

message("Creating the entrezID-symbol table")
idtable <- do.call("rbind", different_id_list)
idx <- which(duplicated(paste(idtable$SYMBOL, idtable$ENTREZID, sep = "-")))
idtable <- idtable[-idx, ]

message("Performing clusters comparison")

!!!!!!!!!!!!!!
functionvec <- c("enrichGO")
        categoryvec <- c("CC", "BP", "MF")
        results <- gofun(functionvec, categoryvec, diffidlistforcomp,
                database_name)
tryCatchFunction <- function(id_list, f, db, currentcat, isread) { # nolint
        return(tryCatch(
                res <- compareCluster(geneClusters = id_list, fun = f,
                        OrgDb = db, ont = currentcat, readable = isread),
                error = function(e) e))
}

gofun <- function(functionvec, categoryvec, diffidlistforcomp, database_name) {
        results <- lapply(functionvec,
                function(currentfunction, catvec, diff_id_list, dbname) {
                        message("\t Performing ", currentfunction)
                        threego <- lapply(catvec, function(currentcat, f,
                                id_list, db, levelnum) {
                                message("\t\t on ", currentcat)
                                comp <- tryCatchFunction(
                                        id_list, f, db, currentcat, TRUE)
                                return(comp)
                        }, currentfunction, diff_id_list, dbname)
                        return(threego)
                }, categoryvec, diffidlistforcomp, database_name)
        results <- unlist(results, recursive = FALSE)
        names(results) <- paste0(functionvec[1], categoryvec)
        return(results)
}




!!!!!!!!!!!!!!!!!!!






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
