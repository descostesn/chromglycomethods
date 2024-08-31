##############
# This script performs gene ontologies from gene lists with cluster profiler.
# It enables also the comparison of ontologies between different samples.
# Performed with R 4.3.2
# Descostes Feb 2020 - April 2023 - update Dec 2023
##############


suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("ggplot2"))
suppressMessages(library("biomaRt"))


################
# PARAMETERS
################


gffvec <- c("data/fig4D_genes_cluster1.gff",
"data/fig4D_genes_cluster2.gff",
"data/fig4D_genes_cluster3.gff",
"data/fig4D_genes_cluster4.gff",
"data/fig4D_genes_cluster5.gff")
expnamesvec <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")
species_name <- "human"
output_folder <- c("../results")
output_format <- "png"
database_name <- "org.Hs.eg.db"
kegg_name <- "hsa"
backgroundpath <- NA
datasetname <- "hsapiens_gene_ensembl"
ensemblversion <- "Ensembl Genes 109"
levelnum <- 3
altmirror <- TRUE

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
                stop("There is a problem of connexion to Ensembl for now. ",
                    "Please retry later or set altmirror=TRUE")
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
    return(clusterProfiler::bitr(fi[, 3], fromType = "SYMBOL",
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
        versionname = ensemblversion, alternativemirrorchoice = altmirror)
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

message("Creating the input list of entrezID") # nolint
idcomplist <- lapply(different_id_list, function(x) return(unique(x$ENTREZID)))

message("Creating the entrezID-symbol table") # nolint
idtable <- do.call("rbind", different_id_list)
idx <- which(duplicated(paste(idtable$SYMBOL, idtable$ENTREZID, sep = "-")))
idtable <- idtable[-idx, ]

message("Performing clusters comparison on molecular function") # nolint
res <- clusterProfiler::compareCluster(geneClusters = idcomplist,
            fun = "enrichGO", OrgDb = database_name, # nolint
            ont = "MF", readable = TRUE) # nolint

message("Output the dotplot of the comparison into ", output_folder, "/")
dpres <- clusterProfiler::dotplot(res, color = "p.adjust",
            showCategory = 10, font.size = 10, title = "")
ggplot2::ggsave(paste0("dotplot_top10.", output_format),
            plot = dpres, device = output_format, path = output_folder)
