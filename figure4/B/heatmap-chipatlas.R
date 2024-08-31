#############################
## This script builds different heatmaps based on the chip-atlas results.
##
## Run on R 4.3.2
## Descostes June 2020 - modified jan 2023
#############################



library("pheatmap")


#############
## PARAMS
#############


resultpathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1nodox1_histone_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1nodox1_pol_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1nodox1_TF_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1nodox2_pol_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1nodox2_TF_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1dox3_histone_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1dox3_pol_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1dox3_TF_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1dox4_histone_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1dox4_pol_all.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/results/hdld1dox4_TF_all.txt") #nolint

resultnamevec <- c("nodox1_histone", "nodox1_pol", "nodox1_TF", "nodox2_pol", "nodox2_TF", "dox3_histone", #nolint
"dox3_pol", "dox3_TF", "dox4_histone", "dox4_pol", "dox4_TF")

repprefixvec <- c("nodox1", "nodox2", "dox3", "dox4")
suffixmerged <- c("hist", "polTFs")
experimentname <- c("peaksPolII")
outputfolder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/human_polIIGlc/heatmaps") #nolint
percentthreshold <- c(20)
ignoreqval <- FALSE

!!!!!!!!!!!!!!!!!!!!!!!


resultpathvec <- c("chipatlas-results/hdld1nodox1_histone_all.txt",
"chipatlas-results/hdld1nodox1_pol_all.txt",
"chipatlas-results/hdld1nodox1_TF_all.txt",
"chipatlas-results/hdld1nodox2_pol_all.txt",
"chipatlas-results/hdld1nodox2_TF_all.txt",
"chipatlas-results/hdld1dox3_histone_all.txt",
"chipatlas-results/hdld1dox3_pol_all.txt",
"chipatlas-results/hdld1dox3_TF_all.txt",
"chipatlas-results/hdld1dox4_histone_all.txt",
"chipatlas-results/hdld1dox4_pol_all.txt",
"chipatlas-results/hdld1dox4_TF_all.txt")

resultnamevec <- c("nodox1_histone", "nodox1_pol", "nodox1_TF", "nodox2_pol",
"nodox2_TF", "dox3_histone", "dox3_pol", "dox3_TF", "dox4_histone", "dox4_pol",
"dox4_TF")

repprefixvec <- c("nodox1", "nodox2", "dox3", "dox4")
suffixmerged <- c("hist", "polTFs")
experimentname <- c("peaksPolII")
outputfolder <- c("results")
percentthreshold <- c(20)
ignoreqval <- FALSE


#############
## FUNCTIONS
#############

.checkparams <- function(outputfolder, resultpathvec, resultnamevec,
        experimentname, repprefixvec) {

            outfold <- file.path(outputfolder, experimentname)
            if (!file.exists(outfold))
                dir.create(outfold, recursive = TRUE)

            if (!isTRUE(all.equal(length(resultpathvec),
                length(resultnamevec))))
                stop("One name should be given for each result file.")

            if (!isTRUE(all.equal(length(experimentname), 1)))
                stop("Experiment name should be unique.")

            if (!isTRUE(all.equal(length(outputfolder), 1)))
                stop("Output folder shoulld be unique.")

            if (!isTRUE(all.equal(length(repprefixvec), 4)))
                stop("This script was designed to handle four replicates.")
}

.readandfilter <- function(resultpathvec, resultnamevec, chipatlascolnames,
    ignoreqval) {

        resultlist <- mapply(function(currentpath, currentname, colnamevec) {

            message("\t Processing ", currentname)
            fi <- read.delim(currentpath, stringsAsFactors = FALSE,
                header = FALSE)
            colnames(fi) <- colnamevec
            initialnb <- nrow(fi)
            nbquery <- as.numeric(unlist(lapply(strsplit(fi$OverlapQuery, "/"),
              "[", 1)))
            nbcontrol <- as.numeric(unlist(lapply(strsplit(fi$OverlapControl,
              "/"), "[", 1)))
            fi$QueryHigher <- nbquery > nbcontrol

            ## Filtering results having more overlap on random control
            if (isTRUE(all.equal(length(which(fi$QueryHigher)), 0)))
                stop("All random controls have more peaks than query.")
            fi <- fi[fi$QueryHigher, ]

            ## Filtering on the Qvalue
            if (!ignoreqval) {
                if (isTRUE(all.equal(length(which(fi$LogQval < -1.3)), 0)))
                    stop("All q-value are not significant. See line 93.")
                fi <- fi[fi$LogQval < -1.3, ]
            }

            ## Filtering the input expriments
            idxinput <- fi$AntigenClass != "Input control"
            if (isTRUE(all.equal(length(which(idxinput)), 0)))
                stop("All results are mock IP. See line 92.")
            fi <- fi[idxinput, ]

            percentkept <- round((nrow(fi) * 100) / initialnb)
            message("\t Keeping ", nrow(fi), "/", initialnb, "(", percentkept,
                "%)")
            if (!isTRUE(all.equal(percentkept, 0)))
                return(fi)
            else
                return(NA)
        }, resultpathvec, resultnamevec, MoreArgs = list(chipatlascolnames),
        SIMPLIFY = FALSE)

    return(resultlist)
}

.removeelements <- function(resultlist, idxremove) {

    if (!isTRUE(all.equal(length(idxremove), 0))) {
        message("\t Removing the following categories because of lack of ",
                "significant results: ", paste(resultnamevec[idxremove],
                        collapse = "-"))
        resultlist <- resultlist[-idxremove]
    }else {
        message("\t Nothing to remove.")
    }
    return(resultlist)
}

.filteroncells <- function(resultlist, looktypespe = TRUE) {

    result <- mapply(function(currentdf, currentname, looktypespe) {

                if (looktypespe)
                    res <- currentdf[which(currentdf$Cell == "DLD-1"), ]
                else
                    res <- currentdf[which(currentdf$Cell != "DLD-1"), ]

                if (isTRUE(all.equal(nrow(res), 0))) {
                    if (looktypespe)
                        stop("\t No looktypespe found in ", currentname)
                    else
                        message("\t No other tissues found in ", currentname)
                    return(NA)
                }else {
                    message("\t Number of results for ", currentname, ": ",
                            nrow(res))
                    return(res)
                }
            }, resultlist, names(resultlist), MoreArgs = list(looktypespe),
            SIMPLIFY = FALSE)
    return(result)
}

.calculatepercentoverlap <- function(dfname, df, pthres) {

    message("\t Processing ", dfname)
    ##Calculate and order by percentage of Overlap
    nbquery <- as.numeric(unlist(lapply(strsplit(df$OverlapQuery, "/"), "[",
        1)))
    nbpeaks <- unique(as.numeric(unlist(lapply(strsplit(df$OverlapQuery, "/"),
        "[", 2))))
    if (!isTRUE(all.equal(length(nbpeaks), 1)))
        stop("Pb with the number of query peaks.")

    percentoverlap <- (nbquery * 100) / nbpeaks
    df <- cbind(df, percentoverlap)

    idxover <- which(df$percentoverlap >= pthres)
    if (isTRUE(all.equal(length(idxover), 0)))
        stop("No lines were kept, decrease the percentage value.\n")
    else
        message("\t\t Keeping ", length(idxover), "/", nrow(df))
    df <- df[idxover, ]
    df <- df[order(df$percentoverlap), ]
    return(df)
}

.buildantigenlist <- function(resultlist, thres = 1, keepmaxonly = FALSE) {

    result <- mapply(function(df, dfname, thres) {

        df <- .calculatepercentoverlap(dfname, df, thres)
        nbquery <- as.numeric(unlist(lapply(strsplit(df$OverlapQuery, "/"),
                "[", 1)))

        if (keepmaxonly) {

            ## Keeping the line with the highest percentage of overlap by
            ## antigen
            message("\t\t\t Keeping max overlap by antigens")
            antigenlist <- split(df, df$Antigen)
            antigenlist <- lapply(antigenlist, function(currentantigen) {
                if (isTRUE(all.equal(nrow(currentantigen), 1)))
                    return(currentantigen)
                idxmax <- which.max(currentantigen$percentoverlap)
                return(currentantigen[idxmax, ])})
            antigendf <- do.call(rbind, antigenlist)
        } else {
            message("\t\t\t Computing percentages only")
            antigendf <- df
        }

        ## Keeping only the name of the antigen and the percentage
        message("\t\t Returning ", nrow(antigendf), "/", length(nbquery))
        result <- data.frame(Antigen = antigendf$Antigen,
                            PercentOverlap = antigendf$percentoverlap,
                            SRA = antigendf$ID)
        return(result)
    }, resultlist, names(resultlist), MoreArgs = list(thres),
            SIMPLIFY = FALSE)
    return(result)
}

.mergepolandtf <- function(reslist, prefvec) {

    res <- lapply(prefvec, function(currentpref, reslist) {

        message("Processing ", currentpref)
        suffvec <- c("histone", "pol", "TF")
        reslistnames <- names(reslist)

        idxres <- sapply(suffvec,
            function(currentsuff, currentpref, reslistnames) {
                idx <- grep(paste(currentpref, currentsuff, sep = "_"),
                    reslistnames)
                if (isTRUE(all.equal(length(idx), 0)) &&
                    ## nodox2 had no match for histones
                    (!isTRUE(all.equal(currentpref, "nodox2"))))
                    stop("The element was not found in the list when merging",
                        " pol and TF")

                if (isTRUE(all.equal(length(idx), 0)))
                    return(NA)
                else
                    return(idx)
            }, currentpref, reslistnames)
        resnames <- c(paste0(currentpref,"hist"), paste0(currentpref, "polTFs"))
        poltfdf <- rbind(reslist[[idxres[2]]], reslist[[idxres[3]]])
        if (isTRUE(all.equal(currentpref, "nodox2")))
            res <- list(NA, poltfdf)
        else
            res <- list(reslist[[idxres[1]]], poltfdf)
        names(res) <- resnames
        return(res)
    }, reslist)
    return(res)
}

.buildreplacementvectors <- function() {

    sra_to_filter <- list(
        "nodox1" = list(
            "nodox1hist" = c("SRX20436769"),
            "nodox1polTFs" = c("SRX7748724", "SRX7748766", "SRX20436833",
                "SRX8699645", "SRX11555395", "SRX10333914", "SRX15147204",
                "SRX11555379", "SRX15147167")),
        "nodox2" = list(
            "nodox2polTFs" = c("SRX7748724", "SRX7748766", "SRX20436833",
                "SRX8699645", "SRX11555395", "SRX10333914", "SRX15147204",
                "SRX11555379", "SRX15147167")),
        "dox3" = list(
            "dox3hist" = c("SRX20436769"),
            "dox3polTFs" = c("SRX7748724", "SRX7748766", "SRX20436833",
                "SRX13784679", "SRX8699645", "SRX11555395", "SRX10333914",
                "SRX15147204", "SRX17152983", "SRX15147167", "SRX7677759")),
        "dox4" = list(
            "dox4hist" = c("SRX20436769"),
            "dox4polTFs" = c("SRX7748724", "SRX7748766", "SRX20436833",
                "SRX8699645", "SRX11555395", "SRX10333914", "SRX15147204",
                "SRX15147167")))

    perc_replace <- list(
        "nodox1" = list(
            "nodox1hist" = c(72.000000),
            "nodox1polTFs" = c(69.536232, NA, NA, NA, NA, NA, 55.623188, NA,
                NA)),
        "nodox2" = list(
            "nodox2polTFs" = c(70.816472, NA, NA, NA, NA, NA, 56.605570, NA,
                NA)),
        "dox3" = list(
            "dox3hist" = c(79.885932),
            "dox3polTFs" = c(83.384030, NA, NA, NA, 21.863118, NA, NA,
                70.456274, NA, NA, NA)),
        "dox4" = list(
            "dox4hist" = c(74.783550),
            "dox4polTFs" = c(82.900433, NA, NA, 19.372294, NA, NA, 68.019481,
                NA)))

    sra_replace <- list(
        "nodox1" = list(
            "nodox1hist" = c("SRX10580017"),
            "nodox1polTFs" = c("SRX10580013", NA, NA, NA, NA, NA, "SRX7677719",
                NA, NA)),
        "nodox2" = list(
            "nodox2polTFs" = c("SRX10580013", NA, NA, NA, NA, NA, "SRX7677719",
                NA, NA)),
        "dox3" = list(
            "dox3hist" = c("SRX10580016"),
            "dox3polTFs" = c("SRX10580013", NA, NA, NA, "SRX8699710", NA, NA,
                "SRX7677719", NA, NA, NA)),
        "dox4"= list(
            "dox4hist" = c("SRX10580017"),
            "dox4polTFs" = c("SRX10580013", NA, NA, "SRX8699710", NA, NA,
                "SRX7677719", NA)))

    antigen_replace <- list(
        "nodox1" = list(
            "nodox1hist" = c("H3K27ac"),
            "nodox1polTFs" = c("RNA polymerase II", NA, NA, NA, NA, NA,
                "NELFCD", NA, NA)),
        "nodox2" = list(
            "nodox2polTFs" = c("RNA polymerase II", NA, NA, NA, NA, NA,
                "NELFCD", NA, NA)),
        "dox3" = list(
            "dox3hist" = c("H3K27ac"),
            "dox3polTFs" = c("RNA polymerase II", NA, NA, NA, "INTS3", NA, NA,
                "NELFCD", NA, NA, NA)),
        "dox4"= list(
            "dox4hist" = c("H3K27ac"),
            "dox4polTFs" = c("RNA polymerase II", NA, NA, "INTS3", NA, NA,
                "NELFCD", NA)))

    return(list(sra_to_filter, perc_replace, sra_replace, antigen_replace))
}

.testrepnames <- function(namesrepsra, namesreppercreplace, namesrepsrareplace,
    namesantigenreplace) {

        if (!isTRUE(all.equal(namesrepsra, namesreppercreplace)) ||
            !isTRUE(all.equal(namesrepsra, namesrepsrareplace)) ||
            !isTRUE(all.equal(namesrepsra, namesantigenreplace)))
                stop("Problem in names of replicates")
}

.testelementvec <- function(namessra, namespercreplace, namessrareplace,
    namesantigenreplace, sra, percreplace, srareplace, antireplace) {

                    if (!isTRUE(all.equal(namessra, namespercreplace)) ||
                    !isTRUE(all.equal(namessra, namessrareplace)) ||
                    !isTRUE(all.equal(length(namessra),
                        length(namesantigenreplace))))
                        stop("Problem in names of elements")

                if (!isTRUE(all.equal(length(sra), length(percreplace))) ||
                    !isTRUE(all.equal(length(sra), length(srareplace))) ||
                    !isTRUE(all.equal(length(sra), length(antireplace)))) {
                        message("sra: ", length(sra))
                        message("percreplace: ", length(percreplace))
                        message("srareplace: ", length(srareplace))
                        message("antireplace: ", length(antireplace))
                        stop("Vectors differs in length")
                }
}

.testvaluesvec <- function(sra_to_filter, perc_replace, sra_replace,
    antigen_replace) {

    invisible(mapply(function(repsra, reppercreplace, repsrareplace,
        repantigenreplace, namesrepsra, namesreppercreplace, namesrepsrareplace,
        namesantigenreplace) {

            .testrepnames(namesrepsra, namesreppercreplace, namesrepsrareplace,
                namesantigenreplace)

        mapply(function(sra, percreplace, srareplace, antireplace, namessra,
            namespercreplace, namessrareplace, namesantireplace) {

                .testelementvec(namessra, namespercreplace, namessrareplace,
                    namesantigenreplace, sra, percreplace, srareplace,
                    antireplace)

            }, repsra, reppercreplace, repsrareplace, repantigenreplace,
                names(repsra), names(reppercreplace), names(repsrareplace),
                names(repantigenreplace))
    }, sra_to_filter, perc_replace, sra_replace, antigen_replace,
        names(sra_to_filter), names(perc_replace), names(sra_replace),
        names(antigen_replace)))
}

.replaceandremovena <- function(repsrafilter, reppercreplace, repsrareplace,
    repantireplace, resrep) {

                    reslist <- mapply(function(currentname, repsrafilter,
                        reppercreplace, repsrareplace, repantireplace, resrep) {

                            if (!isTRUE(all.equal(currentname, "nodox2hist"))) {

                                message("\t\t Elements: ", currentname)
                                elsratofilter <- repsrafilter[[currentname]]
                                newperc <- reppercreplace[[currentname]]
                                newsra <- repsrareplace[[currentname]]
                                newanti <- repantireplace[[currentname]]
                                restomodif <- resrep[[currentname]]

                                idx <- match(elsratofilter, restomodif$SRA)
                                if (any(is.na(idx))) {
                                    idxna <- which(is.na(idx))
                                    stop("Problem in retrieving sra to filter:",
                                        elsratofilter[idxna])
                                }
                                restomodif$Antigen[idx] <- newanti
                                restomodif$PercentOverlap[idx] <- newperc
                                restomodif$SRA[idx] <- newsra

                                ## Removing lines containing NA
                                idxremove <- which(is.na(
                                    restomodif$PercentOverlap))
                                if (!isTRUE(all.equal(length(idxremove), 0)))
                                    restomodif <- restomodif[-idxremove, ]
                                return(restomodif)
                            } else {
                                return(NA)
                            }
                    }, names(resrep), MoreArgs = list(repsrafilter,
                        reppercreplace, repsrareplace, repantireplace, resrep),
                        SIMPLIFY = FALSE)
                    return(reslist)
}


.replaceelementschipatlas <- function(sra_to_filter, perc_replace, sra_replace,
    antigen_replace, resultlistedld1sramerged) {

        res <- mapply(function(currentrepname, sratofilter, percreplace,
            srareplace, antigenreplace, res) {

        ## Retrieving the lists of histones and polTFs for current rep
        message("\t Processing ", currentrepname)
        repsrafilter <- sratofilter[[currentrepname]]
        reppercreplace <- percreplace[[currentrepname]]
        repsrareplace <- srareplace[[currentrepname]]
        repantireplace <- antigenreplace[[currentrepname]]
        resrep <- res[[currentrepname]]

        reslist <- .replaceandremovena(repsrafilter, reppercreplace,
                repsrareplace, repantireplace, resrep)
        return(reslist)
    }, names(sra_to_filter),
        MoreArgs = list(sra_to_filter, perc_replace, sra_replace,
            antigen_replace, resultlistedld1sramerged), SIMPLIFY = FALSE)
    return(res)
}



.retrieveall <- function(res, sectionvec) {

    message("\t Retrieving all ")
    resdflist <- mapply(function(currentsection, currentrep, repname) {

        if (!isTRUE(all.equal(currentsection, "nodox2hist"))) {
            message("\t\t Extracting ", currentsection, " from ", repname)
            sravec <- currentrep[[currentsection]]$SRA
            antigenvec <- currentrep[[currentsection]]$Antigen
            return(data.frame(SRA = sravec,  Antigen = antigenvec))
        }
    }, sectionvec, res, names(res), SIMPLIFY = FALSE)
    idxnull <- which(sapply(resdflist, is.null))
    if (!isTRUE(all.equal(length(idxnull), 0)))
        resdflist <- resdflist[-idxnull]
    resdf <- do.call("rbind", resdflist)

    ## Only missing antigen will be used for completing data, following line is
    ## ok. See function .completemissingdata
    ## The sra of rep1 will be used for missing datas in rep2 and the other
    ## way around
    dupantigen <- which(duplicated(resdf$Antigen))
    if (!isTRUE(all.equal(length(dupantigen), 0)))
        resdf <- resdf[-dupantigen, ]

    return(resdf)
}

.retrievepercentvec <- function(missingsra, repcomplete) {

    return(unlist(lapply(missingsra, function(currentsra, repcomplete) {
                    idxsracomplete <- which(repcomplete$SRA == currentsra)
                    if (!isTRUE(all.equal(length(idxsracomplete), 0))) {
                        message("\t\t\t Complementary sra found")
                        return(repcomplete$PercentOverlap[idxsracomplete])
                    }else {
                        return(0)
                    }
                }, repcomplete)))
}

.completemissingdata <- function(sectionvec, res, rescomplete, allsra,
    allantigen) {

        rescompletedlist <- mapply(function(currentsection, currentrep, repname,
        currentrepcomplete, allsra, allantigen) {

            if (!isTRUE(all.equal(currentsection, "nodox2hist"))) {

                rep <- currentrep[[currentsection]]
                repcomplete <- currentrepcomplete[[currentsection]]
                idx <-  match(allantigen, rep$Antigen)
                idxna <- which(is.na(idx))

                if (!isTRUE(all.equal(length(idxna), 0))) {
                    message("\t Completing missing data for ", currentsection,
                        " from ", repname)
                    missingantigen <- allantigen[idxna]
                    missingsra <- allsra[idxna]
                    message("\t\t The missing data are ",
                        paste(missingantigen, collapse = "-"), " having sra ",
                        paste(missingsra, collapse = "-"))
                    percentvec <- .retrievepercentvec(missingsra, repcomplete)
                    missingdf <- data.frame(Antigen = missingantigen,
                                        PercentOverlap = percentvec,
                                        SRA = missingsra)
                    rep <- rbind(rep, missingdf)
                }
            return(rep)
            }
        }, sectionvec, res, names(res), rescomplete,
            MoreArgs = list(allsra, allantigen), SIMPLIFY = FALSE)
        idxnull <- which(sapply(rescompletedlist, is.null))
        if (!isTRUE(all.equal(length(idxnull), 0)))
            rescompletedlist <- rescompletedlist[-idxnull]
        return(rescompletedlist)
}

.plotheatmap <- function(mat, scalemethod, ignoreqval, experimentname,
    percentthreshold, outputfolder, sectiontype) {

        main_title <- if (!ignoreqval) experimentname else paste(experimentname,
                    "No Qval", sep = "-")
        outfile <- file.path(outputfolder, experimentname, sectiontype,
                paste0(experimentname, "-", percentthreshold,
                        if (ignoreqval) "-noQval", "-", scalemethod, ".pdf"))

        message("\t Plotting to ", outfile)

        pheatmap::pheatmap(mat,
        scale = scalemethod,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        annotation_names_row = TRUE,
        annotation_names_col = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = main_title,
        fontsize_row = 4,
        filename = outfile)
}

.callheatmapgeneration <- function(res, rescomplete, sectionvec, sectiontype,
    outputfolder, experimentname, percentthreshold, ignoreqval) {

    message("Generating heatmap for ", sectiontype)
    alldf <- .retrieveall(res, sectionvec)
    allsra <- alldf$SRA
    allantigen <- alldf$Antigen

    ## Completing missing data if necessary
    rescompletedlist <- .completemissingdata(sectionvec, res, rescomplete,
        allsra, allantigen)

    ## Merge results by antigen (assuming two replicates, verified when checking
    ## parameters)
    message("Merging replicates and create result table")
    resmerged <- do.call("cbind", rescompletedlist)
    idxantigen <- grep("Antigen", colnames(resmerged))
    if (isTRUE(all.equal(length(idxantigen), 0)))
        stop("Antigen col was not found, this should not happen")
    resmerged <- resmerged[, -idxantigen[-1]]
    colnames(resmerged)[idxantigen[1]] <- "Antigen"

    ## Writing the result table
    message("Writing output table result")
    outfold <- file.path(outputfolder, experimentname, sectiontype)
    if (!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
    filename <- paste0(experimentname, "-", percentthreshold,
        if (ignoreqval) "-noQval", ".txt")
    outfile <- file.path(outfold, filename)
    write.table(resmerged, file = outfile, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE)

    ## Creating matrix for plotting
    idxsra <- grep("SRA", colnames(resmerged))
    if (isTRUE(all.equal(length(idxsra), 0)))
        stop("SRA col was not found, this should not happen")
    resmerged <- resmerged[, -idxsra]
    colnames(resmerged) <- gsub(".PercentOverlap", "", colnames(resmerged))
    idxantigen <- grep("Antigen", colnames(resmerged))
    if (!isTRUE(all.equal(length(idxantigen), 1)))
        stop("Problem with Antigen col, this should not happen")
    antigenvec <- resmerged[, idxantigen]
    mat <- as.matrix(resmerged[, -idxantigen])
    rownames(mat) <- antigenvec

    ## Plotting heatmap
    .plotheatmap(mat, "none", ignoreqval, experimentname,
        percentthreshold, outputfolder, sectiontype)
    .plotheatmap(mat, "row", ignoreqval, experimentname, percentthreshold,
        outputfolder, sectiontype)
}



#############
## MAIN
#############

####
## PART 1: Preparing the data
####

## Checking parameters
.checkparams(outputfolder, resultpathvec, resultnamevec, experimentname,
    repprefixvec)

if (ignoreqval)
    warning("!!!!!!!! #############\n Q-val is ignored. Only matches are ",
            "reported without guarantee that they are significant.\n",
            "############# !!!!!!!!")

## Reading and filtering results files
message("Reading ChIP-Atlas results")
chipatlascolnames <- c("ID", "AntigenClass", "Antigen", "CellClass", "Cell",
        "NumOfPeaks", "OverlapQuery", "OverlapControl", "LogPVal", "LogQval",
        "FoldEnrichment")
resultlist <- .readandfilter(resultpathvec, resultnamevec, chipatlascolnames,
    ignoreqval)
names(resultlist) <- resultnamevec
idxremove <- which(is.na(resultlist))
resultlist <- .removeelements(resultlist, idxremove)
nbpeaksvec <- unlist(lapply(resultlist, function(currentcat){
                    nbpeaks <- unique(as.numeric(unlist(lapply(
                        strsplit(currentcat$OverlapQuery, "/"), "[", 2))))
                    if (!isTRUE(all.equal(length(nbpeaks), 1)))
                        stop("nbpeaks should be unique")
                    return(nbpeaks)
                }))


## Retrieving DLD-1 and other cells
message("\n\n Retrieving DLD-1 cells")
resultlistedld1 <- .filteroncells(resultlist)
idxremove <- which(is.na(resultlistedld1))
resultlistedld1 <- .removeelements(resultlistedld1, idxremove)
if (!isTRUE(all.equal(length(idxremove), 0)))
    nbpeaksvec <- nbpeaksvec[-idxremove]

## Building antigen unique list
message("Building antigen unique lists for DLD-1 Cells:")
resultlistedld1sra <- .buildantigenlist(resultlistedld1,
    thres = percentthreshold, keepmaxonly = TRUE)
completeresultlistedld1sra <- .buildantigenlist(resultlistedld1)

if (!isTRUE(all.equal(names(resultlistedld1sra),
    names(completeresultlistedld1sra))))
    stop("The two lists do not contain the same elements.")

## Merging pol and TFs for each replicates
resultlistedld1sramerged <- .mergepolandtf(resultlistedld1sra, repprefixvec)
completeresultlistedld1sramerged <- .mergepolandtf(completeresultlistedld1sra,
    repprefixvec)
names(resultlistedld1sramerged) <- repprefixvec
names(completeresultlistedld1sramerged) <- repprefixvec


####
## PART 2: Filtering the data
####

## The list of sra to filter was determined manually looking at the sra records.
##
## See the top of the script and the function .buildreplacementvectors() for
## details
resultsreplacement <- .buildreplacementvectors()
sra_to_filter <- resultsreplacement[[1]]
perc_replace <- resultsreplacement[[2]]
sra_replace <- resultsreplacement[[3]]
antigen_replace <- resultsreplacement[[4]]

## Verifying elements of each vector that must be of same length
.testvaluesvec(sra_to_filter, perc_replace, sra_replace, antigen_replace)

## Replacing values for each vector
message("Replacing values")
resultlistedld1sramerged <- .replaceelementschipatlas(sra_to_filter,
    perc_replace, sra_replace, antigen_replace, resultlistedld1sramerged)


####
## PART 3: Generating the heatmap
####

## Completing each element of the list with missing sra
## First retrieving all sra of the list
histsections <- paste0(repprefixvec, suffixmerged[1])
poltfsections <- paste0(repprefixvec, suffixmerged[2])
.callheatmapgeneration(resultlistedld1sramerged,
    completeresultlistedld1sramerged, histsections, "Histones", outputfolder,
    experimentname, percentthreshold, ignoreqval)
.callheatmapgeneration(resultlistedld1sramerged,
    completeresultlistedld1sramerged, poltfsections, "PolII-TFs", outputfolder,
    experimentname, percentthreshold, ignoreqval)
