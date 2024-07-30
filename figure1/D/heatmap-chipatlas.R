#############################
## This script builds different heatmaps based on the chip-atlas results.
##
## Run on R 4.3.2
## Descostes June 2020 - modified dec 2023
#############################

library("pheatmap")


#############
## PARAMS
#############

resultpathvec <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG13_hist_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG13_pol_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG13_TF_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG14_hist_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG14_pol_PST.txt", #nolint
"/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/results/mHG14_TF_PST.txt") #nolint

resultnamevec <- c("mHG13hist", "mHG13pol", "mHG13TF", "mHG14hist", "mHG14pol",
    "mHG14TF")
repprefixvec <- c("mHG13", "mHG14")
suffixmerged <- c("hist", "polTFs")
experimentname <- c("peaksHG")
outputfolder <- c("/g/boulard/Projects/O-N-acetylglucosamine/analysis/chipatlas/sept2023/mouse_HG/heatmaps") #nolint
percentthreshold <- c(20)
ignoreqval <- FALSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

resultpathvec <- c("chipatlas_results/mHG1_hist_PST.txt",
"chipatlas_results/mHG1_pol_PST.txt",
"chipatlas_results/mHG1_TF_PST.txt",
"chipatlas_results/mHG2_hist_PST.txt",
"chipatlas_results/mHG2_pol_PST.txt",
"chipatlas_results/mHG2_TF_PST.txt")

resultnamevec <- c("mHG1hist", "mHG1pol", "mHG1TF", "mHG2hist", "mHG2pol",
    "mHG2TF")
repprefixvec <- c("mHG1", "mHG2")
suffixmerged <- c("hist", "polTFs")
experimentname <- c("peaksHG")
outputfolder <- c("result")
percentthreshold <- c(20)
ignoreqval <- FALSE



#############
## FUNCTIONS
#############

.checkparams <- function(outputfolder, resultpathvec, resultnamevec,
        experimentname, repprefixvec) {

    if (!file.exists(file.path(outputfolder, experimentname)))
        dir.create(file.path(outputfolder, experimentname), recursive = TRUE)

    if (!isTRUE(all.equal(length(resultpathvec), length(resultnamevec))))
        stop("One name should be given for each result file.")

    if (!isTRUE(all.equal(length(experimentname), 1)))
        stop("Experiment name should be unique.")

    if (!isTRUE(all.equal(length(outputfolder), 1)))
        stop("Output folder shoulld be unique.")

    if (!isTRUE(all.equal(length(repprefixvec), 2)))
        stop("This script was designed to handle two replicates.")
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
                    stop("All q-value are not significant. See line 73.")
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

.filteroncells <- function(resultlist, lookesc = TRUE) {

    result <- mapply(function(currentdf, currentname, lookesc) {

                if (lookesc)
                    res <- currentdf[which(currentdf$Cell == "ES cells"), ]
                else
                    res <- currentdf[which(currentdf$Cell != "ES cells"), ]

                if (isTRUE(all.equal(nrow(res), 0))) {
                    if (lookesc)
                        stop("\t No lookesc found in ", currentname)
                    else
                        message("\t No other tissues found in ", currentname)
                    return(NA)
                }else {
                    message("\t Number of results for ", currentname, ": ",
                            nrow(res))
                    return(res)
                }
            }, resultlist, names(resultlist), MoreArgs = list(lookesc),
            SIMPLIFY = FALSE)
    return(result)
}

.calculatepercentoverlap <- function(dfname, df, pthres) {

    message("\t Processing ", dfname)
    ##Calculate and order by percentage of Overlap
    nbquery <- as.numeric(unlist(lapply(strsplit(df$OverlapQuery,"/"), "[", 1)))
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
        suffvec <- c("hist", "pol", "TF")
        reslistnames <- names(reslist)

        idxres <- sapply(suffvec,
            function(currentsuff, currentpref, reslistnames){
                idx <- grep(paste0(currentpref, currentsuff), reslistnames)
                if (isTRUE(all.equal(length(idx), 0)))
                    stop("The element was not found in the list when merging",
                        " pol and TF")
                return(idx)
            }, currentpref, reslistnames)
        resnames <- c(paste0(currentpref,"hist"), paste0(currentpref, "polTFs"))
        poltfdf <- rbind(reslist[[idxres[2]]], reslist[[idxres[3]]])
        res <- list(reslist[[idxres[1]]], poltfdf)
        names(res) <- resnames
        return(res)
    }, reslist)
    return(res)
}

.buildreplacementvectors <- function() {

    sra_to_filter <- list(
        "mHG13" = list(
            "mHG13hist" = c("SRX16815323", "SRX10304465", "SRX4681571",
                "SRX2245482", "SRX1560890", "SRX1560888", "SRX13332952",
                "SRX12382446", "SRX2955747", "SRX185848", "SRX10445814",
                "SRX14210149", "SRX13332963", "SRX298193", "ERX245530",
                "SRX10424671"),
            "mHG13polTFs" = c("SRX9486380", "SRX5371695", "SRX2762171",
            "SRX5371690", "SRX11677010", "SRX404068", "SRX11221956",
            "SRX5600782", "SRX9486383", "SRX4178816", "SRX5850150",
            "SRX2486173", "SRX2486173", "SRX2486173", "SRX2486173",
            "SRX2486173", "SRX2486173", "SRX2486173", "SRX9450269",
            "SRX14663260", "SRX13439109", "SRX11677018", "SRX11221882",
            "SRX385377", "SRX385377", "SRX385377", "SRX385377", "SRX385377",
            "SRX385377", "SRX385377", "SRX385377", "SRX385377", "SRX11221945",
            "SRX4178805", "SRX9450552", "SRX11221877", "SRX10157266",
            "SRX11677013", "SRX11221901", "SRX13332957", "SRX3898501",
            "SRX2752337", "SRX2415040", "SRX213793", "SRX5926394", "SRX9450256",
            "SRX390394", "SRX10437500", "SRX3751909", "SRX482879",
            "SRX11676997", "SRX15652647", "SRX1670207", "SRX15464337",
            "SRX12622651", "SRX1688321", "SRX13945545", "SRX13439146",
            "SRX13439086", "SRX9424475", "SRX5850126", "SRX11221923",
            "SRX2545592", "SRX13332921", "SRX868181", "SRX1688323", "SRX334975",
            "SRX5850112", "SRX2388542", "SRX6071847", "SRX13440032",
            "SRX8754229", "SRX128354", "SRX9101843", "SRX191946", "SRX11525492",
            "SRX19028417", "SRX5817869", "SRX19028478", "SRX4980239",
            "SRX186548", "SRX3581849")),
        "mHG14" = list(
            "mHG14hist" = c("SRX10304465", "SRX16815323", "SRX13332952",
                "SRX4681571", "SRX1560890", "SRX1560888", "SRX2245482",
                "SRX185848", "SRX2955747", "SRX12382446", "SRX10445814",
                "SRX14210149", "SRX13332963", "SRX298193", "ERX245530",
                "SRX10424671"),
            "mHG14polTFs" = c("SRX9486380", "SRX2762171", "SRX5371695",
                "SRX11677010", "SRX19391748", "SRX9486383", "SRX5600782",
                "SRX5371688", "SRX404068", "SRX11221956", "SRX9450269",
                "SRX13439109", "SRX11677018", "SRX385377", "SRX385377",
                "SRX385377", "SRX385377", "SRX385377", "SRX385377",
                "SRX385377", "SRX385377", "SRX385377", "SRX11221882",
                "SRX14663260", "SRX11221945", "SRX2486173", "SRX2486173",
                "SRX2486173", "SRX2486173", "SRX2486173", "SRX2486173",
                "SRX2486173", "SRX4178805", "SRX10157266", "SRX11221901",
                "SRX11221877", "SRX5850150", "SRX3898501", "SRX9450552",
                "SRX13332957", "SRX5926395", "SRX2752337", "SRX390394",
                "SRX9450256", "SRX10437500", "SRX213793", "SRX2415040",
                "SRX11676997", "SRX1670207", "SRX12622651", "SRX6071873",
                "SRX15652647", "SRX3751909", "SRX13439086", "SRX1688321",
                "SRX13332921", "SRX9424476", "SRX2545592", "SRX13945545",
                "SRX11221923", "SRX334975", "SRX868181", "SRX13439146",
                "SRX5850126", "SRX2388542", "SRX6071847", "SRX15464337",
                "SRX1688324", "SRX5850112", "SRX11525492", "SRX13440032",
                "SRX17793746", "SRX9101843", "SRX5817869", "SRX128354",
                "SRX19028417", "SRX19028478", "SRX8754229")))

    perc_replace <- list(
        "mHG13" = list(
            "mHG13hist" = c(78.776291, 74.569790, 55.066922, NA, 20.65010,
                23.326960, NA, 34.416826, 41.682600, 38.814532, NA, NA, NA, NA,
                NA, NA),
            "mHG13polTFs" = c(85.468451, NA, 76.099426, NA, 68.068834, NA,
                56.40535, 64.24474, 68.068834, 71.892925, 47.227533, 72.275335,
                40.726577, 39.579350, 36.520076, 33.078394, 30.019120,
                29.827916, 56.787763, 22.944551, NA, NA, NA, 66.347992,
                63.097514, 58.508604, 53.346080, 51.051625, 47.036329,
                35.181644, 33.652008, 26.003824, NA, NA, NA, 58.126195,
                21.032505, NA, NA, NA, 49.521989, NA, NA, NA, 56.59656,
                41.873805, NA, 20.650096, NA, NA, NA, 23.900574, 26.959847,
                27.91587, rep(NA, 28))),
        "mHG14" = list(
            "mHG14hist" = c(73.556797, 75.698324, NA, 60.893855, NA, 22.998138,
                NA, 34.078212, 36.312849, 26.163873, NA, NA, NA, NA, NA, NA),
            "mHG14polTFs" = c(84.823091, 76.536313, NA, 66.945996, 75.884544,
                68.156425, 59.86965, NA, NA, 57.35568, 63.314711, NA, NA,
                65.921788, 64.990689, 59.124767, 57.728119, 54.748603,
                46.461825, 34.264432, 33.240223, 25.884544, NA, NA, NA,
                68.528864, 44.040968, 36.778399, 33.705773, 29.329609,
                28.305400, 23.74302, NA, 17.783985, NA, 58.286778, 47.299814,
                48.789572, NA, NA, 58.47300, NA, NA, 40.875233, 19.553073, NA,
                NA, NA, 29.981378, NA, NA, 21.694600, rep(NA, 14), 17.13222,
                NA, NA, NA, NA, 21.042831, NA, NA, NA, NA, NA, NA)))

    sra_replace <- list(
        "mHG13" = list(
            "mHG13hist" = c("SRX16815331", "SRX5382140", "SRX8818318", NA,
                "SRX4090625", "SRX4090627", NA, "SRX12382379", "SRX7568534",
                "SRX13341956", NA, NA, NA, NA, NA, NA),
            "mHG13polTFs" = c("SRX9486378", NA, "SRX2762151", NA, "SRX9195301",
                NA, "SRX11221932", "SRX2762152", "SRX9486379", "SRX8556273",
                "SRX2873907", "SRX15888978", "SRX147771", "SRX15888985",
                "SRX8832797", "SRX191070", "SRX9195280", "SRX5247634",
                "SRX9450269", "SRX12702093", NA, NA, NA, "SRX373167",
                "SRX11417001", "SRX11417015", "SRX11416997", "SRX11416999",
                "SRX11417014", "SRX11417011", "SRX11417005", "SRX11417007", NA,
                NA, NA, "SRX9195310", "SRX8873444", NA, NA, NA, "SRX4506780",
                NA, NA, NA, "SRX5926394", "SRX14134855", NA, "SRX8994331", NA,
                NA, NA, "SRX022697", "SRX9107990", "SRX15464336", rep(NA, 28))),
        "mHG14" = list(
            "mHG14hist" = c("SRX5382140", "SRX16815331", NA, "SRX8818318", NA,
                "SRX4090627", NA, "SRX13341956", "SRX7568534", "SRX12382379",
                NA, NA, NA, NA, NA, NA),
            "mHG14polTFs" = c("SRX9486378", "SRX2762151", NA, "SRX9195301",
                "SRX8556273", "SRX9486379", "SRX2762152", NA, NA, "SRX11221932",
                "SRX4403317", NA, NA, "SRX11417001", "SRX373167", "SRX11416997",
                "SRX11417015", "SRX11416999", "SRX11417014", "SRX11417005",
                "SRX11417011", "SRX11417007", NA, NA, NA, "SRX15888978",
                "SRX147771", "SRX15888985", "SRX8832797", "SRX9195280",
                "SRX5247634", "SRX191070", NA, "SRX8873444", NA, "SRX9195310",
                "SRX2873907", "SRX4506780", NA, NA, "SRX5926394", NA, NA,
                "SRX14134855", "SRX8994331", NA, NA, NA, "SRX9107990", NA, NA,
                "SRX022697", rep(NA, 14), "SRX15464336", NA, NA, NA, NA,
                "SRX5276471", NA, NA, NA, NA, NA, NA)))

    antigen_replace <- list(
        "mHG13" = list(
            "mHG13hist" = c("H4", "H3K4me3", "H3K4me2", NA, "H3K64ac",
                "H3K122ac", NA, "H3K36me3", "H2A.Z", "H3K9ac", NA, NA, NA, NA,
                NA, NA),
            "mHG13polTFs" = c("Nanog", NA, "Kmt2d", NA, "Tbp", NA, "Taf12",
                "Setd1a", "Sox2", "RNA polymerase II", "Yy1", "Ascl1", "FBXL10",
                "Ngn2", "Chd8", "Kdm2b", "Taf1", "Zbtb11", "Sin3a", "Smad2", NA,
                NA, NA, "Ash2L", "RPB3", "RPB10", "RPB1", "RPB2", "RPB9",
                "RPB8", "RPB5", "RPB6", NA, NA, NA, "Med1", "Banp", NA, NA, NA,
                "Brd4", NA, NA, NA, "Med24", "Nelfe", NA, "Ctcf", NA, NA, NA,
                "Nipbl", "Rad21", "Brca2", rep(NA, 28))),
        "mHG14" = list(
            "mHG14hist" = c("H3K4me3", "H4", NA, "H3K4me2", NA, "H3K122ac", NA,
                "H3K9ac", "H2A.Z", "H3K36me3", NA, NA, NA, NA, NA, NA),
            "mHG14polTFs" = c("Nanog", "Kmt2d", NA, "Tbp", "RNApolII", "Sox2",
                "Setd1a", NA, NA, "Taf12", "Sin3a", NA, NA, "RPB3", "Ash2L",
                "RPB1", "RPB10", "RPB2", "RPB9", "RPB5", "RPB8", "RPB6", NA,
                NA, NA, "Ascl1", "FBXL10", "Ngn2", "Chd8", "Taf1", "Zbtb11",
                "Kdm2b", NA, "Banp", NA, "Med1", "Yy1", "Brd4", NA, NA, "Med24",
                NA, NA, "Nelfe", "Ctcf", NA, NA, NA, "Rad21", NA, NA,
                "Nipbl", rep(NA, 14), "Brca2", NA, NA, NA, NA, "Myc", NA, NA,
                NA, NA, NA, NA)))

    return(list(sra_to_filter, perc_replace, sra_replace, antigen_replace))
}

.testrepnames <- function(namesrepsra, namesreppercreplace, namesrepsrareplace,
    namesantigenreplace) {

        if (!isTRUE(all.equal(namesrepsra, namesreppercreplace)) ||
            !isTRUE(all.equal(namesrepsra, namesrepsrareplace)) ||
            !isTRUE(all.equal(namesrepsra, namesantigenreplace)))
                stop("Problem in names of replicates")
}

.testna <- function(percreplace, srareplace, antireplace) {

    idxnaperc <- which(is.na(percreplace))
    idxnasra <- which(is.na(srareplace))
    idxnaanti <- which(is.na(antireplace))

    if (!isTRUE(all.equal(idxnaperc, idxnasra)) ||
        !isTRUE(all.equal(idxnaperc, idxnaanti)))
            stop("There is an error in NA values. They are not at",
                " the same indexes between vectors")
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

                .testna(percreplace, srareplace, antireplace)
}

.testvaluesvec <- function(sra_to_filter, perc_replace, sra_replace,
    antigen_replace) {

        # repsra <- sra_to_filter[[1]]
        # reppercreplace <- perc_replace[[1]]
        # repsrareplace <- sra_replace[[1]]
        # repantigenreplace <- antigen_replace[[1]]
        # namesrepsra <- names(sra_to_filter)[1]
        # namesreppercreplace <- names(perc_replace)[1]
        # namesrepsrareplace <- names(sra_replace)[1]
        # namesantigenreplace <- names(antigen_replace)[1]
    invisible(mapply(function(repsra, reppercreplace, repsrareplace,
        repantigenreplace, namesrepsra, namesreppercreplace, namesrepsrareplace,
        namesantigenreplace) {

            .testrepnames(namesrepsra, namesreppercreplace, namesrepsrareplace,
                namesantigenreplace)

            # sra <- repsra[[1]]
            # percreplace <- reppercreplace[[1]]
            # srareplace <- repsrareplace[[1]]
            # antireplace <- repantigenreplace[[1]]
            # namessra <- names(repsra)[1]
            # namespercreplace <- names(reppercreplace)[1]
            # namessrareplace <- names(repsrareplace)[1]
            # namesantireplace <- names(repantigenreplace)[1]
            mapply(function(sra, percreplace, srareplace, antireplace, namessra,
            namespercreplace, namessrareplace, namesantireplace) {

                message("Testing ", namespercreplace)

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

                            message("\t\t Elements: ", currentname)
                            elsratofilter <- repsrafilter[[currentname]]
                            newperc <- reppercreplace[[currentname]]
                            newsra <- repsrareplace[[currentname]]
                            newanti <- repantireplace[[currentname]]
                            restomodif <- resrep[[currentname]]

                            idx <- match(elsratofilter, restomodif$SRA)
                            if (any(is.na(idx))) {
                                idxna <- which(is.na(idx))
                                stop("Problem in retrieving sra to filter: ",
                                    elsratofilter[idxna])
                            }
                            restomodif$Antigen[idx] <- newanti
                            restomodif$PercentOverlap[idx] <- newperc
                            restomodif$SRA[idx] <- newsra

                            ## Removing lines containing NA
                            idxremove <- which(is.na(restomodif$PercentOverlap))
                            if (!isTRUE(all.equal(length(idxremove), 0)))
                                restomodif <- restomodif[-idxremove, ]
                            return(restomodif)
                    }, names(resrep), MoreArgs = list(repsrafilter,
                        reppercreplace, repsrareplace, repantireplace, resrep),
                        SIMPLIFY = FALSE)
                    return(reslist)
}

.replaceelementschipatlas <- function(sra_to_filter, perc_replace, sra_replace,
    antigen_replace, resultlistescsramerged) {

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
            antigen_replace, resultlistescsramerged), SIMPLIFY = FALSE)
    return(res)
}



.retrieveall <- function(res, sectionvec) {

    message("\t Retrieving all ")
    resdflist <- mapply(function(currentsection, currentrep, repname) {

        message("\t\t Extracting ", currentsection, " from ", repname)
        sravec <- currentrep[[currentsection]]$SRA
        antigenvec <- currentrep[[currentsection]]$Antigen
        return(data.frame(SRA = sravec,
            Antigen = antigenvec))
    }, sectionvec, res, names(res), SIMPLIFY = FALSE)
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
    }, sectionvec, res, names(res), rescomplete,
        MoreArgs = list(allsra, allantigen), SIMPLIFY = FALSE)
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
    resmerged <- merge(rescompletedlist[[1]], rescompletedlist[[2]],
        by = "Antigen")

    ## Writing the result table
    message("Writing output table result")
    outfold <- file.path(outputfolder, experimentname, sectiontype)
    if(!file.exists(outfold))
        dir.create(outfold, recursive = TRUE)
    filename <- paste0(experimentname, "-", percentthreshold,
        if (ignoreqval) "-noQval", ".txt")
    outfile <- file.path(outfold, filename)
    write.table(resmerged, file = outfile, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE)

    ## Restricting table for plotting
    mat <- cbind(resmerged$PercentOverlap.x, resmerged$PercentOverlap.y)
    rownames(mat) <- resmerged$Antigen
    colnames(mat) <- sectionvec

    ## Plotting heatmap
    .plotheatmap(mat, "none", ignoreqval, experimentname, percentthreshold,
        outputfolder, sectiontype)
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


## Retrieving embryonic stem cells and other cells
message("\n\n Retrieving embryonic stem cells")
resultlistesc <- .filteroncells(resultlist)
idxremove <- which(is.na(resultlistesc))
resultlistesc <- .removeelements(resultlistesc, idxremove)
if (!isTRUE(all.equal(length(idxremove), 0)))
    nbpeaksvec <- nbpeaksvec[-idxremove]

## Building antigen unique list
message("Building antigen unique lists for Embryonic Stem Cells:")
resultlistescsra <- .buildantigenlist(resultlistesc, thres = percentthreshold,
    keepmaxonly = TRUE)
completeresultlistescsra <- .buildantigenlist(resultlistesc)

if (!isTRUE(all.equal(names(resultlistescsra),
    names(completeresultlistescsra))))
    stop("The two lists do not contain the same elements.")

## Merging pol and TFs for each replicates
resultlistescsramerged <- .mergepolandtf(resultlistescsra, repprefixvec)
completeresultlistescsramerged <- .mergepolandtf(completeresultlistescsra,
    repprefixvec)
names(resultlistescsramerged) <- repprefixvec
names(completeresultlistescsramerged) <- repprefixvec


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
resultlistescsramerged <- .replaceelementschipatlas(sra_to_filter, perc_replace,
    sra_replace, antigen_replace, resultlistescsramerged)


####
## PART 3: Generating the heatmap
####

## Completing each element of the list with missing sra
## First retrieving all sra of the list
histsections <- paste0(repprefixvec, suffixmerged[1])
poltfsections <- paste0(repprefixvec, suffixmerged[2])
.callheatmapgeneration(resultlistescsramerged, completeresultlistescsramerged,
    histsections, "Histones", outputfolder, experimentname, percentthreshold,
    ignoreqval)
.callheatmapgeneration(resultlistescsramerged, completeresultlistescsramerged,
    poltfsections, "PolII-TFs", outputfolder, experimentname,
    percentthreshold, ignoreqval)
