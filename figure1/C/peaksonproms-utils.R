checkChromosomes <- function(fi){
    
    mainchromvec <- c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", 
            "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr3", "chr4",
            "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY")
    idxNA <- which(is.na(match(fi$V1, mainchromvec)))
    
    if(!isTRUE(all.equal(length(idxNA),0)))
        return(fi[-idxNA,])
    else
        return(fi)
}


buildGR <- function(currentPath){
	
	message("\t Processing ", currentPath)
	fi <- read.table(currentPath, stringsAsFactors=FALSE)
    fi <- checkChromosomes(fi)
	gr <- GenomicRanges::GRanges(seqnames=fi$V1,
			ranges=IRanges(start=fi$V4, end=fi$V5, names=fi$V9),
			strand=fi$V7)
	return(gr)
}

.retrieveElementsListOfGR <- function(lEx, f, f2){
	return(f(f2(lEx)))
}


.retrieveFirstExonGR <- function(firstExonList){
	
	chrVec <- unlist(lapply(firstExonList, function(x) 
						.retrieveElementsListOfGR(x, as.character, seqnames)))
	startVec <- unlist(lapply(firstExonList, function(x) 
						.retrieveElementsListOfGR(x, as.numeric, start)))
	endVec <- unlist(lapply(firstExonList, function(x) 
						.retrieveElementsListOfGR(x, as.numeric, end)))
	strandVec <- unlist(lapply(firstExonList, function(x) 
						.retrieveElementsListOfGR(x, as.character, strand)))
	exonIdVec <- unlist(lapply(firstExonList, function(x) mcols(x)$exon_id))
	firstExonGR <- GRanges(seqnames=chrVec, ranges=IRanges(start=startVec, 
					end=endVec), strand=strandVec, exonID=exonIdVec)
	
	return(firstExonGR)
}


buildRepeatsTarget <- function(txdb, repeatsList, enhancerspath){
    
    ## Retrieving Promoter, 5' UTR, 3' UTR, Exon, Intron, Downstream
    message("Retrieving genomic features")
    message("\t Retrieving  Promoter, 5' UTR, 3' UTR, Exon, Intron, ",
            "Downstream")
    promotersGR <- unique(GenomicFeatures::promoters(txdb, upstream=1000, 
				    downstream=1000))
    intronsGR <- unique(unlist(GenomicFeatures::intronsByTranscript(txdb)))
    exonsGR <- GenomicFeatures::exons(txdb, columns=NULL)
    fiveUTRsGR <- unique(unlist(GenomicFeatures::fiveUTRsByTranscript(txdb)))
    threeUTRsGR <- unique(unlist(GenomicFeatures::threeUTRsByTranscript(txdb)))
    
    ## Retrieving the first exons and subtracting it to exons
    message("\t Retrieving first exons")
    allExonsByGeneGR <- exonsBy(txdb, by = "gene")
    firstExonList <- lapply(allExonsByGeneGR, function(x) x[1,])
    firstExonGR <- .retrieveFirstExonGR(firstExonList)
    
    ## Removing firstExonGR from exonsGR
    exonsGR <- GenomicRanges::setdiff(exonsGR, firstExonGR)
    
    ## Retrieving enhancers
    enhancersGR <- buildGR(enhancerspath)
    
    ## Building list of annotations. The order of elements define the priorities
    message("Building list of annotations")
    gfNamesVec <- c("promoters", "introns", "firstExons", "otherExons", 
            "fiveUTR", "threeUTR", "enhancers")
    annotationsList <- c(repeatsList, promotersGR, intronsGR, firstExonGR, 
            exonsGR, fiveUTRsGR, threeUTRsGR, enhancersGR)
    annotationsList <- lapply(annotationsList, 
            function(.anno){mcols(.anno)<-NULL; .anno})
    names(annotationsList) <- c(repeatsNameVec, gfNamesVec)
    annotationsGRList <- GRangesList(annotationsList)
    
    ## Computing the other locations
    newAnno <- c(unlist(annotationsGRList))
    newAnno.rd <- reduce(trim(newAnno))
    otherLocations <- gaps(newAnno.rd, end=seqlengths(txdb))
    otherLocations <-  otherLocations[strand(otherLocations)!="*"]
    names(otherLocations) <- NULL
    annotationsGRList$otherLocations <- otherLocations
    
    return(annotationsGRList)
}


performOverlap <- function(annotationsGRList, queryGR){
    
    #############
    ## Performing the overlap and determining preferencce order
    ## The order of repeats matters since a peak could overlap two different
    ## repeats if they are close enough. Note the repeats in mouse do not overlap
    ## between them.
    #############
    
    message("\t Computing overlap")
    annoNamesVec <- names(annotationsGRList)
    mcols(queryGR) <- NULL
    resultOverlap <- findOverlaps(queryGR, annotationsGRList, ignore.strand=FALSE)
    idxKeep <- which(!duplicated(queryHits(resultOverlap)))
    resultOverlapPriority <- resultOverlap[idxKeep,]
    
    if(isTRUE(all.equal(length(resultOverlapPriority),0)) || 
		    length(resultOverlapPriority) > length(resultOverlap))
	    stop("Pb in the script.")
    
    if(!isTRUE(all.equal(length(resultOverlapPriority), length(queryGR))))
	    stop("All peaks of query should have a unique mapping location, pb in ",
			    "the script")
    
    return(list(resultOverlapPriority, resultOverlap, annoNamesVec))
}


performPieChart <- function(annoNamesVec, overlapPriority, pieColorVec, 
        outFold, nameQuery, percentageRepVec){
    
    message("\t Plotting piechart")
    ## Calculate number of overlaps
    subjectHitsNamesPriority <- annoNamesVec[subjectHits(overlapPriority)]
    counts <- table(subjectHitsNamesPriority)
    percentageVec <- 100 * counts / length(queryGR)
    # Line below is commented to keep all annotations in the pie chart of 
    # proportions
    #percentageRepVec <- percentageRepVec[names(percentageVec)]
    if(!isTRUE(all.equal(sum(percentageVec), 100)))
        stop("For the priority piechart, the percentages sum is not equal to ",
                "100: ", sum(percentageVec))
    
    pieColorVecHits <-  pieColorVec[names(percentageVec)]
    
    ## Plotting the priority piechart
    pdf(file=file.path(outFold, paste0(nameQuery, "priorityPiechart.pdf")), 
            width=10, height=10)
    par(mfrow=c(1,2))
    pie(percentageVec, labels = names(percentageVec), col= pieColorVecHits, 
            main=nameQuery)
    pie(percentageRepVec, labels=names(percentageRepVec), col=pieColorVec, 
            main="Proportions")
    dev.off()
    
    return(list(counts, percentageVec, subjectHitsNamesPriority))
}


performUpset <- function(annoNamesVec, overlap, nameQuery, outFold){
    
    message("Plotting upset")
    subjectHitsNames <- annoNamesVec[subjectHits(overlap)]
    regionsPerPeakList <- split(subjectHitsNames, queryHits(overlap))
    regionsList <- lapply(regionsPerPeakList, function(currentPeakRegions,
                    namesCategories){
			    isVec <- rep(FALSE, length(namesCategories))
			    names(isVec) <- namesCategories
			    currentPeakRegions <- unique(currentPeakRegions)
			    isVec[currentPeakRegions] <- TRUE
			    return(isVec)
		    }, annoNamesVec)
    regionsMatrix <- do.call(rbind, regionsList)
    
    if(!isTRUE(all.equal(nrow(regionsMatrix), length(regionsPerPeakList))) ||
		    !isTRUE(all.equal(ncol(regionsMatrix), length(annoNamesVec))))
	    stop("regionsMatrix does not have the correction dimensions.")
    
    namesColRegions <- colnames(regionsMatrix)
    res <- tibble::tibble(anno = lapply(seq_len(nrow(regionsMatrix)), 
				    function(i) namesColRegions[regionsMatrix[i,]]))
    
    g <- ggplot(res, aes_(x = ~anno)) + geom_bar() +
		    xlab(NULL) + ylab(NULL) + theme_minimal() +
		    ggupset::scale_x_upset(n_intersections = 20, order_by = "freq")
    
    ggsave(filename=paste0(nameQuery, "-priorityUpSet.pdf"), plot=g, 
            device=pdf(), path=outFold)
}


savingPeaksPerCategory <- function(overlapPriority, subjectHitsNamesPriority, 
        queryGR, outFold){
    
    message("Saving peaks per category")
    peaksIdxByCatPriorList <- split(queryHits(overlapPriority), 
            subjectHitsNamesPriority)
    
    catVec <- names(peaksIdxByCatPriorList)
    invisible(mapply(function(currentIdxVec, currentCat, queryObj,outputFold){
                        selectQueryGR <- queryObj[currentIdxVec,]
                        towrite <- data.frame(seqname=seqnames(selectQueryGR),
                                source="genomicRepartition", feature=currentCat,
                                start=start(selectQueryGR), end=end(selectQueryGR),
                                score=0, strand=strand(selectQueryGR), frame=".",
                                group=".")
                        filename <- paste0("peaks-", currentCat, ".gff")
                        write.table(towrite, 
                                file=file.path(outputFold, filename), quote=FALSE,
                                sep="\t", row.names=FALSE, col.names=FALSE)
				    }, peaksIdxByCatPriorList, catVec, 
				    MoreArgs=list(queryGR,outFold)))
    return(peaksIdxByCatPriorList)
}


## Code taken from sharedInternals.R of the CONCLUS package
tryUseMart <- function(biomart="ensembl", dataset, host, 
        alternativeMirror){
    
    c <- 1
    
    repeat{
        message("# Attempt ", c, "/5 # ",
                "Connection to Ensembl ... ")
        
        if(!alternativeMirror)
            ensembl <- try(useMart(biomart, dataset=dataset, host=host), 
                    silent=TRUE)
        else
            ensembl <- try(useEnsembl(biomart, dataset=dataset, host=host, 
                            mirror="useast"), silent=TRUE)
        
        if(isTRUE(is(ensembl, "try-error"))){
            c <- c + 1
            error_type <- attr(ensembl, "condition")
            message(error_type$message)
            
            if(c > 5)
                stop("There is a problem of connexion to Ensembl for ",
                        "now. Please retry later or set ",
                        "alternativeMirror=TRUE. ATTENTION: If you use the ", 
                        "alternative mirror at a time a more recent version ",
                        "of FVB or mm39 is availablde, the alternative mirror ",
                        "will pick the more recent version.")
        }else{
            message("Connected with success.")
            return(ensembl)
        }
    }
    
}


tryGetBM <- function(attributes, ensembl, values=NULL, filters=NULL){
    
    c <- 1
    
    repeat{
        
        message("# Attempt ", c, "/5 # ",
                "Retrieving information about genes from biomaRt ...") 
        
        
        if (is.null(values) && is.null(filters))
            res <- try(getBM(attributes=attributes, mart=ensembl), silent=TRUE)
        else
            res <- try(getBM(attributes=attributes, mart=ensembl, values=values,
                            filters=filters), silent=TRUE)
        
        if(isTRUE(is(res, "try-error"))){
            c <- c + 1
            error_type <- attr(res, "condition")
            message(error_type$message)
            
            if(c > 5)
                stop("There is a problem of connexion to Ensembl for ",
                        "now. Please retry later.")
            
        }else{
            message("Information retrieved with success.")
            return(res)
        }
    }
    
}


writePromWithUnique <- function(IDVec, startvec, endvec, symbolsTab, idxTable, 
        outFold, filename){
    
    GFF <- data.frame(seqname=symbolsTab$chromosome_name[idxTable], 
            source= symbolsTab$external_gene_name[idxTable], 
            feature= IDVec, start=startvec, end=endvec,
            score=0, strand=symbolsTab$strand[idxTable], frame=".",
            group=".")
    GFFUnique <- GFF[-which(duplicated(GFF$feature)),]
    write.table(GFF, file=file.path(outFold, paste0(filename, ".gff")), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(GFFUnique, file=file.path(outFold, 
                    paste0(filename, "-unique.gff")), quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=FALSE)
}


retrieveGeneInfo  <- function(ensembl, currentAnnoGR){
    
    attributes <- c('chromosome_name', 'ensembl_gene_id', 
            'ensembl_transcript_id_version', 'external_gene_name', 
            'start_position', 'end_position', 'strand', 'transcript_start', 
            'transcript_end', 'transcription_start_site')                    
    symbolsTab <- tryGetBM(attributes, ensembl, values=names(currentAnnoGR), 
            filters='ensembl_transcript_id_version')
    symbolsTab$strand[which(symbolsTab$strand == 1)] <- '+'
    symbolsTab$chromosome_name <- paste0("chr",symbolsTab$chromosome_name)
    if(!isTRUE(all.equal(length(which(symbolsTab$strand == -1)),0)))
        symbolsTab$strand[which(symbolsTab$strand == -1)] <- '-'
    return(symbolsTab)
}


restrictAnnoOverlap <- function(queryGR, peaksIdxByCatPriorList, currentAnnoGR){
    
    result <- findOverlaps(queryGR[peaksIdxByCatPriorList$promoters,], 
            currentAnnoGR, ignore.strand=FALSE)
    result <- result[which(!duplicated(queryHits(result))),]
    currentAnnoGR <- currentAnnoGR[subjectHits(result),]
    return(currentAnnoGR)
}


outputGFFProm <- function(annotationsGRList, queryGR, peaksIdxByCatPriorList, 
        symbolsTab, ensembl, outFold){
    
    currentAnnoGR <- annotationsGRList$promoters
    currentAnnoGR <- restrictAnnoOverlap(queryGR, peaksIdxByCatPriorList, 
            currentAnnoGR)
    symbolsTab <- retrieveGeneInfo(ensembl, currentAnnoGR)
    idxTable <- match(names(currentAnnoGR), 
            symbolsTab$ensembl_transcript_id_version)
    idxNA <- which(is.na(idxTable))
    currentAnnoGR <- currentAnnoGR[-idxNA,]
    idxTable <- idxTable[-idxNA]
    ## Transcripts
    writePromWithUnique(IDVec = symbolsTab$ensembl_transcript_id_version[idxTable],
            startvec = symbolsTab$transcript_start[idxTable],
            endvec = symbolsTab$transcript_end[idxTable], 
            symbolsTab = symbolsTab, idxTable = idxTable, outFold = outFold, 
            filename = "transcripts_fromProm")
    ## Genes
    writePromWithUnique(IDVec = symbolsTab$ensembl_gene_id[idxTable],
            startvec = symbolsTab$start_position[idxTable],
            endvec = symbolsTab$end_position[idxTable], symbolsTab = symbolsTab,
            idxTable =  idxTable, outFold = outFold, filename = "genes_fromProm")
}


uniformizeCategories <- function(catList, reference){
    return(lapply(catList, function(x,ref){
                        x <- x[names(ref)]
                        x[which(is.na(x))] <- 0
                        names(x) <- names(ref)
                        return(x)
                    }, reference))
}

#catList <- percentagesList
#y_lab <- "percent"
#outFold <- outputFolder
groupedBarplot <- function(catList, y_lab, outFold, ref, colVec){
    
    catList <- uniformizeCategories(catList,ref)
    namevec <- names(catList)
    
    if(!isTRUE(all.equal(y_lab,"Nb_of_overlap"))){
        catList <- append(catList, list(ref))
        names(catList) <- c(namevec, "proportions")
    }else
        colVec <- colVec[-length(colVec)]
    
    catMat <- do.call("rbind",catList)
    catMat <- round(catMat, digits=1)
    catDF <- melt(catMat, values="expName")
    colnames(catDF) <- c("Exp", "Transposon", "value")
    
    g <- ggplot(catDF,aes(x=Transposon,y=value,fill=factor(Exp)))+
            geom_bar(stat="identity",position="dodge")+ 
            scale_fill_manual("Exp", values = colVec)+
            xlab("Compartments")+ylab(y_lab)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()
    ggsave(filename=paste0("all_", y_lab, ".pdf"), plot=g, 
            device=pdf(), path=outFold)
}


