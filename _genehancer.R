### Use the GenehancerDB to analyse promoter/enhancer regions in NGS data according to a
### given set of query genes
### Download 'genehancer.xlsx' from https://genecards.weizmann.ac.il/geneloc/index.shtml
### Use supplied functions to import, match to genes of interest and map to NGS data





## 'importGH' imports the datadump.xlsx from 
## https://genecards.weizmann.ac.il/geneloc/index.shtml and returns adjusted ghDB object
# @param: path, filepath on your system
# @return: ghDB, genehancer-DB import to work with
importGH <- function( path ) {
    genehancer <- openxlsx::read.xlsx(path, cols=c(1:6,9),
		na.strings = c("DO","n.g.","ND","/","NA"," ",".",""))
    names(genehancer) <- c("CHROM","source","type","START","END","score","attributes")
    return(genehancer)
}

## 'evalGH' analyses cols 'CandidateGene' and 'GHID' and returns summary for input-set
# @param: input, annotated list (SV,CNV,INDEL,SNV)
# @param: ghDB, genehancer-DB object 
# @return: summary object with cols 'Gene','scores','sum|max|minscores','nscores'
evalGH <- function( input, ghDB ) {
    # 1. get the GHIDs
    ghIDs <- unlist(strsplit(input$GHID[!is.na(input$GHID)], split=";"))
    # 2. get the unique ones and sort them by occurrence 
    uniqueOnes <- sort(table(ghIDs), decreasing=TRUE)
    # 3. extract associated genes of interest from ghDB
    goi <- unlist(strsplit(ghDB$goi[ghDB$GHID %in% names(uniqueOnes)], ";"))
    # 4. create data.frame with affected genes
    df_goi <- data.frame("Gene"=unique(sapply(goi, function(x) 
        unlist(strsplit(x, ","))[1])), "scores"=NA, "sumscores"=NA, "maxscores"=NA, 
        "minscores"=NA, "nscores"=NA, stringsAsFactors=FALSE)
    # 5. fill up the cols
    # 5.1 'scores'
    #bla <- sapply(df_goi$Gene, function(y) which(sapply(strsplit(goi, ","), function(z) z[1]) %in% y))
    #df_goi$scores <- sapply(bla, function(w) paste0(sapply(strsplit(goi[w], ","), function(x) x[2]), collapse=","))
    df_goi$scores <- sapply(sapply(df_goi$Gene, function(y) which(sapply(strsplit(goi, ","), 
        function(z) z[1]) %in% y)), function(w) paste0(sapply(strsplit(goi[w], ","), 
        function(x) x[2]), collapse=","))
    # 5.2 'sumscores'
    df_goi$sumscores <- sapply(strsplit(df_goi$scores, ","), function(x) sum(as.numeric(x)))
    df_goi$maxscores <- sapply(strsplit(df_goi$scores, ","), function(x) max(as.numeric(x)))
    df_goi$minscores <- sapply(strsplit(df_goi$scores, ","), function(x) min(as.numeric(x)))
    df_goi$nscores <- sapply(strsplit(df_goi$scores, ","), function(x) length(x))
    df_goi$sdscores <- sapply(df_goi$scores, function(x) sd(as.numeric(unlist(strsplit(x, ",")))))
    # replace NAs with zeros
    df_goi$sdscores[df_goi$sdscores %in% NA] <- 0

    return(df_goi)
}

## 'topGGH' takes an evaluation object as input and returns a bar-plot object of the top
# genes of interest associated with GH-regions, i.e., most GH-regions
# @param: evalObj, output of 'evalGH' function
# @param: percent, displayed percentage of genes [0,1], defaults to 10% if left empty
# @param: plottitle, headline for this plot, may be left empty
# @return: p1, a ggplot2 prepped plot
topGGH <- function( evalObj, percent, plottitle ) {
    # calculation
    if( missing(percent) ) {
        percent <- 0.1
    }
    tmpCutoff <- quantile(abs(evalObj$nscores),(1 - percent))
    idx <- which(evalObj$nscores > tmpCutoff)
    df_input <- evalObj[idx,]
    # information
    if( missing(plottitle) ) {
        plottitle <- ""
    }
    annotation <- "" # note this could be adjusted to carry additional information
    xlabel <- paste0("Top ", 100*percent,"% genes of interest", collapse="")
    # plotation
    p1 <- ggplot2::ggplot(data=df_input, ggplot2::aes(x=Gene, y=nscores, fill=Gene)) + ggplot2::guides(fill=FALSE) +
 	    ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge(), alpha=0.5) +
 	    ggplot2::labs(title=plottitle, x=xlabel, y="Number of associated GH-regions") +
 	    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 75, hjust = 1)) +
 	    ggplot2::annotate(geom="text", Inf, Inf, hjust=1, vjust=1, label = annotation)
 	    
    return(p1)
}

## 'topGScore' takes an evaluation object as input and returns a violin-plot object of the top
# genes of interest based on max-score
# @param: evalObj, output of 'evalGH' function
# @param: percent, displayed percentage of genes [0,1]1], defaults to 10% if left empty
# @param: plottitle, headline for this plot, may be left empty
# @return: p1, a ggplot2 prepped plot
topGScore <- function( evalObj, percent, plottitle ) {
    ## DEVEL ##
    #percent <- 0.1
    #evalObj <- snvEV
    #plottitle <- "development"    
    ## DEVEL ##

    #ToDo: 
    # 1.) mark genes that overlap with association plot
    # calculation
    if( missing(percent) ) {
        percent <- 0.1
    }
    tmpCutoff <- quantile(abs(evalObj$sumscores),(1 - percent))
    idx <- which(evalObj$sumscores > tmpCutoff)
    df_input <- data.frame(evalObj[idx,], minimax=NA)
    ## now prep. coloring
    # define used colors
    cmax <- "#E42032"
    cmin <- "#005877"
    cnormal <- "#000000"
    cdist <- "#EC7404"
    cmean <- "#00AEC7"
    # take 10% max scores
    tmpCutoff <- quantile(abs(df_input$sumscores),(0.9))
    df_input$minimax <- ifelse(df_input$sumscores >= tmpCutoff, 'max. sum',
                                   ifelse(df_input$sumscores %in% min(df_input$sumscores), 
                                         'min. sum', 'in between'))    
    df_input <- tidyr::separate_rows(df_input, scores, sep=",", convert=TRUE)
    df_input <- data.frame(df_input, "logscores"=(log10(1 + muh$scores)), stringsAsFactors=FALSE)
    df_input$Gene <- as.factor(muh$Gene)

    # information
    if( missing(plottitle) ) {
        plottitle <- ""
    }
    annotation <- "" # note this could be adjusted to carry additional information
    xlabel <- paste0("Top ", 100*percent,"% genes of interest", collapse="")
    # plotation
    p1 <- ggplot2::ggplot(df_input, aes(x=Gene, y=logscores, color=minimax, fill="score distribution")) + 
        ggplot2::geom_violin(alpha=0.5) + scale_color_manual(values=c(cnormal, cmax, cmin)) + scale_fill_manual(values="black") + theme_bw() + 
        #stat_summary(fun.y=mean, geom="point", shape=23, size=2) + # +mean points
        #stat_summary(fun.y=median, geom="line", size=2, color="red") + # +median points in red
        stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar",width = 0.5, alpha=0.5, color=cmean) + 
        geom_jitter(shape=1, position=position_jitter(0.2), color=cdist, size=1, alpha=0.5, show.legend=TRUE) + 
 	    ggplot2::labs(title=plottitle, x=xlabel, y="Sum of GH-scores") +
 	    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 75, hjust = 1)) +
 	    ggplot2::annotate(geom="text", Inf, Inf, hjust=1, vjust=1, label = annotation)
 	    
    return(p1)
}

## 'topGDA' analyses cols 'CandidateGene' return plot of directly associated Cgenes
# @param: input, annotated list (SV,CNV,INDEL,SNV)
# @param: percent, displayed percentage of genes [0,1]1], defaults to 10% if left empty
# @param: plottitle, headline for this plot, may be left empty
# @return: p1, a ggplot2 prepped plot
topGDA <- function( input, percent, plottitle ) {
    ## top 25% goi matched to SV/CNV/SNV/INDEL.. i.e. direct association
    df_dagoi <- data.frame("Gene"=unique(input$CandidateGene[!is.na(input$CandidateGene)]),
    "n"=sapply(unique(input$CandidateGene[!is.na(input$CandidateGene)]), function(x) length(which(input$CandidateGene %in% x))))
    if( missing(percent) ) {
        percent <- 0.1
    }
    tmpCutoff <- quantile(abs(df_dagoi$n),(1 - percent))
    idx <- which(df_dagoi$n > tmpCutoff)
    df_input <- df_dagoi[idx,]
    # information
    if( missing(plottitle) ) {
        plottitle <- ""
    }
    annotation <- "" # note this could be adjusted to carry additional information
    xlabel <- paste0("Top ", 100*percent,"% genes of interest", collapse="")
    # plotation
    p1 <- ggplot2::ggplot(data=df_input, ggplot2::aes(x=Gene, y=n, fill=Gene)) + ggplot2::guides(fill=FALSE) +
        ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge(), alpha=0.5) +
 	    ggplot2::labs(title=plottitle, x=xlabel, y="Number of directly matched regions") +
 	    ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 75, hjust = 1)) +
 	    ggplot2::annotate(geom="text", Inf, Inf, hjust=1, vjust=1, label = annotation)

    return(p1)
}





## 'matchGenehancers' given a list of names for gois this function finds all the regions
## in ghDB that are associated with those genes and outputs a reduced version of ghDB to
## work with in your project
# @param: query, list of gene names to match against
# @param: genehancer, imported dataframe of ghDB
# @return: tmpGH, ghDB reduced to goi
matchGenehancers <- function( query, genehancer ) {
	tmp <- genehancer$attributes
	idxGH <- NULL
	
	for( idxQuery in query ) {
	    ##print(idxQuery)                                         # print out were we are
	    idxTmp <- which(grepl(paste("=",idxQuery,";",sep=""), tmp, ignore.case=TRUE)) # find matches in GHdb
	    idxGH <- union(idxGH, idxTmp)	                        # join indicies together
	}
	print("Finished GH to query mapping! Proceed refinement.")
	tmpGH <- genehancer[idxGH,c("CHROM","START","END","type","score","attributes")]
	tmpGH <- data.frame("GHID"=NA, tmpGH, "goi"=NA)

    ## extract genehancerIDs from 'attributes' col and insert into 'GHID' col
    # OLD CODE IS OLD
    #idx_CG <- names(which(sapply(unlist(strsplit(tmpGH$attributes, ";")),
	#    function(x)any(grepl("genehancer_id", x, ignore.case=TRUE)))))	
	#ID <- strsplit(idx_CG, "=")
	#ID <- sapply(ID,function(x) x[2])
	## in a single fancy cmd =D
	tmpGH$GHID <- sapply(strsplit(names(which(sapply(unlist(strsplit(tmpGH$attributes, 
	    ";")),function(x)any(grepl("genehancer_id", x, ignore.case=TRUE))))), "="),
	    function(x) x[2])
	## also fill 'goi' col with associated values for easier post-processing
	for( idxQuery2 in query ) {
	    #print(idxQuery2)
	    idx_CG <- which(sapply(strsplit(tmpGH$attributes, ";connected_gene="),function(x)any(grepl(paste("^",idxQuery2,";",sep=""), x, ignore.case=TRUE))))
        # catch empty-list error
        if( length(idx_CG) > 0 ) {
            names(idx_CG) <- sapply(strsplit(names(which(sapply(unlist(strsplit(tmpGH$attributes[idx_CG], ";connected_gene=")),function(x)any(grepl(paste("^",idxQuery2,";",sep=""), x, ignore.case=TRUE))))), ";score="),function(x) paste(x[1],x[2],sep=","))
            tmpGH$goi[idx_CG] <- paste(tmpGH$goi[idx_CG], names(idx_CG), sep=";")
        }
	}
	## remove 'NA;' string in front
	tmpGH$goi <- sub("NA;","",tmpGH$goi)	
    ## return GHdb reduced to match queried genes
    return(tmpGH)
}

## 'queryGHdb' given a list of GH-regions with START/END coordinates this function maps 
## variants in set to this list and creates a new column "GHID"
# @param: set, annotated list to be matched
# @param: GHdb, list of GH-regions with cols "GHID","CHROM","START","END"
# @return: set, input-set extended by "CandidateGene"-column
queryGHdb <- function( set, GHdb ) {
	# insert new column for candidates
	idx_ColPosition <- which(names(set) %in% c("CandidateGene"))
	set <- data.frame(set[,1:idx_ColPosition[1]],"GHID"=NA,set[,(idx_ColPosition[1]+1):ncol(set)])
	# now check for every 'POS' wether it lies in between the range of a candidate-gene
	
	### ToDo: decide if.. 
	# 1. in case of 'POS' col -> it is changed into START/END
	# 2. 'POS' is kept and subsequent analyses switch based on case?	
	
	range_Query <- dim(GHdb)[1]
	for( i in 1:range_Query ) {
		idx_Query <- which(set$CHROM %in% GHdb[i,"CHROM"] & 
            ((set$START <= GHdb[i,"START"] & set$END >= GHdb[i,"START"]) |        # 1.
            (set$START >= GHdb[i,"START"] & set$END <= GHdb[i,"END"]) |           # 2.
            (set$START <= GHdb[i,"END"] & set$END >= GHdb[i,"END"]) |             # 3.
            (set$START <= GHdb[i,"START"] & set$END >= GHdb[i,"END"]) ) )         # 4.
		set$GHID[idx_Query] <- paste(set$GHID[idx_Query], GHdb$GHID[i], sep=";")
	}
	## remove 'NA;' string in front
	set$GHID <- sub("NA;","",set$GHID)	
	return(set)
} 