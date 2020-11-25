##########################################################################################
############################### handle import and preparation of data-sets ###############
##########################################################################################

## W.I.P.!!


### ToDo:
# 1. candidate query for both 'POS' and 'START/END' cols 
# 2. gh-query -||-
# 3. fixAnnotation with adaptive 'chr' adding               DONE!
# 4. fix Annotation adaptive splitting of cols
# 5. VAF-calculation
# 6. merge with *.avinput INFO 
# 7. queryCandidates func adjust selector for 'CandidateGene' col


## 'mergeInfo' provided that 'convert2annovar.pl' as been executed with flags:
## --includeinfo --withzyg all information from original *.vcf is still included and can 
## be merged into output of 'table_annovar.pl' by this alpha version
# @param: set, annotated list that has been pre-processed with 'fixAnnotation'
# @param: info, imported *.avinput
# @return: set, merged file extended by available information cols
mergeInfo <- function( set, info ) {
    ## DEBUG ##
    #set <- test
    #info <- setINDELav
    ## DEBUG ##
    
    ### ToDo: 
    ## make this function adaptive
    ## 
    info <- dplyr::select(info, "V1","V2","V3","V4","V5","V6","V7","V8","V16","V18")
    names(info) <- c("CHROM","START","END","REF","ALT","shared","QUAL","COV","INFO1","GT:AD:DP:GQ:PL")
    set <- plyr::join(set, info, by = c("CHROM","START","END","REF","ALT"), type = "left", match = "all")
    
    return(set)
}

## 'queryCandidates' given a list of genes with START/END coordinates this function maps 
## variants in set to this list and creates a new column "CandidateGene"
# @param: set, annotated list to be matched
# @param: query, list of query-genes with cols "CandidateGene","CHROM","START","END"
# @return: set, input-set extended by "CandidateGene"-column
queryCandidates <- function( set, query ) {
	# insert new column for candidates
	idx_ColPosition <- which(names(set) %in% c("GeneName","Gene"))
	set <- data.frame(set[,1:idx_ColPosition[1]],"CandidateGene"=NA,set[,(idx_ColPosition[1]+1):ncol(set)])
	# now check for every 'POS' wether it lies in between the range of a candidate-gene
	
	### ToDo: decide if.. 
	# 1. in case of 'POS' col -> it is changed into START/END
	# 2. 'POS' is kept and subsequent analyses switch based on case?	
	
	range_Query <- dim(query)[1]
	for( i in 1:range_Query ) {
		idx_Query <- which(set$CHROM %in% query[i,"CHROM"] & 
            ((set$START <= query[i,"START"] & set$END >= query[i,"START"]) |        # 1.
            (set$START >= query[i,"START"] & set$END <= query[i,"END"]) |           # 2.
            (set$START <= query[i,"END"] & set$END >= query[i,"END"]) |             # 3.
            (set$START <= query[i,"START"] & set$END >= query[i,"END"]) ) )         # 4.
		set$CandidateGene[idx_Query] <- query$CandidateGene[i]
	}
	return(set)
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


## 'fixAnnotation' adjusts col-names of annotated lists to my convention
# @param: set, annotated list
# @return: set, fixed colnames, i.e., "CHROM,REF,ALT,POS,START,END,ID,QUAL,FILTER"; 
# removed '_refGene' suffix
fixAnnotation <- function( set ) {
    names(set)[which(grepl("\\bchr", names(set), ignore.case=TRUE))] <- "CHROM"
    names(set)[which(grepl("\\bref", names(set), ignore.case=TRUE))] <- "REF"
    names(set)[which(grepl("\\balt", names(set), ignore.case=TRUE))] <- "ALT"
    names(set)[which(grepl("\\bpos\\b", names(set), ignore.case=TRUE))] <- "POS"
    names(set)[which(grepl("\\bstart\\b", names(set), ignore.case=TRUE))] <- "START"
    names(set)[which(grepl("\\bend\\b", names(set), ignore.case=TRUE))] <- "END"
    names(set)[which(grepl("\\bid\\b", names(set), ignore.case=TRUE))] <- "ID"
    names(set)[which(grepl("\\bqual\\b", names(set), ignore.case=TRUE))] <- "QUAL"
    names(set)[which(grepl("\\bfilter\\b", names(set), ignore.case=TRUE))] <- "FILTER"
    names(set) <- gsub("_refgene\\b","",names(set), ignore.case=TRUE)
    # prepend "chr" in CHROM-col if missing
    set$CHROM[!grepl("chr", set$CHROM)] <- paste("chr", set$CHROM[!grepl("chr", set$CHROM)], sep="")
    # SPECIFIC CADD & DANN denominator adjustment
    names(set)[which(grepl("\\bcadd13_raw", names(set), ignore.case=TRUE))] <- "CADD_raw"
    names(set)[which(grepl("\\bcadd13_phred", names(set), ignore.case=TRUE))] <- "CADD_phred"
    names(set)[which(grepl("\\bdann\\b", names(set), ignore.case=TRUE))] <- "DANN"
    
	return(set)
}

## 'getHash' creates a new 'Hash' column with unique identifiers for each line
# @param: set, annotated list
# @return: set, extended list
getHash <- function( set ) {
	if( any(grepl("\\bpos\\b", names(set), ignore.case=TRUE)) ) {
	    idxHash <- c("CHROM","POS","ID","REF","ALT")
	} else {
	    idxHash <- c("CHROM","START","END","REF","ALT")
	}
	set <- cbind("Hash"=paste(set[[idxHash[1]]],set[[idxHash[2]]],set[[idxHash[3]]],
		set[[idxHash[4]]],set[[idxHash[5]]],sep="_"), set)
	set$Hash <- as.character(set$Hash)
	return(set)
}

## 'allFunctions'
fixIT <- function( set, query, GHdb, info ) {
    set <- fixAnnotation( set )
    set <- getHash( set ) 
    set <- queryCandidates( set, query)
    set <- queryGHdb( set, ghDB )
    set <- mergeInfo( set, info )
    
    return(set)
}







