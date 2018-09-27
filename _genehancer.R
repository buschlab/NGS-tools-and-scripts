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