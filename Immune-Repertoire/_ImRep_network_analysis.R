require(stringdist)
require(igraph)
require(devtools)
require(SpiecEasi)
require(parallel)
require(igraph)
require(NetSwan)
library(grid)
library(gridExtra)
library(ggpubr)
require(data.table)


### Import all the MIXCR files in a directory
## Files should be named <Sample-ID>_<whatever>.<filetype>
## If there is a file named 'metadata' with a column 'Sample' it will be imported as well
#'@param path, STRING with absolute path to directory
#'@param rmNACol, TRUE/FALSE if columns that contain only NA values are removed
#'@return list containing all imported files indexed by sample names
importMIXCR <- function(path, rmNACol=TRUE) {
  Glist <- list()
  ## list files in directory
  input_files <- file.path(file_path, dir(file_path))
  ## reduce names to IDs - assuming that string before first underscore contains sample-IDs
  names(input_files) <- paste0(sapply(dir(file_path), function(x) unlist(strsplit(x, "_", fixed = FALSE))[1])) 
  
  if(any(grepl("metadata",names(input_files), ignore.case = TRUE))) {
    writeLines("Metadata was found!")
    Glist[["metadata"]] <- data.table::fread(input_files[grepl("metadata",names(input_files))], 
                                             data.table=FALSE)
    input_files <- input_files[!grepl("metadata",names(input_files))]
    # check if there are leading digits in the sample names
    if(any(grepl("^[[:digit:]]",names(input_files)))) {
      writeLines("Added 'S_' to Sample names because of leading digits!")
      # attach an 'S_' to the names and add a new column 'sID' to the metadata
      names(input_files) <- paste("S_",names(input_files),sep="")
      Glist[["metadata"]]$sID <- paste("S_",Glist[["metadata"]]$Sample,sep="")
    }
  } else {
    writeLines("No metadata was found!")
  }
  
  writeLines("The following samples were found and will be imported now:\n")
  print(names(input_files))
  ifelse(rmNACol, print("Removing columns that contain only NA values!"), 
         print("Columns containing only NA values will be kept!"))
  
  for(idx in names(input_files)) {
    # LOAD a single file and start processing - then add to list
    tmp_table <- data.table::fread(input_files[idx], data.table=FALSE)
    # if parameter is set - remove all columns that only contain NAs
    if(rmNACol) {
      
      idx.NAcols <- sapply(colnames(tmp_table), function(column) ifelse(all(is.na(tmp_table[,eval(column)])), FALSE, TRUE))
      tmp_table <- tmp_table[,idx.NAcols]
    }
    
    # remove leading/trailing white spaces or blank spaces in colnames 
    colnames(tmp_table) <- trimws(colnames(tmp_table), which="both")
    colnames(tmp_table) <- gsub(" ",".",colnames(tmp_table),fixed=TRUE)
    
    Glist[[idx]] <-  tmp_table 
  }
  
  writeLines("The dimension of the imported files are as follows:\n")
  import.stats <- sapply(names(Glist), function(x) dim(Glist[[x]]))
  rownames(import.stats) <- c("Rows","Columns")
  print(import.stats)
  
  return(Glist)
}


## Constructs a graph by calculating Levenshtein distance between sequences and inserting an 
# edge whenever chosen distance threshold is met.
#'@param sample.obj, a vdjtools formatted matrix of a single sample or a list of those
#'@param seqField, index or name of the column that is used for distance comparison
#'@param maxDist, maximum allowd distance/difference between sequences
#'@param group.vector, contains group denominators for list elements to set 'gID' graph attribute (default: NULL)
#'@param snames.vector, contains sample denominators for list elements to set 'sID' graph attribute (default: NULL)
#'@return igraph object (or list of those) containing sequences as nodes that are connected if distance > maxDist
getCTGraph <- function(sample.obj, seqField, maxDist=1, group.vector=NULL, snames.vector=NULL) {
  require(stringdist); require(igraph)
  if( class(sample.obj) == "data.frame" ) {
    ## add 'sID' and 'gID' column to df - if it does not exist yet
    if( !any(grepl("sID", names(sample.obj))) ) {
      sample.obj$sID <- snames.vector
    }
    if( !any(grepl("gID", names(sample.obj))) ) {
      sample.obj$gID <- group.vector
    }
    
    ## calculate levenshtein distances between strings
    dmatrix <- stringdist::stringdistmatrix(sample.obj[,eval(seqField)], 
                                            sample.obj[,eval(seqField)], method="lv")
    ## retain sequences
    colnames(dmatrix) <- sample.obj[,eval(seqField)]
    ## in case of duplicate sequences, i.e., distance == 0 we nedd a workaround
    dmatrix[(dmatrix == 0)] <- 0.5      # set all zeros to some arbitrary value that does not occur otherwise
    dmatrix[(dmatrix > maxDist)] <- 0   # set all values above threshold to 0
    diag(dmatrix) <- 0                  # reset matrix diagonal to 0 to prevent self-pointing edges
    ## set all values above maximum distance to 0
    dmatrix[(dmatrix > maxDist)] <- 0
    ## construct complete graph - because I might need it later
    g.full <- igraph::graph_from_adjacency_matrix(dmatrix, mode="undirected", weighted=TRUE)
    ## set edge weights according to number of mismatches, i.e., 1/#mismatches
    g.full <- igraph::set.edge.attribute(g.full, "weight", index=igraph::E(g.full),
                                         value=1/igraph::get.edge.attribute(g.full, "weight", 
                                                                            index=igraph::E(g.full)))
    # set shape according to intersection status
    if(any(grepl("intersectStatus",names(PFnT.unique.hyper)))) {
      sample.obj$node.shape <- sapply(sample.obj$intersectStatus, function(int.state) 
        ifelse(int.state, "circle", "square"))
      
    } else {
      # if sequence comes from multiple samples shape is circle, else square
      sample.obj$node.shape <- sapply(sample.obj$sID, function(sample.heritage)
        ifelse(length(unique(unlist(strsplit(sample.heritage, split=",")))) > 1, "circle", "square"))
      # if sequence is hyperexpanded it becomes a sphere
      if(any(colnames(sample.obj) %in% "hyperCT")) {
        sample.obj$node.shape[sample.obj$hyperCT] <- "sphere"
      }
      
    }
    
    # color based on node fraction
    cRamp <- colorRampPalette(c("#FFBB78","#D62728"))(10)  
    sample.obj$node.color <- cRamp[cut(sample.obj$cloneFraction,breaks = 10)]

    # # transfer all the vdjdata to the graph object
    # for( idx in names(sample.obj) ) {
    #   g.full <- igraph::set_vertex_attr(graph=g.full, name=eval(idx), index=igraph::V(g.full), 
    #                                     value=sample.obj[,eval(idx)])
    # } 
    # 
    # originally creating subgraph of connected components, but here we use it to add node data
    g.full <- RedeR::subg(g=g.full, dat=sample.obj, refcol=which(names(sample.obj) %in% eval(seqField)), maincomp=FALSE, connected=FALSE)
    
    # include 'sID' as graph attribute
    g.full <- igraph::set.graph.attribute(g.full, "sID", paste(
      unique(
        unlist(
          strsplit(
            igraph::get.vertex.attribute(
              graph=g.full, name="sID", index=igraph::V(g.full)), split=",",fixed=TRUE))), collapse = ","))
    # include 'gID' as graph attribute
    g.full <- igraph::set.graph.attribute(g.full, "gID", paste(
      unique(
        unlist(
          strsplit(
            igraph::get.vertex.attribute(
              graph=g.full, name="gID", index=igraph::V(g.full)), split=",",fixed=TRUE))), collapse = ","))
    ### ToDo: set edge-weigth to 1/error-rate
    #g.full <- igraph::set_edge_attr(graph=g.full, name="weight", )
    return(g.full)
    
  } else if( class(sample.obj) == "list" ) {
    ## check if group &| sname vectors have been supplied and if there same length as the sample.list
    if( !is.null(group.vector) & !(length(group.vector) == length(sample.obj)) ) {
      stop("Length of 'group.vector' and 'sample.list' is unequal!")
    }
    if( !is.null(snames.vector) & !(length(snames.vector) == length(sample.obj)) ) {
      stop("Length of 'snames.vector' and 'sample.list' is unequal!")
    }
    nw.list <- list()
    for( idx in c(1:length(sample.obj)) ) {
      ## construct AA-graph for samples
      nw.list[[names(sample.obj)[eval(idx)]]] <- getCTGraph(sample.obj[[names(sample.obj)[eval(idx)]]],
                                                            seqField,
                                                            maxDist=maxDist, group.vector[eval(idx)], 
                                                            snames.vector[eval(idx)])
    }
    return(nw.list)
  }
}


## extract Clusters of connected nodes bigger than given minimal size from a graph and 
# return new graph object
#'@param input.obj, graph (or a list of those) to work on
#'@param minSize, minimal cluster size, default 2
#'@param nodeCutOff, set number of nodes to extract, default 0 to take cluster size as metric
#'@return igraph object (or list of those) made-up of found clusters
getClusterGraph <- function(input.obj, minSize=2, nodeCutOff=0) {
  
  ## if input is a graph --> run clustering 
  if( class(input.obj) == "igraph" ) {
    print("is igraph")
    if( class(nodeCutOff) == "logical" ) {
      print("Warning: 'nodeCutOff' is a boolean - please specify numerical value for graph input")
      return(NULL)
    } else {
      if(length(igraph::degree(input.obj)) < nodeCutOff) {
        print("Warning: 'nodeCutOff' is larger than size of the input graph")
        return(NULL)
      }
    }
    ## get graphs connected components 
    tmp <- igraph::components(input.obj)
    ## construct new graph with clusters of given minimal size
    print(paste("minSize is: ",minSize))
    tmp_clusters <- which(tmp$csize >= minSize)
    print(length(tmp_clusters))
    ## if 'nodeCutOff' is 0 then take only cluster that have 'minSize', otherwise sort and take top 
    ## 'nodeCutOff' number of nodes
    if(nodeCutOff == 0) {
      print("chosen cutoff is zero")
      res <- igraph::induced_subgraph(input.obj, names(tmp$membership[which(tmp$membership %in% tmp_clusters)]))
      print(res)
    } else {
      print("chosen cutoff is:"); print(nodeCutOff)
      # 1. get all the nodes that are already included through cluster-extraction
      tmp_names <- names(tmp$membership[which(tmp$membership %in% tmp_clusters)])
      # 2. get all the nodes degrees and sort them from high to low
      tmp_degree <- sort(igraph::degree(input.obj), decreasing=TRUE)
      # 3. check if number of nodes is smaller than 'nodeCutOff'
      if(length(tmp_names) < nodeCutOff) {
        ## finish construction or continue building graph
        # 1. drop all the nodes that are already part of the clusters
        tmp_degree <- tmp_degree[which(names(tmp_degree) %in% setdiff(names(tmp_degree), tmp_names))]
        # 2. get the difference between 'nodeCutOff' and current size to fill-up remaining nodes
        tmp_names <- c(tmp_names, names(tmp_degree)[1:(nodeCutOff - length(tmp_names))])
        # 3. build graph and done.
        res <- igraph::induced_subgraph(input.obj, tmp_names)
      } else if(length(tmp_names) > nodeCutOff) {
        ## prune until it fits 'nodeCutOff'
        sdiff <- length(tmp_names) - nodeCutOff
        # 1. limit degree list to nodes that are included in clusters
        tmp_degree <- tmp_degree[which(names(tmp_degree) %in% tmp_names)]
        # 2. remove required amount of nodes from the bottom up
        tmp_names <- names(tmp_degree[1:nodeCutOff])
        # 3. build graph and done.
        res <- igraph::induced_subgraph(input.obj, tmp_names)
      } else {
        res <- igraph::induced_subgraph(input.obj, names(tmp$membership[which(tmp$membership %in% tmp_clusters)]))
      }
    }
    return(res)
    ## else if it is a list then extract graph and call yourself again
  } else if( class(input.obj) == "list" ) {
    print("is list")
    if( class(nodeCutOff) == "logical" ) {
      print("cutOff is logical")
      if( nodeCutOff == TRUE ) {
        print("is true")
        ## determine maximum cutoff value by selecting minimum length list
        nodeCutOff <- as.numeric(min(unlist(lapply(input.obj, function(x) length(igraph::degree(x))))))
        print(nodeCutOff)
      } else if( nodeCutOff == FALSE) {
        print("is false")
        nodeCutOff <- 0
        print(nodeCutOff)
      }
    }
    if( class(nodeCutOff) == "numeric" ) {
      print("is numeric - going to start list iteration now")
      cl.list <- list()
      for(idx in names(input.obj)) {
        print(paste("start: ",idx))
        cl.list[[eval(idx)]] <- getClusterGraph(input.obj[[eval(idx)]], minSize, nodeCutOff)
        print(paste(idx, " done"))
      }
      
    }
    return(cl.list)
  }
}


## join sample lists into their respective experimental groups
#'@param sample.list, a MIXCR formatted lit of samples to join
#'@param aggr.Col, either "aaSeqCDR3" or "nSeqCDR3"
#'@param group.vector string vector denoting the respective groups
#'@return list containing joined experimental groups
joinSamples <- function(sample.list, aggr.Col="nSeqCDR3", group.vector) {
  return.list <-list()
  
  for(idx.group in unique(group.vector)) {
    tmp <- NULL
    print(paste("Sample group: ", idx.group, sep=""))
    for(idx.sample in which(group.vector %in% eval(idx.group))) {
      ### DEBUG: print(names(sample.list)[idx.sample])
      tmp <- rbind.data.frame(tmp, data.frame(sample.list[[eval(idx.sample)]], 
                                              "sID"=names(sample.list)[idx.sample],
                                              "gID"=eval(idx.group)), stringsAsFactors=FALSE)
    }
    ## group equal AA sequences
    tmp$cloneId <- paste(tmp$cloneId, tmp$sID, sep="_")
    tmp <- data.table(tmp)
    if(eval(aggr.Col) == "aaSeqCDR3"){
      tmp <- tmp[, list(  cloneId = paste(cloneId, collapse=","),
                          sID = paste(sID, collapse=","),
                          gID = as.character(gID[1]),
                          hyperCT = any(hyperCT),
                          cloneClass = paste(cloneClass, collapse=","),
                          cloneCount = sum(cloneCount),
                          cloneFraction = sum(cloneFraction),
                          targetSequences = targetSequences[1],
                          targetQualities = targetQualities[1],
                          allVHitsWithScore = allVHitsWithScore[1],
                          allDHitsWithScore = allDHitsWithScore[1],
                          allJHitsWithScore = allJHitsWithScore[1],
                          allCHitsWithScore = allCHitsWithScore[1],
                          allVAlignments = allVAlignments[1],
                          allDAlignments = allDAlignments[1],
                          allJAlignments = allJAlignments[1],
                          allCAlignments = allCAlignments[1],
                          nSeqCDR3 = nSeqCDR3[1],
                          minQualCDR3 = minQualCDR3[1],
                          refPoints = refPoints[1]
      ),                               
      by = list(aaSeqCDR3)]
    } else {
      tmp <- tmp[, list(  cloneId = paste(cloneId, collapse=","),
                          sID = paste(sID, collapse=","),
                          gID = as.character(gID[1]),
                          hyperCT = any(hyperCT),
                          cloneClass = paste(cloneClass, collapse=","),
                          cloneCount = sum(cloneCount),
                          cloneFraction = sum(cloneFraction),
                          targetSequences = targetSequences[1],
                          targetQualities = targetQualities[1],
                          allVHitsWithScore = allVHitsWithScore[1],
                          allDHitsWithScore = allDHitsWithScore[1],
                          allJHitsWithScore = allJHitsWithScore[1],
                          allCHitsWithScore = allCHitsWithScore[1],
                          allVAlignments = allVAlignments[1],
                          allDAlignments = allDAlignments[1],
                          allJAlignments = allJAlignments[1],
                          allCAlignments = allCAlignments[1],
                          aaSeqCDR3 = aaSeqCDR3[1],
                          minQualCDR3 = minQualCDR3[1],
                          refPoints = refPoints[1]
      ),                               
      by = list(nSeqCDR3)]
    }
    
    # fix duplicated sID entries
    tmp$sID <- sapply(tmp$sID, function(sID.idx) paste(unique(unlist(strsplit(sID.idx, split=","))), collapse = ","))
    #tmp$hyperCT <- sapply(tmp$hyperCT, function(hyperCT.idx) any(unlist(strsplit(hyperCT.idx, split=","))))
    tmp$cloneClass <- sapply(tmp$cloneClass, function(cloneClass.idx) paste(unique(unlist(strsplit(cloneClass.idx, split=","))), collapse = ","))
    
    tmp <- tmp[with(tmp, order(-cloneCount)),]
    tmp$cloneFractionNew <- tmp$cloneCount / sum(tmp$cloneCount)
    ## transform back to data.frame
    tmp <- data.frame(tmp, stringsAsFactors=FALSE)
    ## select column in given order
    return.list[[eval(idx.group)]] <- tmp[,c("cloneId","cloneCount","cloneFraction","cloneFractionNew","targetSequences","targetQualities","allVHitsWithScore","allDHitsWithScore",
                                             "allJHitsWithScore","allCHitsWithScore","allVAlignments","allDAlignments","allJAlignments","allCAlignments","nSeqCDR3",         
                                             "minQualCDR3","aaSeqCDR3","refPoints","sID","gID","hyperCT","cloneClass")]
  }
  return(return.list)
}


### Extract hyperexpanded CTs by applying Z-score/IQR outlier detection
#'@param sample.obj, a MIXCR formatted samples or list of samples
#'@param cnt.Col, STRING name of the column to work with
#'@param method, STRING either "Z" or "IQR" or "TOP" or "PERCENT"
#'@param cutoff, threshold for Z-score method or percent I=[0,1] for PERCENT/TOP, IQR calculates its own cutoff
#'@return matrix/list with hyperexpanded CTs
selectHyperexpanded <- function(sample.obj, cnt.Col, method="Z", cutoff=2, pruneCT=TRUE) {
  
  if( class(sample.obj) == "data.frame" ) {
    if( method == "Z" ) {
      out.idx <- sapply((sample.obj[,eval(cnt.Col)] - mean(sample.obj[,eval(cnt.Col)])) / sd(sample.obj[,eval(cnt.Col)]), 
                        function(idx) ifelse(idx >= cutoff, TRUE, FALSE))
    } else if( method == "IQR" ) {
      # since we look for hyperexpanded we only need Q3 + (1.5 * IQR)
      cutoff <- quantile(sample.obj[,eval(cnt.Col)], 0.75) + (1.5 * IQR(sample.obj[,eval(cnt.Col)]))
      out.idx <- sapply(sample.obj[,eval(cnt.Col)], 
                        function(idx) ifelse(idx >= cutoff, TRUE, FALSE))
    } else if( method == "TOP" ) {
      # order just in case
      sample.obj <- sample.obj[order(sample.obj[eval(cnt.Col),])]
      # how much is xx-percent of this sample?
      cutoff <- floor(dim(sample.obj)[1] * cutoff)
      out.idx <- sapply(sample.obj[,"cloneId"], 
                        function(idx) ifelse(idx <= cutoff, TRUE, FALSE))
    } else if( method == "PERCENT" ) {
      # order just in case
      sample.obj <- sample.obj[order(sample.obj[eval(cnt.Col),])]
      
      first <- TRUE
      for (idx in 1:dim(sample.obj)[1]) {
        if(first) {
          out.idx <- sample.obj[eval(idx),eval(cnt.Col)]
          first = FALSE
        } else {
          out.idx[idx] <- out.idx[(idx-1)] + sample.obj[eval(idx),eval(cnt.Col)]
        }
      }
      out.idx[1] <- TRUE
      out.idx[2:length(out.idx)] <- sapply(out.idx[2:length(out.idx)], 
                        function(idx) ifelse(idx <= cutoff,TRUE, FALSE))
    }
    
    sample.obj$hyperCT <- out.idx
    
    if(pruneCT) {
      ret <- sample.obj[sample.obj$hyperCT,]
    } else {
      # we just return the df
      ret <- sample.obj
    }
    
    return(ret)
    
  } else if( class(sample.obj) == "list" ) {
    hyper.list <- list()
    for( idx in c(1:length(sample.obj)) ) {
      ## construct AA-graph for samples
      hyper.list[[names(sample.obj)[eval(idx)]]] <- selectHyperexpanded(sample.obj[[names(sample.obj)[eval(idx)]]],
                                                                        cnt.Col,method,cutoff,pruneCT)
    }
    return(hyper.list)
  }
}

## 'computeExpansion' adds columns clonePercent and cloneClass to single dataframes 
# or lists of MIXcR tables
#'@param input.obj, either list or single dataframe with MIXcR results
#'@param prop.col, column containing the CT proportions, default "cloneFraction"
#'@return list or dataframe with additional columns clonePercent and cloneClass
computeExpansion <- function(input.obj, prop.col="cloneFraction") {
  if( class(input.obj) == "data.frame" ) {
    input.obj$clonePercent <- input.obj[,eval(prop.col)] *100
    input.obj$cloneClass <- NA
    input.obj$cloneClass <- sapply(input.obj[,eval(prop.col)], 
                                   function(prop) ifelse( prop < 0.001, "c1",
                                                          ifelse(prop < 0.01, "c2",
                                                                 ifelse(prop < 0.1, "c3", 
                                                                        ifelse( prop < 1, "c4","c5")))))
    return(input.obj)
    
  } else if( class(input.obj) == "list" ) {
    for(list.idx in 1:length(input.obj)) {
      input.obj[[list.idx]] <- computeExpansion(input.obj[[list.idx]])
    }
    
    return(input.obj)
  }
}


# THIS ONLY RETURNS SEQUENCES FROM x.obj INPUT
## 'getIntersect' calculates distance between two sequence columns in a data.frame
# returns sequences in x that overlap/don't overlap with sequences in y (within cutoff)
#'@param x.obj, df with sequence column to work on
#'@param y.obj, df to compare against
#'@param seqField, name of sequence column, must be the same in both frames
#'@param maxDist, threshold to be considered overlap
#'@param discardOverlap, TRUE returns non-overlapping seq. in x, FALSE returns everything
#'@param names, sample naming
getIntersect <- function(x.obj, y.obj, seqField, maxDist=1, discardOverlap=TRUE, names = c("S1","S2")) {
  ## compare to itself and return every sequence that overlaps within maxDist
  require(data.table)
  
  print("Start calculating distances.")
  start.time <- Sys.time()
  dmatrix <- stringdist::stringdistmatrix(x.obj[,eval(seqField)], 
                                          y.obj[,eval(seqField)], 
                                          method="lv", useNames="string")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste("Finished distance calculation in ",time.taken, " units of time.\n", sep=""))
  
  x.obj$intersectStatus <- sapply(rownames(dmatrix), function(row.idx)  any(dmatrix[row.idx,] <= maxDist))
  
  
  
  if(discardOverlap) {
    # now check for every sequence in x if there is an overlap <= cutoff
    x.obj  <- x.obj[!(x.obj$intersectStatus),]
  } 
  
  return(x.obj)
}

## calculate comparative metrics for input graphs
#'@param input.obj, graph (or list of those) to work on
#'@param cnt.Col, STRING attribute name for counts of some sort - default: "cloneCount"
#'@param frac.Col, STRING attribute name for relative counts - default: "cloneFraction"
#'@param group.vector, vector of STRINGs containing names to attach to the output - default: NULL
#'@return list containing fields: degree, degree.dist, page_rank, betweeness, closeness
getGraphMetrics <- function(input.obj, cnt.Col="cloneCount", frac.Col="cloneFraction",group.vector=NULL) {
  ## if input is a single graph - extract metrics and return
  if( class(input.obj) == "igraph" ) {
    tmp_metric <- list()
    # CLUSTER-METRICS ---------------------------------------------------------
    tmp <- igraph::components(input.obj)
    ## 1. reads and proportions per cluster
    tmp_cmetrics <- data.frame("cindex"=c(1:length(tmp$csize)),"csize"=tmp$csize, "dummy"=NA, stringsAsFactors=FALSE)
    tmp_cmetrics$dummy <- sapply(c(1:tmp$no), function(idx.cluster) {
      idx.member <- names(tmp$membership[which(tmp$membership %in% idx.cluster)]); 
      paste(sum(igraph::get.vertex.attribute(input.obj, eval(cnt.Col),index=idx.member)),
            sum(igraph::get.vertex.attribute(input.obj, eval(frac.Col),index=idx.member)),
            paste(unique(igraph::get.vertex.attribute(input.obj, "sID",index=idx.member)), collapse=";"),
            paste(unique(igraph::get.vertex.attribute(input.obj, "gID",index=idx.member)), collapse=";"),sep="|")
    }
    )
    tmp_cmetrics <- tmp_cmetrics %>% tidyr::separate(dummy, into=c("creads","cproportion","sID","gID"), 
                                                     sep="\\|", convert=TRUE)
    ## 1.1. in case a group-denominator is given, adjust list accordingly
    if( !is.null(group.vector) ) {
      tmp_cmetrics$gID <- group.vector[1]
    }
    ## 1.2. add cluster-size to reads ratio
    tmp_cmetrics$creadsratio <- -log2(tmp_cmetrics$csize / sum(tmp_cmetrics$csize) * tmp_cmetrics$cproportion)
    
    # NETWORK-METRICS ---------------------------------------------------------
    ## 2. calculate closeness for each cluster
    tmp_nmetrics <- data.frame("sID"=igraph::get.vertex.attribute(input.obj, "sID",index=igraph::V(input.obj)),
                               cnt.Col=igraph::get.vertex.attribute(input.obj, eval(cnt.Col),index=igraph::V(input.obj)),
                               frac.Col=igraph::get.vertex.attribute(input.obj, eval(frac.Col),index=igraph::V(input.obj)),
                               "gID"=igraph::get.vertex.attribute(input.obj, "gID",index=igraph::V(input.obj)),
                               "nodeID"=names(tmp$membership), "clusterID"=tmp$membership, 
                               "clusterSize"=NA,
                               "degree"=igraph::degree(input.obj), 
                               "page_rank"=igraph::page.rank(input.obj)$vector,
                               "betweeness"=igraph::betweenness(input.obj),
                               "closeness"=NA, stringsAsFactors = FALSE)
    for( idx.cluster in c(1:tmp$no) ) {
      # get nodes in cluster
      idx.member <- names(tmp$membership[which(tmp$membership %in% idx.cluster)])
      # construct temporary cluster
      tmp_cluster <- igraph::induced_subgraph(input.obj, idx.member)
      # calculate closeness for cluster and add to list
      tmp_closeness <- igraph::closeness(tmp_cluster, mode="out")
      tmp_nmetrics[tmp_nmetrics$nodeID %in% names(tmp_closeness),"closeness"] <- tmp_closeness
      # insert clustersize as well
      tmp_nmetrics[tmp_nmetrics$clusterID == idx.cluster,"clusterSize"] <- tmp$csize[idx.cluster]
    }
    tmp_nmetrics$closeness[is.nan(tmp_nmetrics$closeness)] <- 0
    ## 2.1. in case a group-denominator is given, adjust list accordingly
    if( !is.null(group.vector) ) {
      tmp_cmetrics$gID <- group.vector[1]
    }
    
    # DEGREE-DISTRIBUTION -----------------------------------------------------
    tmp_degree_dist <- data.frame("sID"=NA,
                                  "gID"=NA,
                                  "node_degree"=NA,
                                  "frequency"=igraph::degree.distribution(input.obj))
    
    tmp_degree_dist$sID <- paste(unique(igraph::get.vertex.attribute(input.obj, "sID",index=idx.member)), collapse=";")
    tmp_degree_dist$gID <- paste(unique(igraph::get.vertex.attribute(input.obj, "gID",index=idx.member)), collapse=";")
    tmp_degree_dist$node_degree <- c(0:(length(tmp_degree_dist$frequency)-1))
    # CONSTRUCT RETURN LIST ---------------------------------------------------
    tmp_metric[["cmetrics"]] <- tmp_cmetrics
    tmp_metric[["nmetrics"]] <- tmp_nmetrics
    tmp_metric[["degree_dist"]] <- tmp_degree_dist
    return(tmp_metric)
    
    # HANDLE LIST INPUT -------------------------------------------------------
  } else if( class(input.obj) == "list" ) {
    metric.list <- list()
    for(idx in names(input.obj)) {
      le_tmp <- getGraphMetrics(input.obj[[eval(idx)]])
      metric.list[["cmetrics"]] <- rbind.data.frame(metric.list[["cmetrics"]], le_tmp[["cmetrics"]])
      metric.list[["nmetrics"]] <- rbind.data.frame(metric.list[["nmetrics"]], le_tmp[["nmetrics"]])
      metric.list[["degree_dist"]] <- rbind.data.frame(metric.list[["degree_dist"]], le_tmp[["degree_dist"]])
    }
    ### ToDo: fixed varying lengths in degree distribution list
    return(metric.list)
  }
}


## Takes a ClonoType igraph object from "getCTGraph" or "getClusterGraph", colors connected
# components and adjusts node size according to degree (log2 transformed)
#'@param igraph.obj, CT igraph object containing vertex parameters 'node.shape', 'sID', 'node.color' and edge attribute 'weight'
#'@param layoutType, switch 1:5 layout types for the graph
#'@param edge.label, name of the attribute to label edges with, default is empty ("")
plotCTNet <- function(igraph.obj, layoutType=1, edge.label="") {
  #layoutType <- 1
  
  deg <- degree(igraph.obj, mode="all")
  V(igraph.obj)$node.size <- as.numeric(factor(deg))
  gr <- cluster_edge_betweenness(igraph.obj)
  gr$com.color <- colorRampPalette(c("#AEC7E8", "whitesmoke","#9EDAE5"))(length(gr))[cut(sizes(gr),breaks = length(gr))]
  
  l <- switch(layoutType,
              layout_with_fr(igraph.obj, niter = 2000),
              layout.fruchterman.reingold(igraph.obj),
              layout_in_circle(igraph.obj),
              layout_on_sphere(igraph.obj),
              layout_randomly(igraph.obj))
  
  dnet::visNet(igraph.obj, newpage=F,
               vertex.shape=get.vertex.attribute(igraph.obj, "node.shape", V(igraph.obj)), 
               vertex.label.cex=0.55, vertex.label=V(igraph.obj)$sID,
               vertex.size = log2(V(igraph.obj)$node.size+1),
               vertex.color = get.vertex.attribute(igraph.obj, "node.color", V(igraph.obj)),
               vertex.frame.color = "black",
               vertex.label.family="serif", 
               vertex.label.dist = 0.5,
               # polygon around connected components
               mark.groups=gr, mark.shape=0.5, mark.expand=10, 
               mark.border="#C7C7C7", mark.col=gr$com.color,
               # edge parameters 
               edge.color="#7F7F7F",
               #edge.lty = E(carnival.graph)$edge.type,
               edge.arrow.size=.2,
               edge.curved=.2,
               edge.label.cex=0.5,
               edge.label.family="serif",
               edge.label.color="black",
               edge.label=get.edge.attribute(igraph.obj, name = eval(edge.label)),
               edge.width=get.edge.attribute(igraph.obj)$weight,
               center=F,glayout = l)
  
}



# EXPERIMENTAL ------------------------------------------------------------
## display network statistics
#'@param input.obj, nw.metrics or graph (or list of those) to work on
#'@return list containing plots - statistic output to come
plotNetworkMetrics <- function(input.obj, compare.groups=TRUE) {
  cols <- pals::tableau20(20)
  ## if input is a graph or a list of those, calculate metrics and re-call this function
  if( class(input.obj) == "igraph" | (class(input.obj) == "list" & class(input.obj[[1]]) == "igraph")) {
    print("List of graphs as input - calculating metrics first")
    return.list <- plotNetworkMetrics(getGraphMetrics(input.obj), compare.groups)
    ## otherwise continue plotting
  } else if( class(input.obj) == "list" & class(input.obj[[1]]) == "data.frame") {
    
    
    # TEST-CODE ---------------------------------------------------------------
    comp.select <- ifelse(compare.groups, "gID","sID")
    comp.text <- ifelse(compare.groups, "Group","Samples")
    ## ToDo: adjust colvector selection
    if(length(unique(input.obj$cmetrics[,eval(comp.select)])) <= 3) {
      c.vector <- cols[seq(1,(length(unique(input.obj$cmetrics[,eval(comp.select)]))+2),2)]
    } else {
      c.vector <- cols
    }
    
    ### CLUSTER METRICS
    tmp <- input.obj$cmetrics
    
    p_density_crr <- ggplot(tmp, aes(x = creadsratio, fill = tmp[,eval(comp.select)])) + geom_density(alpha = 0.65) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = comp.text)) + ggtitle("Distribution of cluster size - scaled to #nodes and #reads") +
      theme_bw() + scale_fill_manual(values=c.vector)
    
    my_comparisons <- combn(unique(tmp$gID), m=2, function(x) c(x), simplify = FALSE)
    
    p_violin_crr <- ggpubr::ggviolin(tmp, x = eval(comp.select), y = "creadsratio", trim = T,add.params = list(fill = "white"),
                                     fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                     ylab = "Clustersize / Reads", xlab = comp.text,
                                     title = "Cluster size to reads ratio") 
    kruskal.test( c(tmp$creadsratio), tmp[,eval(comp.select)] )
    # Add p-values comparing groups
    p_violin_crr <- p_violin_crr + #ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
    ## CONSTRUCT PLOT
    p_crr <- grid_arrange_shared_legend(p_density_crr, p_violin_crr, position="bottom")
    
    
    ### NETWORK METRICS
    tmp <- input.obj$nmetrics
    my_comparisons <- combn(unique(tmp$gID), m=2, function(x) c(x), simplify = FALSE)
    
    #my_comparisons <- list(c("skin_K289","skin_K297","skin_K300"))
    
    #min(tmp$closeness[tmp$closeness > 0]) / 2
    
    tmp$logcloseness <- NA
    tmp$logpage_rank <- NA
    tmp$logbetweeness <- NA
    tmp$logdegree <- NA
    
    for(g.idx in unique(tmp$gID)) {
      tmp$logcloseness[which(tmp$gID %in% g.idx)] <- -log2((tmp$closeness[which(tmp$gID %in% g.idx)] / (tmp$Read.proportion[which(tmp$gID %in% g.idx)]) ) + 0.0001)
      tmp$logpage_rank[which(tmp$gID %in% g.idx)] <- -log2((tmp$page_rank[which(tmp$gID %in% g.idx)] / (tmp$Read.proportion[which(tmp$gID %in% g.idx)]) ) + 0.0001)
      tmp$logbetweeness[which(tmp$gID %in% g.idx)] <- -log2((tmp$betweeness[which(tmp$gID %in% g.idx)] / (tmp$Read.proportion[which(tmp$gID %in% g.idx)]) ) + 0.0001)
      tmp$logdegree[which(tmp$gID %in% g.idx)] <- -log2((tmp$degree[which(tmp$gID %in% g.idx)] / (tmp$Read.proportion[which(tmp$gID %in% g.idx)]) ) + 0.0001)
      
    }
    
    ## closeness
    p_density_lcl <- ggplot(tmp, aes(x = logcloseness, fill = tmp[,eval(comp.select)])) + geom_density(alpha = 0.65) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = comp.text)) + ggtitle("Distribution of cluster size - scaled to #nodes and #reads") +
      theme_bw() + scale_fill_manual(values=c.vector)
    ## statistical tests on graph metrics
    p_violin_lcl <- ggpubr::ggviolin(tmp, x = eval(comp.select), y = "logcloseness", trim = T,add.params = list(fill = "white"),
                                     fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                     ylab = "Clustersize / Reads", xlab = comp.text,
                                     title = "Log2 closeness") 
    kruskal.test( c(tmp$logcloseness), tmp[,eval(comp.select)] )
    # Add p-values comparing groups
    p_violin_lcl <- p_violin_lcl + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
    ## CONSTRUCT PLOT
    p_lcl <- grid_arrange_shared_legend(p_density_lcl, p_violin_lcl, position="bottom")
    
    ## page_rank
    p_density_lpr <- ggplot(tmp, aes(x = logpage_rank, fill = tmp[,eval(comp.select)])) + geom_density(alpha = 0.65) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = comp.text)) + ggtitle("Distribution of page-rank - scaled to read proportion") +
      theme_bw() + scale_fill_manual(values=c.vector)
    ## statistical tests on graph metrics
    p_violin_lpr <- ggpubr::ggviolin(tmp, x = eval(comp.select), y = "logpage_rank", trim = T,add.params = list(fill = "white"),
                                     fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                     ylab = "Clustersize / Reads", xlab = comp.text,
                                     title = "Log2 page-rank") 
    kruskal.test( c(tmp$logpage_rank), tmp[,eval(comp.select)] )
    # Add p-values comparing groups
    p_violin_lpr <- p_violin_lpr + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
    ## CONSTRUCT PLOT
    p_lpr <- grid_arrange_shared_legend(p_density_lpr, p_violin_lpr, position="bottom")
    
    ## degree
    p_density_ldeg <- ggplot(tmp, aes(x = logdegree, fill = tmp[,eval(comp.select)])) + geom_density(alpha = 0.65) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = comp.text)) + ggtitle("Distribution of node degree") +
      theme_bw() + scale_fill_manual(values=c.vector)
    ## statistical tests on graph metrics
    kruskal.test( c(tmp$logdegree), tmp[,eval(comp.select)] )
    p_violin_ldeg <- ggpubr::ggviolin(tmp, x = eval(comp.select), y = "logdegree", trim = T,add.params = list(fill = "white"),
                                      fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                      ylab = "Clustersize / Reads", xlab = comp.text,
                                      title = "Log2 node degree") 
    kruskal.test( c(tmp$logdegree), tmp[,eval(comp.select)] )
    # Add p-values comparing groups
    p_violin_ldeg <- p_violin_ldeg + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
    ## CONSTRUCT PLOT
    p_ldeg <- grid_arrange_shared_legend(p_density_ldeg, p_violin_ldeg, position="bottom")
    
    ## betweeness
    p_density_lbet <- ggplot(tmp, aes(x = logbetweeness, fill = tmp[,eval(comp.select)])) + geom_density(alpha = 0.65) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = comp.text)) + ggtitle("Distribution of node betweeness") +
      theme_bw() + scale_fill_manual(values=c.vector)
    ## statistical tests on graph metrics
    kruskal.test( c(tmp$logbetweeness), tmp[,eval(comp.select)] )
    p_violin_lbet <- ggpubr::ggviolin(tmp, x = eval(comp.select), y = "logbetweeness", trim = T,add.params = list(fill = "white"),
                                      fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                      ylab = "Clustersize / Reads", xlab = comp.text,
                                      title = "Log2 betweeness") 
    kruskal.test( c(tmp$logbetweeness), tmp[,eval(comp.select)] )
    # Add p-values comparing groups
    p_violin_lbet <- p_violin_lbet + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
    ## CONSTRUCT PLOT
    p_lbet <- grid_arrange_shared_legend(p_density_lbet, p_violin_lbet, position="bottom")
    
    
    ### DEG>rEE-DISTRIBUTION
    tmp <- input.obj$degree_dist
    tmp$logfrequency <- NA
    for(g.idx in unique(tmp$gID)) {
      tmp$logfrequency[which(tmp$gID %in% g.idx)] <- -log2((tmp$frequency[which(tmp$gID %in% g.idx)]) + 0.0001 )
    }
    
    ## betweeness
    p_density_degdist <- ggplot(tmp, aes(x = logfrequency, fill = tmp[,eval(comp.select)])) + geom_density(alpha = 0.65) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = comp.text)) + ggtitle("Distribution of node-degree") +
      theme_bw() + scale_fill_manual(values=c.vector)
    ## statistical tests on graph metrics
    kruskal.test( c(tmp$logfrequency), tmp[,eval(comp.select)] )
    p_violin_degdist <- ggpubr::ggviolin(tmp, x = eval(comp.select), y = "logfrequency", trim = T,add.params = list(fill = "white"),
                                         fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                         ylab = "Clustersize / Reads", xlab = comp.text,
                                         title = "Log2 node degree distribution") 
    kruskal.test( c(tmp$logfrequency), tmp[,eval(comp.select)] )
    # Add p-values comparing groups
    p_violin_degdist <- p_violin_degdist + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
    ## CONSTRUCT PLOT
    p_degdist <- grid_arrange_shared_legend(p_density_degdist, p_violin_degdist, position="bottom")
    
    
    
    
    
    
    # STATISTICAL TESTING ------------------------------------------------------
    ### ToDo: return list with results for statistical tests
    
    #kruskal.test( c(plot_df$clusterSize), plot_df$group )
    #plot_df$Deg2Read <- plot_df$degree / plot_df$Read.count
    
    
    #ggplot(plot_df, aes(x=logcloseness, fill=group)) + geom_density(alpha=.3)
    
    #ggplot(plot_df, aes(x=degree, fill=group)) +
    #geom_histogram(binwidth=.5, position="dodge")
    
    #ggplot(plot_df, aes(x=page_rank, fill=group)) + geom_density(alpha=.3)
    
    # p_dd_sample <- ggplot(plot_df, aes(x=degree, y=page_rank, colour=sID)) + 
    #   geom_line() + 
    #   geom_point(size=2.5) + 
    #   xlim(0,max(degdist_df$node_degree)) + 
    #   #scale_color_manual(values=cols[1:length(unique(degdist_df$sID))]) + 
    #   theme_bw()
    # 
    # p_dd_group <- ggplot(degdist_df, aes(x=node_degree, y=frequency, colour=group)) + 
    #   geom_line() + 
    #   geom_point(size=2.5) + 
    #   xlim(0,max(degdist_df$node_degree)) + 
    #   scale_color_manual(values=cols[seq(1,(length(unique(group.vector))+1),2)]) + 
    #   theme_bw()
    
    return.list <- list()
    return.list[["crr"]] <- p_crr
    return.list[["lcl"]] <- p_lcl
    return.list[["lpr"]] <- p_lpr
    return.list[["ldeg"]] <- p_ldeg
    return.list[["lbet"]] <- p_lbet
    return.list[["degdist"]] <- p_degdist
  }
  
  pdf(file=paste("Network.metrics",comp.text,".pdf",sep=""), width=18, height=9)
    lapply(return.list, function(x) {grid.newpage(); grid.draw(x)})
  dev.off()
  return(return.list)
}


### ToDo: a lot
plotRobustness <- function(fobj, ptitle) {
  pdf(file=paste("data/plots/","Network_Properties_Robustness.",ptitle,".pdf",sep=""), width=12, height=6)
  
  plot(fobj[,1],fobj[,5], type='o', col='yellow',xlab="Fraction of nodes removed",
       ylab="Connectivity loss", main=paste("Robustness Analysis - ", ptitle, sep=""))
  lines(fobj[,1],fobj[,3], type='o', col='red')
  lines(fobj[,1],fobj[,4], type='o', col='orange')
  lines(fobj[,1],fobj[,2], type='o', col='blue')
  legend('bottomright',c("Random", "Betweenness", "Degree", "Cascading"),
         lty=c(1,1,1,1), pch=c(1,1,1,1),
         col=c("yellow","blue","red", "orange"))
  dev.off() 
}





















