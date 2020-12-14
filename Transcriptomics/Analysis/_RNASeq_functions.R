# INITIALIZATION ----------------------------------------------------------

# This is a nice color palette.
cols <- cols <- pals::tableau20(20)
# Take a look at it.
scales::show_col(cols)
# Create a variable to our working directory for CARNIVAL - Cplex output
WORKDIR <- getwd()

# We also want to have a certain directory structure to store output and saved
# data objects. For me this is 1. data - 1.1 output; 1.2 RDS Show warnings is 
# disabled so that we are not bothered if the directories already exist.
dir.create("data", showWarnings = FALSE)
dir.create("data/output", showWarnings = FALSE)
dir.create("data/RDS", showWarnings = FALSE)
dir.create("data/RAW", showWarnings = FALSE)


# REQUIRED PACKAGES -------------------------------------------------------

# We will need these libraries. Check out 'install.packages.R' if you want
# some help with the installation.
library(tidyverse)
library(DESeq2) # This is to normalize the data and to compare count differences according to a certain statistic
library(EnhancedVolcano) # This is to plot nice volcano plots. THese plot fold change of genes vs. the significance of the differnce
library(limma)
library(MSigDB)
library(GSVA)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(dnet)
library(igraph)
library(progeny)
library(dorothea)
library(readr)
library(CARNIVAL)
library(patchwork)

# HELPER FUNCTIONS --------------------------------------------------------

## 'quantile_breaks' optimizes color range for heatmap plotting.
#'@param xs, vector of values to adjust palette breaks to
#'@param n, number of steps in color palette
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


# CARNIVAL FUNCTIONS ------------------------------------------------------

## 'networkCARNIVAL' constructs a network object out of CARNIVALs output
#'@param input.obj, unformatted CARNIVAL output
#'@param weightCut, threshold to dismiss edges with low weight
#'@param clusterSize, only allow clusters of given size in resulting graph
networkCARNIVAL <- function(input.obj, weightCut=0, clusterSize=2) {
  
  cRamp <- colorRampPalette(c("#1F77B4", "whitesmoke","#D62728"))(10)  
  
  C.net <- input.obj$weightedSIF %>%
    # transform to data.frame
    as.data.frame() %>%
    # make Node1 and Node2 the first columns of the data.frame
    dplyr::relocate(Node1) %>%
    dplyr::relocate(Node2, .after = Node1) %>%
    # transform Weight attribute to numeric and apply cutoff to remove unwanted edges
    dplyr::mutate( Weight = as.numeric(Weight) ) %>% 
    dplyr::filter( Weight >= weightCut ) %>%
    # edge colors - light blue: activation & light red: inhibition
    dplyr::mutate( edge.color = ifelse(Sign == 1, "#AEC7E8", ifelse(Sign == -1, "#FF9896","#C7C7C7"))) %>%
    dplyr::mutate( edge.label = ifelse(Sign == 1, "activation", ifelse(Sign == -1, "inhibition","connection"))) %>%
    dplyr::mutate( edge.type = ifelse(Weight >= 50, "solid", "dotted")) %>%
    #dplyr::mutate( edge.width = as.numeric(cut(Weight, breaks=c(0,25,50,75,100))))
    dplyr::mutate( edge.width = log(Weight, 5))
  
  
  nodes <- input.obj$nodesAttributes %>% 
    as.data.frame() %>% 
    dplyr::mutate(act_sign = sign(as.numeric(AvgAct))) %>%
    dplyr::mutate(AvgAct = as.double(AvgAct)) %>%
    #dplyr::filter(as.numeric(ZeroAct) != 100) %>%
    dplyr::mutate(NodeType = ifelse(nchar(NodeType)==0, "Protein",ifelse(NodeType == "T","TF","Perturbed"))) %>%
    dplyr::mutate( node.title = paste0(Node,"\n","Up activity: ",UpAct,"\n", "Down activity: ",DownAct,sep="")) %>%
    dplyr::mutate( node.color = cRamp[cut(AvgAct,breaks = 10)]) %>%
    dplyr::mutate( node.shape = ifelse(NodeType == "TF", "circle", ifelse(NodeType == "Perturbed", "square", "sphere")))
  
  ## this is tricky .. so its separate 
  # find all the nodes that ARE TARGETED by perturbed nodes (S)
  target.by.perturbed <- unique(C.net$Node2[which(C.net$Node1 %in% unique(nodes$Node[which(nodes$NodeType=="Perturbed")]))])
  # find all the nodes that TARGET TFs (T)
  target.tfs <- unique(C.net$Node1[which(C.net$Node2 %in% unique(nodes$Node[which(nodes$NodeType=="TF")]))])
  # And then apply
  nodes$node.level <- 3
  nodes$node.level[which(nodes$Node %in% target.by.perturbed)] = 2
  nodes$node.level[which(nodes$Node %in% target.tfs)] = 4
  nodes$node.level[which(nodes$NodeType=="Perturbed")] = 1
  nodes$node.level[which(nodes$NodeType=="TF")] = 5
  nodes$node.level[which(nodes$NodeType=="Pathway")] = 6 # i don't get this .. fix later maybe
  
  # this produces a network of all the nodes in the CARNIVAL output
  # so afterwards we will need to filter in some way
  C.net <- C.net %>%
    igraph::graph_from_data_frame(directed = TRUE, vertices = nodes)  
  
  
  deg <- degree(C.net, mode="all")
  V(C.net)$size <- as.numeric(factor(deg))
  
  # get component statistics and create subgraph that only contains connected nodes
  tmp <- igraph::components(C.net)
  tmp_clusters <- which(tmp$csize >= clusterSize)
  C.net <- igraph::induced_subgraph(C.net, names(tmp$membership[which(tmp$membership %in% tmp_clusters)]))
  
  return(C.net)
}


## 'plotCnet' visualizes the given graph using 'visNet' 
#'@param graph.obj, needs to be output of 'networkCARNIVAL' to ensure all parameters are set.
#'@param layoutType, 1:Fruchtermann-Reingold, 2:GraphOpt, 3:Circle, 4:Sphere, 5:Hierarchy
#'@param f.title, title of output generated in workdir, can also be full or relative path.
#'@param draw.legend, default TRUE
#'@param draw.description, default TRUE
#'@ToDo Hierarchy layout needs some work
plotCnet <- function(graph.obj, layoutType=1, f.title="CTF",draw.legend=TRUE, draw.description=TRUE) {
  # switch between different preset layouts
  l <- switch(layoutType,
              layout_with_fr(graph.obj, niter = 2000),
              layout_with_graphopt(graph.obj),
              layout_in_circle(graph.obj),
              layout_on_sphere(graph.obj),
              layout_with_sugiyama(graph.obj, layers = V(graph.obj)$node.level, hgap = 2,vgap = 2))
  
  ## ToDo: this needs work..
  if(layoutType == 5) {
    #graph.obj <- l$extd_graph
    l <- l$layout
  }
  
  if(draw.description) {
    s.title <- "Node color indicates average activation I=[-100,100] from blue to red. \n Edge width indicates strength of interaction I=[0,100]. \n Interaction strength lower 50 produces dotted lines."
  } else { 
    s.title = NULL}
  
  ## Proceed with plotation
  pdf(file = paste(f.title,".pdf",sep=""), height = 16, width = 16, useDingbats = F, onefile=TRUE)
  # plot graph and transform to ggplot-type to add legends
  plot.igraph(graph.obj, newpage=F,
              # vertex parameters
              vertex.shape=V(graph.obj)$node.shape, 
              vertex.size = scales::rescale(V(graph.obj)$size, to = c(2,7)),
              vertex.color = V(graph.obj)$node.color, vertex.frame.color = "black",
              vertex.label.cex=0.75, vertex.label.dist = 0.5, vertex.label.color = "black",
              vertex.label.family = "Times", vertex.label.font = 1,
              #vertex.label = V(graph.obj)$node.level,
              # edge parameters
              edge.lty = E(graph.obj)$edge.type,
              edge.color=E(graph.obj)$edge.color,
              edge.arrow.size=E(graph.obj)$edge.width / 5, #.99,
              edge.curved=.2,
              edge.width=E(graph.obj)$edge.width, #scales::rescale(E(graph.obj)$Weight, to = c(1,4)), #E(graph.obj)$edge.width, center=F,
              # Nope
              #edge.label = E(graph.obj)$edge.label,
              #edge.label.cex=0.25, edge.label.dist = 0.5, edge.label.color = "black",
              # grouping and layout
              main = "CARNIVAL TF-network", sub = s.title,
              mark.shape=1, mark.expand=10, #mark.groups = list(V(graph.obj)$name[V(graph.obj)$node.level == 1]),
              layout = l)
  if(draw.legend) {
    legend("topleft", 
           #title = "Node types and status",
           legend=c("TF","Perturbed","Protein","activation","inactivation"), 
           box.col="transparent", pch=c(1,0,16,NA,NA), 
           col = c("black", "black", "black","black","black"),
           cex=1, bty="n",
           #horiz=TRUE,
           #xpd=TRUE, inset=c(0, -.15), cex=.8,
           ncol=1)
    par(font = 5) #change font to get arrows
    legend("topleft", 
           legend = c(NA,NA,NA,NA,NA), pch = c(NA,NA,NA, 174, 174),
           #lwd = 1, 
           #lty = c(1,1,1,NA,NA),
           col = c(NA, NA, NA,"#AEC7E8","#FF9896"),
           cex=1, bty="n",
           ncol=1)
    par(font = 1) #back to default
  }
  # and stop plotation
  dev.off()
  
}



# GSE FUNCTIONS -----------------------------------------------------------

# Code for hypergeometric test
GSE_analysis <- function(geneList,Annotation_DB){
  
  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]
  
  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")
  
  DB_genecontent = length(unique(unlist(Annotation_DB)))
  
  GenesDB = DB_genecontent 
  SelectedGenes = length(geneList)
  
  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))
    
    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    #ResultsDF[gset,"p_value"] = HypergeometricTest(overRepres = TRUE,N = N,K = K,n = n,k = k)
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, 
                                       k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }
  
  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]
  
  return(ResultsDF)
  
}

## 'plotGSE' plot significantly enriched pathways
#'@param input.obj, unformatted CARNIVAL output
#'@param p.select, either "p_value" or "corr_p_value"
#'@param p.threshold, determines which pathways to include
plotGSE <- function(input.obj, p.title = "", p.select = "corr_p_value", p.threshold = 0.05) {
  
  res.plot <- input.obj %>% 
    tibble::rownames_to_column(var = "set") %>% 
    tibble %>% 
    dplyr::mutate( corr_p_value = as.numeric(corr_p_value) ) %>% 
    plyr::mutate( p_value = as.numeric(p_value) ) %>% 
    #dplyr::arrange( dplyr::desc(abs(corr_p_value)) ) %>% 
    dplyr::mutate( set = gsub(pattern = "REACTOME_", "", set)) %>% 
    dplyr::mutate( set = gsub(pattern = "_", " ", set)) %>% 
    dplyr::mutate( set = stringr::str_trunc( set, 40 ) ) %>% 
    dplyr::filter( get(p.select) < p.threshold ) %>% 
    dplyr::mutate( y.value = -log10(get(p.select)) ) %>% 
    #head(25) %>% 
    #dplyr::mutate( set = substr(set, 1, 40)) %>% 
    ggpubr::ggdotchart(., x = "set", y = "y.value", xlab = "", ylab = paste("-log10(",p.select,")", sep=""),
                       dot.size = 4, title = eval(p.title),
                       add.params = list(color = "lightgray", size = 2),
                       add = "segments",
                       rotate = TRUE
    ) + theme(axis.text=element_text(size=8), 
              legend.position = "none",
              axis.text.x = element_text(size = 12),
              axis.title=element_text(size=14,face="bold")) 
  
  return(res.plot)
}


# GSVA FUNCTIONS ----------------------------------------------------------

dat=res.REACTOME
sourceDB=REACTOME.DB
setOverlap=0.995
fieldPvalue="P.Value"
cutPvalue=0.05 
clusterSize=2
p.title="REACTOME.ssGSVA.Thyroid"
cutString="REACTOME_"


## 'plotGSVAGraph' plot network of significantly enriched pathways
#'@param dat, output of limma 'topTable' function 
#'@param soureDB, MSigDb that enrichment was performed on
#'@param setOverlap, overlap in gene-sets to draw an edge I=[0,1]
#'@param fieldPvalue, colname of p-value in dat
#'@param cutPvalue, significance threshold
#'@param clusterSize, minimum allowed nodes in clusters
#'@param cutString, relates to DB-prefix, e.g. "HALLMARK_"
#'@param p.title, title for output pdf
plotGSVAGraph <- function(dat, sourceDB, fieldPvalue="P.Value",cutPvalue=0.05, clusterSize=2, p.title="none", cutString, setOverlap=0.5) {
  require(igraph, RedeR)

  # We also require a colorramp to visualize PWs regulatory state
  cRamp <- colorRampPalette(c("#1F77B4", "whitesmoke","#D62728"))(10)  
  
  # prep distance matrix based on gene set overlap for significant sets
  mygroups <- sourceDB[rownames(dat[which(dat[[eval(fieldPvalue)]] <= cutPvalue),])]
  mat <- sapply(mygroups,
                function(y)
                  sapply(mygroups,
                         function(x) length(as.vector(
                           intersect(as.vector(unlist(y)),
                                     as.vector(unlist(x)))))/min(length(as.vector(unlist(y))),
                                                                 length(as.vector(unlist(x))))))
  
  # Draw an edge, if they have min. overlap
  mat[mat < setOverlap] <- 0
  # Create a graph
  groups.net <- igraph::graph.adjacency(mat, weighted=T,mode="undirected",diag=F)
   
  # Now add information to the nodes
  nodes <- dat[which(dat[[eval(fieldPvalue)]] <= cutPvalue),] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "PWID") %>%
    dplyr::mutate( node.color = cRamp[cut(logFoldChange,breaks = 10)]) %>%
    dplyr::mutate( node.shape = ifelse(get(fieldPvalue) <= cutPvalue, "sphere", "square")) %>%
    dplyr::mutate( node.size = sapply( mygroups, length))
  
  # originally creating subgraph of connected components, but here we use it to add node data
  gt3 <- RedeR::subg(g=groups.net, dat=nodes, refcol=1, maincomp=FALSE, connected=FALSE)
  
  # include edge labels - denoting percentage of overlaps
  E(gt3)$edge.label <- paste(round(get.edge.attribute(gt3)$weight * 100, digits=2), "%","\n","overlap",sep="")
  
  
  # rename node names, i.e. if necessary cut of DB prefix
  V(gt3)$name <- gsub(cutString, "", V(gt3)$name)
  V(gt3)$name <- gsub("_", "\n", V(gt3)$name)
  
  # remove singletons - IF they are not significant nodes
  #Pcut <- V(gt3)$name[get.vertex.attribute(gt3, name=eval(fieldPvalue)) <= cutPvalue]
  tmp <- igraph::components(gt3)
  tmp_clusters <- which(tmp$csize >= clusterSize)
  node.keep <- names(tmp$membership[which(tmp$membership %in% tmp_clusters)])
  ## error message
  try(if(length(node.keep) == 0) stop("Cluster threshold too stringent - no nodes left."))
  
  gt3 <- igraph::induced_subgraph(gt3, node.keep)
  
  # calculate connected components for nice polygone drawing
  gr <- cluster_edge_betweenness(gt3)
  gr$com.color <- colorRampPalette(c("#AEC7E8", "whitesmoke","#9EDAE5"))(length(gr))[cut(sizes(gr),breaks = length(gr))]
  
  l <- layout_with_fr(gt3, niter = 2000)
  
  pdf(file = paste(p.title,".pdf",sep=""), height = 16, width = 16, useDingbats = F, onefile=TRUE)
  dnet::visNet(gt3, newpage=F,
               vertex.shape=get.vertex.attribute(gt3, "node.shape", V(gt3)), 
               vertex.label.cex=0.55,  vertex.size = log2(V(gt3)$node.size+1),
               vertex.color = get.vertex.attribute(gt3, "node.color", V(gt3)),
               vertex.frame.color = "black",
               vertex.label.family="serif", 
               vertex.label.dist = 1.1,
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
               edge.label=get.edge.attribute(gt3, name = "edge.label"),
               edge.width= scales::rescale(get.edge.attribute(gt3)$weight, to=c(1,4)), center=F,
               glayout = l)
  dev.off()
  
}


