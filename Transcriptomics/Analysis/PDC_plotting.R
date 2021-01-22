# SET-UP ------------------------------------------------------------------
# rm(list = ls())
# Set-up directory for analysis and load required packages and functions.
source("scripts/_RNASeq_functions.R")

set.seed(16341)



# CARNIVAL & GENES OF INTEREST --------------------------------------------

# Load all results
cres.LvE <- readRDS(file="data/RDS/CARNIVAL_results_LvE.rds", refhook = NULL)
cres.LvH <- readRDS(file="data/RDS/CARNIVAL_results_LvH.rds", refhook = NULL)
cres.EvH <- readRDS(file="data/RDS/CARNIVAL_results_EvH.rds", refhook = NULL)

# Create vector with genes to investigate
genes.of.interest <- c("MAPK6","MAPKBP1","MAPK10","MAPK8IP3","MAPK3","MAPKAPK5",
                       "MAPKAPK2","MAPK8IP1","MAPK4","MAPK8IP2","MAPK9","MAPK13",
                       "MAPK14","MAPK11","MAPK1IP1L","MAPKAPK3","MAPKAP1","MAPK15",
                       "MAPK1","MAPKAPK5-AS1","MAPK8","MAPK7","MAPK12")


# Construct networks from CARNIVAL results
# parameters work as follows:
# input.obj: a Carnival result list
# weightCut: weighted edges I=[0,100] are removed if lower than threshold
# clusterSize: minimum size of connected clusters to be kept, IF set to 0 while providing GoI vector - this will return GoI and neighborhood only
# GoI: vector with gene denominators, these genes and their neighbors will be kept irrespective of their cluster size
cnet.LvE <- networkCARNIVAL(cres.LvE,                # which results to work with
                            weightCut = 0,           # keep all edges
                            clusterSize = 0,         # in combination with GoI this will return those and their neighborhood only
                            GoI = genes.of.interest) # in this case this are the MAPK-genes
# And plot the network. The plot is exported to pdf Format.
plotCnet(cnet.LvE, layoutType=1, f.title="data/output/CARNIVAL_GoI_MAPK_LvE",draw.legend=F, draw.description=F)

## In case we want to extract all clusters with at least 3 nodes in addition 
# to our genes of interest and remove edges with weak evidence.
cnet.LvE.C3 <- networkCARNIVAL(cres.LvE,                  # which results to work with
                              weightCut = 75,             # drop edges with weight less than 75
                              clusterSize = 3,            # remove clusters smaller than 3
                              GoI = genes.of.interest)    # BUT keep genes of interest and their neighborhood
# And plot the network. The plot is exported to pdf Format.
plotCnet(cnet.LvE.C3, layoutType=1, f.title="data/output/CARNIVAL_GoI_MAPK_LvE_C3",draw.legend=F, draw.description=F)


## And repeat for LvH results
cnet.LvH <- networkCARNIVAL(cres.LvH,                # which results to work with
                            weightCut = 0,           # keep all edges
                            clusterSize = 0,         # in combination with GoI this will return those and their neighborhood only
                            GoI = genes.of.interest) # in this case this are the MAPK-genes
# And plot the network. The plot is exported to pdf Format.
plotCnet(cnet.LvH, layoutType=1, f.title="data/output/CARNIVAL_GoI_MAPK_LvH",draw.legend=F, draw.description=F)

## In case we want to extract all clusters with at least 3 nodes in addition 
# to our genes of interest and remove edges with weak evidence.
cnet.LvH.C3 <- networkCARNIVAL(cres.LvH,                  # which results to work with
                               weightCut = 75,             # drop edges with weight less than 75
                               clusterSize = 3,            # remove clusters smaller than 3
                               GoI = genes.of.interest)    # BUT keep genes of interest and their neighborhood
# And plot the network. The plot is exported to pdf Format.
plotCnet(cnet.LvH.C3, layoutType=1, f.title="data/output/CARNIVAL_GoI_MAPK_LvH_C3",draw.legend=F, draw.description=F)


## And repeat for EvH results
cnet.EvH <- networkCARNIVAL(cres.EvH,                # which results to work with
                            weightCut = 0,           # keep all edges
                            clusterSize = 0,         # in combination with GoI this will return those and their neighborhood only
                            GoI = genes.of.interest) # in this case this are the MAPK-genes
# And plot the network. The plot is exported to pdf Format.
plotCnet(cnet.EvH, layoutType=1, f.title="data/output/CARNIVAL_GoI_MAPK_EvH",draw.legend=F, draw.description=F)

## In case we want to extract all clusters with at least 3 nodes in addition 
# to our genes of interest and remove edges with weak evidence.
cnet.EvH.C3 <- networkCARNIVAL(cres.EvH,                  # which results to work with
                               weightCut = 75,             # drop edges with weight less than 75
                               clusterSize = 3,            # remove clusters smaller than 3
                               GoI = genes.of.interest)    # BUT keep genes of interest and their neighborhood
# And plot the network. The plot is exported to pdf Format.
plotCnet(cnet.EvH.C3, layoutType=1, f.title="data/output/CARNIVAL_GoI_MAPK_EvH_C3",draw.legend=F, draw.description=F)



















