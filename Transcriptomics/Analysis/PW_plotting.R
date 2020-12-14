# SET-UP ------------------------------------------------------------------
# rm(list = ls())
# Set-up directory for analysis and load required packages and functions.
source("scripts/_RNASeq_functions.R")

# STRINGdb package is required.
# BiocManager::install("STRINGdb")

set.seed(16341)

# Prepare a list to collect our results.
myres_list <- list()

# PREPARATION -------------------------------------------------------------

# Set-up databases for pathway mapping
HALLMARK.DB <- MSigDB::MSigDB[["HALLMARK"]]

Hs.c2 <- MSigDB[["C2_CURATED"]]
KEGG.DB <- Hs.c2[grep("KEGG_", names(Hs.c2))]
REACTOME.DB <- Hs.c2[grep("REACTOME_", names(Hs.c2))]
# Hs.c3 <- MSigDB[["C3_MOTIF"]]
# Hs.c5 <- MSigDB[["C5_GENE_ONTOLOGY"]]

# load the results of the expression analysis
Normalised_counts_matrix <- readRDS(file="data/RDS/Thyroid_normalized_counts.rds", refhook = NULL)
cova <- readRDS(file="data/RDS/Thyroid_covariates.rds", refhook = NULL)
result <- readRDS(file="data/RDS/Thyroid_full_result.rds", refhook = NULL)


# ssGSEA KEGG -------------------------------------------------------------
# We perform the actual enrichment
ssgsea.KEGG <- gsva(
  as.matrix(Normalised_counts_matrix), KEGG.DB, min.sz=10, max.sz=500,
  method="ssgsea", ssgsea.norm=TRUE,
  verbose=TRUE)

# We rename our samples - if we want to.
colnames(ssgsea.KEGG) <- c("healthy_1","healthy_2","healthy_3","early_1","early_2","early_3","late_1","late_2","late_3")
# Now we create the design matrix that allows us to perform comparisons.
design <- model.matrix(~ 0+factor(c(rep(1,3),rep(2,3),rep(3,3))))
colnames(design) <- c("healthy","early","late")             
# Here we specify which comparison we want to make
# In this case we take 'early' as reference level, because that's the direction in which
# we calculated the differential expression. To that end we use the column names of out design
# and say that it is 'late-stage' minus 'early-stage'
contrast.matrix <- makeContrasts(early-healthy,
                                 late-healthy,
                                 late-early, levels=design)
# fit the linear model 
fit  <- limma::lmFit(as.matrix(ssgsea.KEGG), design)
# and get results for our chosen contrast
fit2 <- limma::contrasts.fit(fit, contrast.matrix)
# calculate Empirical Bayes Statistics for Differential Expression
fit2 <- limma::eBayes(fit2)
# extract the top ranked pathways and perform multiple testing correction
res.KEGG <- limma::topTableF(fit2, adjust="BH",number=Inf,sort.by="F", p.value=1)
# adjust colnames for later use
colnames(res.KEGG) <- c("EvH","LvH","LvE","AveExpr","F","P.Value","adj.P.Val")

# take a look at the activity scores for significant pathways
pheatmap::pheatmap(ssgsea.KEGG[rownames(res.KEGG[which(res.KEGG$adj.P.Val <= 0.001),]),],fontsize=8,cellwidth=10,cellheight=10)
# Summary of UP/DOWN regulated pathways
res_HB <- limma::decideTests(fit2, p.value=0.01)
summary(res_HB)



# VISUALIZATION -----------------------------------------------------------

# This will plot a network of significant pathways connected by edges according to the
# overlap in gene sets. Pathway regulation, i.e. are they UP/DOWN regulated in comparison 
# will be shown on a colorRamp blue-white-red. You can select which p-value and cutoff to use.

# 'exp.change' determines the column in 'dat' that contains the lFC values. This plot shows the
# differentially regulated pathways according to the comparison "Early vs Healthy".
plotGSVAGraph(dat=res.KEGG, KEGG.DB, setOverlap=0.35, fieldPvalue="P.Value",cutPvalue=0.05, exp.change = "EvH",
              clusterSize=1, p.title="KEGG.ssGSVA.Thyroid.EvH", cutString="KEGG_" )

# This plot shows the Differentially regulated pathways according to the comparison "Late vs Healthy".
plotGSVAGraph(dat=res.KEGG, KEGG.DB, setOverlap=0.35, fieldPvalue="P.Value",cutPvalue=0.05, exp.change = "LvH",
              clusterSize=1, p.title="KEGG.ssGSVA.Thyroid.LvH", cutString="KEGG_" )

# This plot shows the Differentially regulated pathways according to the comparison "Late vs Healthy".
plotGSVAGraph(dat=res.KEGG, KEGG.DB, setOverlap=0.35, fieldPvalue="P.Value",cutPvalue=0.05, exp.change = "LvE",
              clusterSize=1, p.title="KEGG.ssGSVA.Thyroid.LvE", cutString="KEGG_" )


## Now we will get more specific information about pathways. To that end, we will map 
# all the genes in a pathway to a Protein-Protein-Interaction network (PPI). If we add
# the results of a Differential Expression Analysis, the genes will be colored according
# to the lFC and shaped according to their significance (circle/sphere).

# So first lets have a look at a pathway without DEA information.
# Here we will use the plotting function of stringDB. This will produce a nice plot,
# but can at times be convoluted
plotPWGenes(KEGG.DB$KEGG_GALACTOSE_METABOLISM,        # take a look at the galactose metabolism
            pw.name="Galactose_fancy_noDEA",          # name the plot
            fancy=T)                                  # use stringDB plotting

# Repeat BUT with DEA information - note how significant genes have a red/green aura
# according to their expression.
plotPWGenes(KEGG.DB$KEGG_GALACTOSE_METABOLISM,        
            result$EvH,                               # include DEA results
            pw.name="Galactose_fancy_withDEA",          
            exp.change="log2FoldChange",
            p.select="pvalue", p.cutoff=0.05,
            fancy=T)                                  

## NOW follows the additional plotting functionality. STRINGdb plots can be unwieldy 
# at times, so we can determine the percentage of edges to keep by the PPI interaction 
# strength.

# Here we will use my plotting function. Connected components of the graph 
# will be marked by colored polygons.
plotPWGenes(KEGG.DB$KEGG_GALACTOSE_METABOLISM,        # take a look at the galactose metabolism
            pw.name="Galactose_fancy_noDEA_plain",          # name the plot
            fancy=F)                                  # use stringDB plotting

# Repeat BUT with DEA information - note how genes are colored according to their
# lFC. Significant genes will be displayed as spheres
plotPWGenes(KEGG.DB$KEGG_GALACTOSE_METABOLISM,        
            result$EvH,                               # include DEA results
            pw.name="Galactose_fancy_withDEA_plain",          
            exp.change="log2FoldChange",
            p.select="pvalue", p.cutoff=0.05,
            fancy=F) 

## FINALLY we can 'automatically' produce plots of significant pathways with the 
# 'plotSigPW' wrapper function.
# To that end we select the significant PWs from the results and give them to 
# the plotting function.

res.KEGG <- res.KEGG[order(res.KEGG$adj.P.Val),] 
sig.pathways <- res.KEGG[1:5,] # select top 5 pathways

# And then execute the wrapper function.
plotSigPW(sig.pathways, KEGG.DB, result$EvH, title.prefix="EvH", fancy=T,
          edge.cut=0.9, exp.change="log2FoldChange",
          p.select="pvalue", p.cutoff=0.05)

# Now we repeat the process for the normal plotting method. Note how the color-palette
# for the connected components changes according to their lFC (blue/red). The 'edge.cut'
# parameter determines the percentage of edges to remove based on interaction strength.
plotSigPW(sig.pathways, KEGG.DB, result$EvH, title.prefix="EvH", fancy=F,
          edge.cut=0.9, exp.change="log2FoldChange",
          p.select="pvalue", p.cutoff=0.05)

