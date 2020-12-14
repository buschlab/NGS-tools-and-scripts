# SET-UP ------------------------------------------------------------------
# rm(list = ls())
# Set-up directory for analysis and load required packages and functions.
source("scripts/_RNASeq_functions.R")

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

# ssGSEA HALLMARK ---------------------------------------------------------
## IQR based filtering - distribution in this data does not warrant filtering.
# iqr <- apply(Normalised_counts_matrix,1,IQR)
# plot(ecdf(iqr))
# keep <- iqr >= quantile(iqr, 0.25)
# Normalised_counts_matrix <- Normalised_counts_matrix[keep,]

# perform actual enrichment
ssgsea.HALLMARK <- gsva(
  as.matrix(Normalised_counts_matrix), HALLMARK.DB, min.sz=10, max.sz=500,
  method="ssgsea", ssgsea.norm=TRUE,
  verbose=TRUE)

# have a look at the output
head( ssgsea.HALLMARK[order(ssgsea.HALLMARK[,1], decreasing=T),], 10 )

# We rename our samples - if we want to.
colnames(ssgsea.HALLMARK) <- c("early_1","early_2","early_3","late_1","late_2","late_3")
# construct our design to facilitate comparisons
design <- model.matrix(~ 0+cova$condition)
colnames(design) <- c("early","late")             
# here we specify which comparison we want to make
# In this case we take 'early' as reference level, because that's the direction in which
# we calculated the differential expression. To that end we use the column names of out design
# and say that it is 'late-stage' minus 'early-stage'
contrast.matrix <- limma::makeContrasts(late-early, levels=design)
# fit the linear model 
fit  <- limma::lmFit(as.matrix(ssgsea.HALLMARK), design)
# and get results for our chosen contrast
fit2 <- limma::contrasts.fit(fit, contrast.matrix)
# calculate Empirical Bayes Statistics for Differential Expression
fit2 <- limma::eBayes(fit2)
# extract the top ranked pathways and perform multiple testing correction
res.HALLMARK <- limma::topTableF(fit2, adjust="BH",number=Inf,sort.by="F", p.value=1)
# adjust colnames for later use
colnames(res.HALLMARK) <- c("logFoldChange","AveExpr","F","P.Value","adj.P.Val")

# take a look at the activity scores for significant pathways
pheatmap::pheatmap(ssgsea.HALLMARK[rownames(res.HALLMARK[which(res.HALLMARK$adj.P.Val <= 0.1),]),],fontsize=8,cellwidth=10,cellheight=10)
# Summary of UP/DOWN regulated pathways
res_HB <- limma::decideTests(fit2, p.value=0.1)
summary(res_HB)

## And plot the network
plotGSVAGraph(res.HALLMARK, HALLMARK.DB, setOverlap=0.2, fieldPvalue="P.Value",cutPvalue=0.1, 
              clusterSize=1, p.title="data/output/HALLMARK.ssGSVA.Thyroid", cutString="HALLMARK_" )


# ssGSEA KEGG -------------------------------------------------------------
# perform actual enrichment
ssgsea.KEGG <- gsva(
  as.matrix(Normalised_counts_matrix), KEGG.DB, min.sz=10, max.sz=500,
  method="ssgsea", ssgsea.norm=TRUE,
  verbose=TRUE)

# have a look at the output
head( ssgsea.KEGG[order(ssgsea.KEGG[,1], decreasing=T),], 10 )

adjPvalueCutoff <- 0.1
logFCcutoff <- log2(2)

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

## And plot the network
plotGSVAGraph(res.KEGG, KEGG.DB, setOverlap=0.35, fieldPvalue="P.Value",cutPvalue=0.05, 
              clusterSize=1, p.title="KEGG.ssGSVA.Thyroid", cutString="KEGG_" )

res.KEGG <- res.KEGG[1:5,]


plotSigPW(res.KEGG[1:5,], KEGG.DB, result$EvH)


plotPWGenes(KEGG.DB$KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS, 
            result$EvH,
            pw.name="UBIQUITIN_MEDIATED_PROTEOLYSIS",
            fancy=T)

plotPWGenes(KEGG.DB$KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS,
            pw.name="test_2",
            fancy=F)




#

# ssGSEA REACTOME ---------------------------------------------------------
# perform actual enrichment
ssgsea.REACTOME <- gsva(
  as.matrix(Normalised_counts_matrix), REACTOME.DB, min.sz=10, max.sz=500,
  method="ssgsea", ssgsea.norm=TRUE,
  verbose=TRUE)

# have a look at the output
head( ssgsea.REACTOME[order(ssgsea.REACTOME[,1], decreasing=T),], 10 )

# We rename our samples - if we want to.
colnames(ssgsea.REACTOME) <- c("early_1","early_2","early_3","late_1","late_2","late_3")
# construct our design to facilitate comparisons
design <- model.matrix(~ 0+cova$condition)
colnames(design) <- c("early","late")             
# here we specify which comparison we want to make
# In this case we take 'early' as reference level, because that's the direction in which
# we calculated the differential expression. To that end we use the column names of out design
# and say that it is 'late-stage' minus 'early-stage'
contrast.matrix <- limma::makeContrasts(late-early, levels=design)
# fit the linear model 
fit  <- limma::lmFit(as.matrix(ssgsea.REACTOME), design)
# and get results for our chosen contrast
fit2 <- limma::contrasts.fit(fit, contrast.matrix)
# calculate Empirical Bayes Statistics for Differential Expression
fit2 <- limma::eBayes(fit2)
# extract the top ranked pathways and perform multiple testing correction
res.REACTOME <- limma::topTableF(fit2, adjust="BH",number=Inf,sort.by="F", p.value=1)
# adjust colnames for later use
colnames(res.REACTOME) <- c("logFoldChange","AveExpr","F","P.Value","adj.P.Val")

# take a look at the activity scores for significant pathways
pheatmap::pheatmap(ssgsea.REACTOME[rownames(res.REACTOME[which(res.REACTOME$adj.P.Val <= 0.1),]),],fontsize=8,cellwidth=10,cellheight=10)
# Summary of UP/DOWN regulated pathways
res_HB <- limma::decideTests(fit2, p.value=0.1)
summary(res_HB)

## And plot the network
plotGSVAGraph(res.REACTOME, REACTOME.DB, setOverlap=0.995, fieldPvalue="P.Value",cutPvalue=0.05, 
              clusterSize=2, p.title="REACTOME.ssGSVA.Thyroid", cutString="REACTOME_" )















