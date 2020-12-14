# SET-UP ------------------------------------------------------------------
# rm(list = ls())
# Set-up directory for analysis and load required packages and functions.
source("scripts/_RNASeq_functions.R")

set.seed(16341)

# DIFFERENTIAL EXPRESSION ANALYSIS ----------------------------------------
# Load in the Count data into the variable cts
# For this you must start your R session from the same directory where your 
# File resised in, or give a path to the file or set the working directory 
# to the directory of the file with the command setwd("Path/to/your/file")
cts <- as.matrix(read.csv("data/RAW/ASMA_Thyroid NGS Tumours Only.txt",sep="\t",row.names="Gene",header=T))

# Now we have to define which columns belong to which tumor type
# First create a data frame with the number of columns in cts as number of rows
# and 1 column for - guessed right - the solid and cystic cancer type
# Take a look what rep("Early",3) does
# Factor is special in R, denoting that this is a condition, and not just a character
samples <- data.frame(row.names = colnames(cts), sID = colnames(cts), condition = as.factor(c(rep("Early",3),rep("ThyPap",3))))
# Check true?
all(rownames(samples) == colnames(cts))

# This does all the magic of initialization
# First reads in the count data from cts
# Gets the information on the smaples from "samples"
# And defines an "experimental design" along the different conditions
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = samples,
                              design = ~ condition)

# Now have a look what you got                              
dds

# Next we delete all Genes that have too little counts across all samples
# We get the count matrix back and calculate the sums of every row
# if this is less than 10 (or whatever you set it), we will discard, else will keep it
keep <- rowSums(counts(dds)) >= 10

# Just keep the rows from keep ;-)
dds <- dds[keep,]

# Look at dds and the dim variable? Did the number of rows change?
dds

# Now we make the statitical test - spits out some steps it performs, we can explain later 
dds <- DESeq(dds)

# Now have a look what we caluclated and which variables were defined
resultsNames(dds)

# And extract the results, but do a shrinkage of the fold change first. will explain later
# We extract exactly the result name
result <- lfcShrink(dds, coef="condition_ThyPap_vs_Early", type="apeglm")

# Order the result
result <- result[order(result$pvalue),]

# Look at the most differentially regulated genes
# The columns denote the mean expression, the fold change, the standard error of the fold change, 
# the p-value and the corrected p-value for multiple testing
head(result)

### A COMMENT ON THE DIRECTION OF DIFFERENTIAL EXPRESSION ###
# In the previous code we see that the condition is ThyPap versus Early, 
# i.e., we calculate the difference in expression with ThyPap/Early
# That means:
# For a HIGHER expression of a particular gene in ThyPap compared to Early we obtain a POSITIVE FoldChange
# For a LOWER expression of a particular gene in ThyPap compared to Early we obtain a NEGATIVE FoldChange
# With this knowledge we take a look at the ordered results.
head(result)
#           baseMean  log2FoldChange  lfcSE     pvalue      padj
#           <numeric> <numeric>       <numeric> <numeric>   <numeric>
#   STAR    887.4095  -2.576015       0.532963  5.68891e-08 0.000183524
# Since 'STAR' has a negative FoldChange, we know that it is expressed to a lesser degree in the
# ThyPap samples --> have a look at the normalized counts
counts(dds, normalized=T)[which(rownames(counts(dds)) %in% c("STAR")),]
# We can see that the first three samples have much higher expression than the latter.
### A COMMENT ON THE DIRECTION OF DIFFERENTIAL EXPRESSION ###


# Finally we save our data for easy reloading
saveRDS(dds, file="data/RDS/Thyroid_dds.rds", compress = TRUE)
saveRDS(result, file="data/RDS/Thyroid_result.rds", compress = TRUE)

# And export results to excel table
openxlsx::write.xlsx(result,file = paste("data/output/Thyroid.differential.expression.LATEvsEARLY.xlsx",sep=""),
                     rowNames=TRUE,colNames=TRUE, colWidths = c(NA, "auto", "auto"))

# VOLCANO PLOT ------------------------------------------------------------
# Now we can start visualizing our results.
# If you start from here and don't want to repeat your analyses, you can load the results.
result <- readRDS(file="data/RDS/Thyroid_result.rds", refhook = NULL)

# pdf() and dev.off() will export your plot to a pdf-file in your working directory
pdf(file=paste("data/output/volcano.Thyroid_Late_vs_Early.pdf",sep=""), width=8, height=8, onefile=T)
EnhancedVolcano::EnhancedVolcano(result,                   # our data
                                 lab = rownames(result),   # how to name the dots later on
                                 x = 'log2FoldChange',     # where do we find the log fold column
                                 y = 'pvalue',             # where the p-value
                                 xlim = c(-4, 4),          # plot from where to where on the x axis
                                 ylim = c(0, 8),           # and on the y axis
                                 legendPosition = 'bottom', # where to put the legend
                                 pCutoff = 0.001,          # what do we consider significant (p-value)
                                 FCcutoff = 1,             # what do we consider significant (fold change)     
                                 pointSize = 1.0,          # Some other parameters
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.2,
                                 colConnectors = 'grey30',
                                 title = 'Papillary thyroid cancer vs. early cancer', # plot title
                                 subtitle = 'pvalue  < 0.001',
                                 xlab = "logFC") # how to name the x axis
dev.off()



# HEATMAP -----------------------------------------------------------------
# In a next step we create a heatmap for the gene expression.
# A heatmap is a grid in which the you can display the 
# genes/samples row or column-wise. 
# The gene expression is then encoded in color.
# You can cluster rows and columns so that similar samples/gene 

# Load the data 
dds <- readRDS(file="data/RDS/Thyroid_dds.rds", refhook = NULL)
result <- readRDS(file="data/RDS/Thyroid_result.rds", refhook = NULL)

# Next we prepare the data we want to show in the heatmap. 
# We must extract the expression data of the genes from the dds object
# In order to compare the results, we must first "normalize" the data
# Here with a Variance stabilization method
vsd <- DESeq2::vst(dds, blind=FALSE,fitType='local')

# Then we can extract the data with the "assay" command
mygeneexpression <- SummarizedExperiment::assay(vsd)

# It is no use to show all 12294 genes in a heatmap. 
# Therefore we restrict ourselves to the most differentially regulated genes
# Let us select genes having an adjusted p-value < 0.1
# The diff. genes are stored in the results data frame
idx <- rownames(result[which(result$padj < 0.1),])  # the variable idx now contains the names of significant genes

# How many genes did we find? 
length(idx)

# reduce list to those significant genes
tf_data.filter <- mygeneexpression[idx,]  # because the rownames of our matrix are the gene names, we cann just select the ones we would like to keep

# transform to log-scale to 
plot.df <- log10(tf_data.filter)

# how much nuance is the colormap supposed to have
col_cnt <- 10
mat_breaks <- quantile_breaks(seq(min(plot.df), max(plot.df), length.out = col_cnt), n = col_cnt + 1)

# produce the plot and put into the variable 'complete'
complete <- pheatmap::pheatmap(plot.df, color = viridis::inferno(length(mat_breaks) - 1),
                               breaks = mat_breaks, scale = "none",
                               cluster_rows = T, clustering_distance_rows = "correlation",
                               cluster_cols = F,
                               fontsize_row = 5, fontsize_col = 12,
                               main = "Papillary thyroid cancer vs. early cancer - log10(normalized)")
# save to pdf-file
pdf(file=paste("data/output/HeatMap.Thyroid.DEA.pdf",sep=""), width=6, height=12, onefile=T)
print(complete)
dev.off()

# PRINCIPAL COMPONENT ANALYSIS --------------------------------------------
# load the data
dds <- readRDS(file="data/RDS/Thyroid_dds.rds", refhook = NULL)
result <- readRDS(file="data/RDS/Thyroid_result.rds", refhook = NULL)

# retrieve the metadata 
samples <- SummarizedExperiment::colData(dds)

vsd <- vst(dds, blind=FALSE,fitType='local')    # normalization
mygeneexpression <- assay(vsd)                  # get counts

# select significant genes
idx <- rownames(result[which(result$padj < 0.1),])  # the variable idx now contains the names of significant genes

# How many genes did we find? 
length(idx)

# reduce list to those significant genes
tf_data.filter <- mygeneexpression[idx,]  # because the rownames of our matrix are the gene names, we cann just select the ones we would like to keep

# sort table in decreasing order
iqr <- apply(mygeneexpression,1,IQR)
tpm_norm <- mygeneexpression[order(iqr,decreasing=T),]
# 1.3 PCA on 5000 most variable genes
PCA <- prcomp(t(tpm_norm), scale = T)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dataGG = data.frame(PC1=PCA$x[,1],PC2=PCA$x[,2],PC3=PCA$x[,3],PC4=PCA$x[,4])    

pca1 <-ggplot2::ggplot() + 
  geom_point(data=dataGG, aes(PC1, PC2, group=samples$condition, colour=samples$condition,
                              shape=samples$condition), size=2.5) + 
  #geom_line(aes(x=dataGG$PC1, y=dataGG$PC2, group=plotdata$response), size=0.25, , alpha=0.5) + 
  labs(title=paste(unique(samples$condition), "pval < 0.05",sep=" "),
       x=paste("PC1 [", percentVar[1], "%]", sep=""), y = paste("PC2 [", percentVar[2], "%]", sep="")) +
  ggrepel::geom_text_repel(data=dataGG, ggplot2::aes(PC1, PC2,label=rownames(samples)))

# save to file - with a different method this time
ggsave(filename="data/output/PCA_Thyroid_DEA_Late_vs_Early_sigGenes.pdf", plot=pca1, width=6, height=6)

# the easy transformation
rld <- vst(dds, blind=FALSE,fitType='local')
plotPCA(rld, ntop=500)


# DOWNSTREAM PREPARATION --------------------------------------------------
## To facilitate more convenient downstream analyses, we will prepare the required data.
# First we will prepare our variables so that we can progress in the investigation.
dds <- readRDS(file="data/RDS/Thyroid_dds.rds", refhook = NULL)

# For downstream analysis we normalize our count matrix.
vsd <- DESeq2::vst(dds, blind=FALSE,fitType='local')   # normalization
Normalised_counts_matrix <- SummarizedExperiment::assay(vsd)   # get counts

## We also extract our covariate information from the desseq object and store them separately for later use.
cova <- SummarizedExperiment::colData(dds)

# Sanity check - we just test if sample IDs match in both datasets
nrow(cova); ncol(Normalised_counts_matrix)
Normalised_counts_matrix <- Normalised_counts_matrix[ , colnames(Normalised_counts_matrix) %in% cova$sID ]
cova            <- cova[ cova$sID %in% colnames(Normalised_counts_matrix),  ]
nrow(cova) == ncol(Normalised_counts_matrix)
rownames( cova ) == colnames( Normalised_counts_matrix )

# And if everything is ok we can save out files
saveRDS(Normalised_counts_matrix, file="data/RDS/Thyroid_normalized_counts.rds", compress = TRUE)
saveRDS(cova, file="data/RDS/Thyroid_covariates.rds", compress = TRUE)
