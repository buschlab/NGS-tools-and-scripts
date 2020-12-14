# SET-UP ------------------------------------------------------------------
# rm(list = ls())
# Set-up directory for analysis and load required packages and functions.
source("scripts/_RNASeq_functions.R")

set.seed(16341)

# Prepare a list to collect our results.
myres_list <- list()

# C6: Onco -> 189 gene sets
onco_gsets <- msigdf::msigdf.human %>% 
  dplyr::filter(category_code == "c6" ) %>% 
  dplyr::select(geneset, symbol) %>% 
  dplyr::group_by(geneset) %>% 
  dplyr::summarize(symbol=list(symbol)) %>% 
  deframe() 


# PROGENy - 1. STEP: Computation ------------------------------------------
# First we load our Expression data set and the results of the differential expression
# analysis. In this case we are working with output from DESeq2, but this script can
# be adjusted to work with most DEA frameworks.
Normalised_counts_matrix <- readRDS(file="data/RDS/airways_normalized_counts.rds", refhook = NULL)
cova <- readRDS(file="data/RDS/airways_covariates.rds", refhook = NULL)
result <- readRDS(file="data/RDS/airways_result.rds", refhook = NULL)

# Pathway activities can be computed sample-wise on normalized counts or for the
# t-statistic of the case/control comparison. To that end we need to calculate t - 
# which is the ratio of fold change to its standard error.
deseq_t <- as.matrix(data.frame(row.names=rownames(result), t=result$log2FoldChange / result$lfcSE))

# We are already prepared to execute PROGENy on our count data with the Human organism data set.
PathwayActivity_counts <- progeny::progeny(Normalised_counts_matrix, scale=TRUE, 
                                           organism="Human", top = 100)
## Progeny somehow messes up the PW-naming - we will fix that now.
colnames(PathwayActivity_counts) <- gsub("-",".",colnames(PathwayActivity_counts))

PathwayActivity_stats <- progeny::progeny(deseq_t, 
                          scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t ()

colnames(PathwayActivity_stats) <- "t"

## As a result we obtain the  matrix'PathwayActivity_counts' that contains sample-wise PW-activity 
# scores. We will now use limma to fit a linear model on the basis of the study design. This facilitates
# two things: 
# 1. The statistical model helps us in the evaluation of each pathways significance.
# 2. The moderated t-statistc can be used as input for the final computation in CARNIVAL.

## To execute the limma steps we transpose our output and change sample names (the latter step is optional)
PathwayActivity_counts <- t(PathwayActivity_counts)
colnames(PathwayActivity_counts) <- c("early_1","early_2","early_3","late_1","late_2","late_3")
# Now we create the design matrix that allows us to perform comparisons.
design <- model.matrix(~ 0+factor(c(rep(1,3),rep(2,3))))
colnames(design) <- c("early","late")             
# Here we specify which comparison we want to make
# In this case we take 'early' as reference level, because that's the direction in which
# we calculated the differential expression. To that end we use the column names of out design
# and say that it is 'late-stage' minus 'early-stage'
contrast.matrix <- makeContrasts(late-early, levels=design)
# Now we fit the linear model.
fit.PPA  <- lmFit(as.matrix(PathwayActivity_counts), design)
# Compute results for our chosen comparison.
fit.PPA.contrast <- contrasts.fit(fit.PPA, contrast.matrix)
# And calculate Empirical Bayes Statistics for this computation.
fit.PPA.contrast <- eBayes(fit.PPA.contrast)
# Extract the top ranked pathways and perform multiple testing correction - this gives us the t-statistic
# We are not using a p.value cutoff (p.value=1) because we want toretain all the information at this point.
res.progeny <- topTable(fit.PPA.contrast, adjust="BH",number=Inf,sort.by="t",p.value=1) 

# We can also use 'decideTests' to the direction of expression of the Pathways.
res.progeny.direction <- decideTests(fit.PPA.contrast, p.value=0.1)
# And summary gives us a qick look at the outcome.
summary(res.progeny.direction)

# The TF activity enrichment results provided by Viper are used as an input in the 
# CARNIVAL method. CARNIVAL tries to infer the most likely upstream signaling events 
# leading to the current TF activity results.
## We save the results for later use, before continuing with the analysis
# For this we put together all computed outputs.
resList.progeny <- list()
resList.progeny[["counts"]] <- PathwayActivity_counts
resList.progeny[["stats"]] <- PathwayActivity_stats
resList.progeny[["results"]] <- res.progeny
saveRDS(resList.progeny, file="data/RDS/PathwayActivity_results_Late_vs_Early.rds", compress = TRUE)

## And put it in our final table
myres_list[["PW counts"]] <- PathwayActivity_counts
myres_list[["PW t-statistic"]] <- PathwayActivity_stats
myres_list[["PW results"]] <- res.progeny

# PROGENy - 2. STEP: Visualization & Investigation ------------------------
# As the input for CARNIVAL goes we are done with PROGENy. Naturally, we can investigate
# these results a little further. A first step is to create a heatmap of calculated Pathway
# activity scores per sample to check it for conspicuous clustering.

## pHEATMAP
# rotate the matrix to facilitate usage in pheatmap
plot.df <- PathwayActivity_counts
# how much nuance is the colormap supposed to have
col_cnt <- 10
#myColor <- colorRampPalette(c("darkblue", "whitesmoke","indianred"))(col_cnt)
mat_breaks <- quantile_breaks(seq(min(plot.df), max(plot.df), length.out = col_cnt), n = col_cnt + 1)

# produce the plot and put into the variable 'complete'
PW.heat <- pheatmap::pheatmap(plot.df, color = viridis::inferno(length(mat_breaks) - 1),
                               breaks = mat_breaks, scale = "none",
                               cluster_rows = T, clustering_distance_rows = "correlation",
                               cluster_cols = F,
                               fontsize_row = 10, fontsize_col = 10,
                               treeheight_col = 0,  border_color = NA,
                               main = "PROGENy (100)")
# save to pdf-file
pdf(file=paste("data/output/HeatMap.Progeny_Late_vs_Early.pdf",sep=""), width=6, height=12, onefile=T)
print(PW.heat)
dev.off()

# Because we already fitted a model and calculated fold-changes and p-values, we can now
# plot those statistics for the pathway activity as well. We could for example use a barplot
# to display the fold-change or t-statistic.

## BARPLOT pathway activity
PathwayActivity_zscore_df <- as.data.frame(res.progeny) %>% 
  dplyr::rename(NES = "t") %>%
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(logFC) %>%
  dplyr::mutate(Pathway = factor(Pathway))

pdf(file=paste("data/output/BarPlot.Progeny.PWactivity.statistic_Late_vs_Early.pdf",sep=""), width=6, height=6, onefile=T)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, logFC), y = logFC)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = cols[1], high = cols[3], 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")
dev.off()


# Dorothea - 1. STEP: Computation -----------------------------------------
# The package Dorothea calculates activity scores for Transcription Factors (TF) based 
# on expression data. To be precise, the package provides a Wrapper for 'Viper' which 
# actually does the work. So we will now compute TF activity scores based on the normalized
# counts. Subsequently employ limma to fit a linear model and perform statistic analyses.
# This provides us with the means to further evaluate the output of Dorothea, as well as
# with input for the CARNIVAL analysis step.

# First we load our Expression data set and the results of the differential expression
# analysis. In this case we are working with output from DESeq2, but this script can
# be adjusted to work with most DEA frameworks.
Normalised_counts_matrix <- readRDS(file="data/RDS/airways_normalized_counts.rds", refhook = NULL)
cova <- readRDS(file="data/RDS/airways_covariates.rds", refhook = NULL)
result <- readRDS(file="data/RDS/airways_result.rds", refhook = NULL)

# Transcription Factor activities can be computed sample-wise on normalized counts or 
# for the t-statistic of the case/control comparison. To that end we need to calculate 
# t - which is the ratio of fold change to its standard error.
deseq_t <- as.matrix(data.frame(row.names=rownames(result), t=result$log2FoldChange / result$lfcSE))

# We estimate the transcription factor (TF) activity using the DoRothEA R package. 
# We select interactions with confidence level A, B and C.
## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

# Execute dorotheas 'run_viper' on t-statistic to calculate overall TF-activity scores
TFActivity_stat <- dorothea::run_viper(deseq_t, regulons,
                                          options =  list(minsize = 15, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

# Execute dorotheas 'run_viper' to calculate sample-wise TF-activity scores
TFActivity_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 15, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

## 'TFActivity_counts' contains sample-wise TF-activity scores and we require a single input value
# for the TF-activites - so now we will use limma to fit a linear model and calculate the t-statistic
# somehow
colnames(TFActivity_counts) <- c("early_1","early_2","early_3","late_1","late_2","late_3")
# construct our design to facilitate comparisons
design <- model.matrix(~ 0+factor(c(rep(1,3),rep(2,3))))
colnames(design) <- c("early","late")             
# here we specify which comparison we want to make
# In this case we take 'early' as reference level, because that's the direction in which
# we calculated the differential expression. To that end we use the column names of out design
# and say that it is 'late-stage' minus 'early-stage'
contrast.matrix <- makeContrasts(late-early, levels=design)
# fit the linear model 
fit.DTA  <- lmFit(as.matrix(TFActivity_counts), design)
# and get results for our chosen contrast
fit.DTA.contrast <- contrasts.fit(fit.DTA, contrast.matrix)
# calculate Empirical Bayes Statistics for Differential Expression
fit.DTA.contrast <- eBayes(fit.DTA.contrast)
# extract the top ranked pathways and perform multiple testing correction - this gives us the t-statistic
res.dorothea <- topTable(fit.DTA.contrast, adjust="BH",number=Inf,sort.by="t",p.value=1) # keep p.cutoff at one because we want to keep everything

# We can also use 'decideTests' to the direction of expression of the Pathways.
res.dorothea.direction <- decideTests(fit.DTA.contrast, p.value=0.1)
# And summary gives us a qick look at the outcome.
summary(res.dorothea.direction)

# The TF activity enrichment results provided by Viper are used as an input in the 
# CARNIVAL method. CARNIVAL tries to infer the most likely upstream signaling events 
# leading to the current TF activity results.
## We save the results for later use, before continuing with the analysis

resList.dorothea <- list()
resList.dorothea[["counts"]] <- TFActivity_counts
resList.dorothea[["stats"]] <- TFActivity_stat
resList.dorothea[["results"]] <- res.dorothea
saveRDS(resList.dorothea, file="data/RDS/TFActivity_CARNIVAL_results_Late_vs_Early.rds", compress = TRUE)

## And put it in our final table
myres_list[["TF counts"]] <- TFActivity_counts
myres_list[["TF t-statistic"]] <- TFActivity_stat
myres_list[["TF results"]] <- res.dorothea

# Dorothea - 2. STEP: Visualization & Investigation -----------------------
# As the input for CARNIVAL goes we are done with Dorothea. Naturally, we can investigate
# these results a little further. A first step is to create a heatmap of calculated TF
# activity scores per sample to check it for conspicuous clustering.

# We use the computed statistics to limit the selection to TFs with significant pvalues.
TF_activity <- rownames(res.dorothea[res.dorothea$P.Value <= 0.05,])
# And create our plot-data-frame.
plot.df <- TFActivity_counts %>%
  as.data.frame() %>% 
  dplyr::filter(rownames(TFActivity_counts) %in% TF_activity)

# how much nuance is the colormap supposed to have
col_cnt <- 10
mat_breaks <- quantile_breaks(seq(min(plot.df), max(plot.df), length.out = col_cnt), n = col_cnt + 1)

# produce the plot and put into the variable 'complete'
TF.heat <- pheatmap::pheatmap(plot.df, color = viridis::inferno(length(mat_breaks) - 1),
                               breaks = mat_breaks, scale = "none",
                               cluster_rows = T, clustering_distance_rows = "correlation",
                               cluster_cols = F,
                               fontsize_row = 5, fontsize_col = 12,
                               main = "airways cancer - TF-activity scores",
                               subtitle = 'significant TFs - p <= 0.05')

# save to pdf-file
pdf(file=paste("data/output/HeatMap.airways.TF.activity.scores_Late_vs_Early.pdf",sep=""), width=6, height=12, onefile=T)
print(TF.heat)
dev.off()


# And we create a bar-plot of fold-change or NES (normalized enrichment scores; t-statistic)
# for significant TFs.
tf_activities_stat <- res.dorothea %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t") %>%
  dplyr::filter(P.Value < 0.05) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

TF.scores <- ggplot(tf_activities_stat,aes(x = reorder(GeneID, logFC), y = logFC)) + 
  geom_bar(aes(fill = logFC), stat = "identity") +
  scale_fill_gradient2(low = cols[1], high = cols[3], 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

pdf(file=paste("data/output/Barplot.airways.TF.NES.significant_Late_vs_Early.pdf",sep=""), width=6, height=12, onefile=T)
print(TF.scores)
dev.off()

# To interpret the results, we can look at the expression of targets of one of the 
# most deregulated TFs, such as MYC. To that end we select all the genes that are
# associated with the MYC transcription factor and display their DEA results in a
# volcano plot. Another idea would be to produce a heatmap of expression values for
# all samples in the set - but to that end some more filtering is required, as the 
# number of genes affected by any one TF can be quite large.
targets_MYC <- regulons$target[regulons$tf == "MYC"]
# subset for chosen targets to create our plotting data.frame
plot.df <- result %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% targets_MYC)

# Visualize the expresssion of associated genes in a volcano plot
vp_dorothea <- EnhancedVolcano::EnhancedVolcano(plot.df,                   # our data
                                                lab = plot.df$GeneID,   # how to name the dots later on
                                                x = 'log2FoldChange',     # where do we find the log fold column
                                                y = 'pvalue',             # where the p-value
                                                xlim = c(-4, 4),          # plot from where to where on the x axis
                                                ylim = c(0, 8),           # and on the y axis
                                                legendPosition = 'bottom', # where to put the legend
                                                pCutoff = 0.05,          # what do we consider significant (p-value)
                                                FCcutoff = .5,             # what do we consider significant (fold change)     
                                                pointSize = 1.0,          # Some other parameters
                                                drawConnectors = TRUE,
                                                widthConnectors = 0.2,
                                                colConnectors = 'grey30',
                                                title = 'Dorothea - TF:', # plot title
                                                subtitle = 'MYC target gene expression',
                                                xlab = "logFC") # how to name the x axis

pdf(file=paste("data/output/Volcano.airways.TF.MYC.target.expression_Late_vs_Early.pdf",sep=""), width=6, height=12, onefile=T)
print(vp_dorothea)
dev.off()

## This concludes all the steps to facilitate the execution of CARNIVAL, but the output
# of PROGENy and Dorothea could be investigated in more detail!

# CARNIVAL - STEP 1: Omnipath scaffold ------------------------------------
## CARNIVAL is an approach that works on networks which represent the interactions
# between genes. In this step we will retrieve one of such networks from Omnipath
# to use it as a basis for the execution of CARNIVAL.
AllInteractions = OmnipathR::import_Omnipath_Interactions(select_organism = 9606) 

## We transform to the format needed by CARNIVAL. We just keep signed and 
## directed interactions 
SignedDirectedInteractions <- filter(AllInteractions, is_directed == 1) %>%
  filter(is_stimulation == 1 | is_inhibition == 1)

InputCarnival <- bind_rows(
  (SignedDirectedInteractions %>%
     filter(is_stimulation == 1 & is_inhibition == 0) %>%
     transmute(source_genesymbol, interaction = 1, target_genesymbol)),   
  (SignedDirectedInteractions %>%
     filter(is_stimulation == 0 & is_inhibition == 1) %>%
     transmute(source_genesymbol, interaction = -1, target_genesymbol))) %>%  
  distinct() 

## We have to be careful with the gene names with a "-". CPLEX gets crazy. 
InputCarnival$source_genesymbol <- 
  gsub("-","_",InputCarnival$source_genesymbol)
InputCarnival$target_genesymbol <- 
  gsub("-","_",InputCarnival$target_genesymbol)

bad_int <- which(duplicated(paste(InputCarnival[,1],InputCarnival[,3])))

if ( length(bad_int) == 0 ) {
  InputCarnival <- InputCarnival
} else {
  InputCarnival = InputCarnival[-bad_int,]
}

# This network only changes if the underlying database has been updated. To facilite a
# faster workflow and have a fallback in case of problems with the internet connection,
# we will save the finished network-object.
saveRDS(InputCarnival, file="data/RDS/CARNIVAL.network.rds", compress = TRUE)


# CARNIVAL - STEP 2: prepare PD-input -------------------------------------
# In addition to the interaction network, we need actual study related input.
# The transcription factor activity computed with dorothea is a required input for CARNIVAl.
# The PW activity scores of PROGENy are an optional input- that we will happily provide here.
# 1. Dorothea output
TFActivity_results <- readRDS(file="data/RDS/TFActivity_CARNIVAL_results_Late_vs_Early.rds", refhook = NULL)
# 2. OPTIONAL Progeny output
PathwayActivity_results <- readRDS(file="data/RDS/PathwayActivity_results_Late_vs_Early.rds", refhook = NULL)

## the data requires some adjustment to work with CARNIVAL, we need vectors with gene names 
# and relevant scores
# starting with the pathway activity scores in the progeny output
PathwayActivity_CARNIVALinput <- PathwayActivity_results$stats %>%
  as.data.frame() %>%
  dplyr::rename(score = "t") %>%
  rownames_to_column(var = "Pathway")
  
# load progeny pathways and associated genes/proteins and respective identifiers
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))
# we only require the gene names 
progenyMembers <- progenyMembers$gene
# most pathways have multiple associated genes - we want a vector that contains all the genes and
# the score of their respective pathway
# So we transform the 'progenyMembers' list into a data.frame where the genes are comma separated
# use separate_rows() to split up the commas and insert new rows
# use left_join() to insert the calculated scores we loaded earlier
tmp <- data.frame(dplyr::left_join(
  tidyr::separate_rows(
    data.frame(row.names = NULL,
               Pathway= names(progenyMembers), 
               genes=sapply(progenyMembers, function(pw) paste(pw,collapse = ",")), 
               stringsAsFactors = F), 
    genes, sep=","), 
  PathwayActivity_CARNIVALinput, by="Pathway"),stringsAsFactors = F)
# Then just extract the scores and add gene name to each element =)
progenyScores <- data.frame(t(sort(tmp[,"score"])))
#rownames(progenyScores) = "NES"
colnames(progenyScores) = tmp$genes


# Now we will do something similar with the TF-scores we calculated in Dorothea
# This is much easier - we just need to select our top 25,50,100,.. genes
## sorting by absolute value
TFActivity_CARNIVALinput <- TFActivity_results$results %>%
  as.data.frame() %>% 
  rownames_to_column(var = "TF") %>%
  dplyr::filter(P.Value < 1.05)

# And select the TopXX scores from the 't' column - also transpose because apparently CARNIVAL has the shittiest code ever..
dorotheaScores <- data.frame(t(sort(TFActivity_CARNIVALinput[,"t"])))
  rownames(dorotheaScores) = "NES"
  colnames(dorotheaScores) = TFActivity_CARNIVALinput$TF



# CARNIVAL - STEP 3: Run the Algorithm ------------------------------------
InputCarnival <- readRDS(file="data/RDS/CARNIVAL.network.rds", refhook = NULL)

# get initial nodes - I'm not quite sure what this is supposed to do, because the code produces
# an empty named vector?!
iniMTX = base::setdiff(InputCarnival$source_genesymbol, InputCarnival$target_genesymbol)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX


# We produced all the input required to run CARNIVAL - depending on the complexity of the
# task at hand this can take quite some time.
carnival_result = runCARNIVAL( inputObj= iniciators,
                               measObj = dorotheaScores,        # dorothea output
                               netObj = InputCarnival,          # Network of interacting genes
                               weightObj = progenyScores,       # (OPTIONAL) progeny output
                               solverPath = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex", 
                               solver = "cplex",
                               timelimit=10000,                   # if no solution can be found, enlarge this value
                               dir_name = paste(WORKDIR,"data/output", sep="/"),              # path to the output directory
                               mipGAP=0,
                               poolrelGAP=0 )
# And as always we will save our results
saveRDS(carnival_result, file="data/RDS/CARNIVAL_results_late_vs_early.rds", compress = TRUE)

# Also save it to file for later inspection
write.table(carnival_result$weightedSIF,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_results_late_vs_early_network.csv", quote = F)
write.table(carnival_result$nodesAttributes,col.names = T,row.names = F,
            sep = ",", file = "CARNIVAL_results_late_vs_early_nodes.csv", quote = F)


# This concludes the mostly computational steps in the Progeny, Dorothea, Carnival (PDC) pipeline.
# Next we will continue with the investigation and visualization of our results.



# CARNIVAL - STEP 4: Visualization & Investigation ------------------------
carnival_result <- readRDS(file="data/RDS/CARNIVAL_results_late_vs_early.rds", refhook = NULL)

## Here we extract the network topology from CARNIVAL's output, i.e., this is a list of 
# node pairs (Node1, Node2) connected by directed edges that have a weight attribute (Weight) as
# well as an indication for the type of relationship (Sign), i.e., activation (1) or inhibition (-1).
### Display the inferred transcription factor activity in a network.
# Parameter 'weightCut' removes edges, i.e., transcription factor activities lower than the chosen
# threshold I=[0,100].
# Parameter 'clusterSize' only allows clusters with size equal or bigger than threshold.
carnival.net <- networkCARNIVAL(carnival_result, weightCut = 75, clusterSize = 15)

# And plot the network. The plot is exported to pdf Format.
plotCnet(carnival.net, layoutType=1, f.title="data/output/CARNIVAL_C15",draw.legend=F, draw.description=F)

### Now we will investigate the activity of TFs in a functional context with gene set enrichment
# analysis.

## Here we retrieve the attributes of all nodes in the network, i.e., node type, activation status.
nodes <- data.frame(carnival_result$nodesAttributes,
                    stringsAsFactors = F) %>% 
  mutate(act_sign = sign(as.numeric(AvgAct)))

# Are active nodes enriched for a function?
# Extract the active nodes and perform gene-set enrichment (GSE) to determine significantly enriched pathways.
active_nodes = nodes %>% filter(act_sign == 1)
# GSE and Hypergeometric test for the active nodes.
gse_active <- GSE_analysis(geneList = active_nodes[[1]],
                           Annotation_DB = onco_gsets)

## Are inactive nodes enriched for a function?
# Extract the inactive nodes and perform gene-set enrichment (GSE) to determine significantly enriched pathways.
inactive_nodes = nodes %>% filter(act_sign == -1)
# GSE and Hypergeometric test for the inactive nodes.
gse_inactive <- GSE_analysis(geneList = inactive_nodes[[1]],
                             Annotation_DB = onco_gsets)

## Set variables with type of p-value and threshold for significance.
p.select <- "p_value"
p.threshold <- 0.05
# And have a look, i.e., how many significant gene sets remain with these parameters.
gse_active %>% dplyr::filter( get(p.select) < p.threshold ) %>% nrow
gse_inactive %>% dplyr::filter( get(p.select) < p.threshold ) %>% nrow

## After we have chosen the correct values we will visualize the most significant pathways.
p.active <- plotGSE(gse_active, p.title = "Active nodes ", p.select, p.threshold)
p.inactive <- plotGSE(gse_inactive, p.title = "Inactive nodes", p.select, p.threshold)
# And write to file
pdf(file = "data/output/CARNIVAL_TF_OncogenicGeneSets_late_vs_early.pdf", height = 18, width = 18, useDingbats = F, onefile=TRUE)
p.active | p.inactive
dev.off()

# We extract our TF (Dorothea) results to take a look at the distribution of active and inactive nodes.
TF.stats <- TFActivity_results$results %>%
  rownames_to_column(var = "Gene")
  
# Are active nodes upregulated?
useful_nodes <- nodes %>% filter(act_sign != 0) %>%
  left_join(TF.stats, by = c("Node"="Gene"))

# prepare Violin plot with statistical test on active/inactive nodes distribution.
p.nodes <- ggpubr::ggviolin(useful_nodes, x="act_sign", y = "logFC", fill = "act_sign",
                            alpha = 0.75,
                            palette = c("#D62728", "#7F7F7F"),
                            add = c("jitter", "boxplot"),
                            add.params = list(fill = "white"),
                            trim = T,
                            title	= "Dorothea activation status",
                            subtitle = "airways Late vs Early",
                            xlab = "Activation direction",
                            ylab = "log2 FoldChange") +
  ggpubr::stat_compare_means() +
  theme_bw()
# Move the plots legend to the bottom for a ncer plot.
p.nodes <- ggpubr::ggpar(p.nodes, legend="bottom")

# Only save this plot to file.
pdf(file = "data/output/CARNIVAL_TF_regulation_late_vs_early.pdf.pdf", height = 6, width = 6, useDingbats = F, onefile=TRUE)
p.nodes
dev.off()

# Save all plots to file
pdf(file = "data/output/CARNIVAL_TF_regulation_OncoGeneSets_late_vs_early.pdf.pdf", height = 6, width = 6, useDingbats = F, onefile=TRUE)
ggpubr::ggarrange(p.nodes,ggpubr::ggarrange(p.active, p.inactive, nrow=2))
dev.off()

summary( lm(logFC ~ factor(act_sign), data = useful_nodes) )

## And export results to table format.
myres_list[["ONCO activated Late vs Early"]] <- gse_active
myres_list[["ONCO in-activated Late vs Early"]] <- gse_inactive

## And finally export all the results in an excel table.
openxlsx::write.xlsx(myres_list,file = paste("data/output/airways.PDC.LATEvsEARLY.xlsx",sep=""),
                     rowNames=TRUE,colNames=TRUE, colWidths = c(NA, "auto", "auto"))



















