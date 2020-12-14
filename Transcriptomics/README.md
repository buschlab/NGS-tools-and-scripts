# Transcriptomics


##### A collection of pipelines to *process* RNAseq raw-data and *analyse* the output.

#### Analyses
* **Differential Expression Analysis** (DEA) from count table (<a href="https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html" target="_blank">DESeq2</a>) an visualization in Heatmap, PCA and Volcano plot.
* **Single Sample Gene Set Enrichment Analysis** (ssGSEA) on normalized count data (<a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7" target="_blank">gsva</a>) and linear model fitting (<a href="https://kasperdanielhansen.github.io/genbioconductor/html/limma.html" target="_blank">limma</a>). Visualize interactions (gene set overlap) of enriched pathways as networks (<a href="https://igraph.org" target="_blank">igraph</a>).
* **PROGENy, DoRothEA, CARNIVAL pipeline** on normalized count data.
	1. <a href="https://saezlab.github.io/progeny/" target="_blank">PROGENy: Pathway RespOnsive GENes for activity inference</a> and <a href="https://kasperdanielhansen.github.io/genbioconductor/html/limma.html" target="_blank">limma</a> on activity scores.
	2. <a href="https://saezlab.github.io/dorothea/articles/single_cell_vignette.html">TF activity inference from scRNA-seq data with DoRothEA as regulon resource.
</a> Subsequent linear modeling (<a href="https://kasperdanielhansen.github.io/genbioconductor/html/limma.html" target="_blank">limma</a>).
	3. <a href="https://saezlab.github.io/CARNIVAL/">CARNIVAL: CAusal Reasoning for Network identification using Integer VALue programming
</a> on normalized count data and the outputs of PROGENy and DoRothEA. Linear model fitting (<a href="https://kasperdanielhansen.github.io/genbioconductor/html/limma.html" target="_blank">limma</a>) and network visualization (<a href="https://igraph.org" target="_blank">igraph</a>).
+ Visualization of enriched pathways on the basis of STRING Protein-Protein-Interaction (PPI) networks with the option to display significant genes from DEA.