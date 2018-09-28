# NGS_scripts

## _genehancer.R

**'importGH'** *imports the datadump.xlsx from 
https://genecards.weizmann.ac.il/geneloc/index.shtml and returns adjusted ghDB object*
@param: path, filepath on your system
@return: ghDB, genehancer-DB import to work with

**'queryGHdb'** *given a list of GH-regions with START/END coordinates this function maps 
variants in set to this list and creates a new column "GHID"*
@param: set, annotated list to be matched
@param: GHdb, list of GH-regions with cols "GHID","CHROM","START","END"
@return: set, input-set extended by "CandidateGene"-column

**'matchGenehancers'** *given a list of names for gois this function finds all the regions
in ghDB that are associated with those genes and outputs a reduced version of ghDB to
work with in your project*
@param: query, list of gene names to match against
@param: genehancer, imported dataframe of ghDB
@return: tmpGH, ghDB reduced to goi

**'evalGH'** *analyses cols 'CandidateGene' and 'GHID' and returns summary for input-set*
@param: input, annotated list (SV,CNV,INDEL,SNV)
@param: ghDB, genehancer-DB object 
@return: summary object with cols 'Gene','scores','sum|max|minscores','nscores'

**'topGGH'** *takes an evaluation object as input and returns a bar-plot object of the top
genes of interest associated with GH-regions, i.e., most GH-regions*
@param: evalObj, output of 'evalGH' function
@param: percent, displayed percentage of genes [0,1], defaults to 10% if left empty
@param: plottitle, headline for this plot, may be left empty
@return: p1, a ggplot2 prepped plot

**'topGScore'** *takes an evaluation object as input and returns a bar-plot object of the top
genes of interest based on max-score*
@param: evalObj, output of 'evalGH' function
@param: percent, displayed percentage of genes [0,1]1], defaults to 10% if left empty
@param: plottitle, headline for this plot, may be left empty
@return: p1, a ggplot2 prepped plot

**'topGDA'** *analyses cols 'CandidateGene' return plot of directly associated Cgenes*
@param: input, annotated list (SV,CNV,INDEL,SNV)
@param: percent, displayed percentage of genes [0,1]1], defaults to 10% if left empty
@param: plottitle, headline for this plot, may be left empty
@return: p1, a ggplot2 prepped plot



