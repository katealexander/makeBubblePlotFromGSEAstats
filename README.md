# Description
This repository describes how I made bubble plots of GSEA-derived statistics for each cancer type analyzed for speckle signature patient group gene set enrichment in The Cancer Genome Atlas. 

# Input files
To generate gene rank files for GSEA analysis, I used the Python script "makeRnkFiles.py" to calculate the median expression of each gene in speckle Signature I and II patients and exported a files "CANCER_signatureRatio.rnk" with a list of genes and the log2ratio of median expression between Signature I and II. These ".rnk" files, included in this repository, were used as inputs for GSEA analysis

# GSEA analysis
I used the Desktop version of GSEA_4.3.2 to perform GSEA analysis on the Hallmarks (H), Positional (C1), and KEGG (C2:CP:KEGG) gene sets. In hindsight, it would have been more efficient to do this in the command line version. 

# Make a bubble plot in R
The following Rscripts create bubble plots for all the cancer types to summarize the GSEA statistics. They refer to where I stored my GSEA output as well as my naming convention for the folders. These parts will have to be edited based on where GSEA results are stored and what the folder names are.

```Rscript GSEA_KEGG_dotplot.R``` makes "KEGG_dotplot.pdf"

```Rscript GSEA_hallmarks_dotplot.R``` makes "HALLMARK_dotplot.pdf" 

```Rscript chromosomeGSEAs.R``` creates a new directory with a bubble plot for each chromosome arm.


