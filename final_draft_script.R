### Link to paper ###
#https://pubmed.ncbi.nlm.nih.gov/24926665/

### Link to Data ###
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

### HOW TO INSTALL PACKAGES ###

# TIDYVERSE
install.packages("tidyverse") 

# DESEQ2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# AIRWAY
BiocManager::install("airway")

# LOAD IN LIBRARIES
library(tidyverse)
library(DESeq2)
library(airway)
library(BiocGenerics)

# SET YOUR WORKING DIRECTORY
setwd("~/bioinfo_final_project/") #MAKE SURE YOUR FILES AND SCRIPT ARE IN THE SAME DIRECTORY (FOLDER)

# READ IN DATAFRAMES

# READ IN COUNTS DF
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# READ IN SAMPLE INFO DF
colData <- read.csv('sample_info.csv')

# MAKE SURE THAT COLUMNS IN 'COUNTS' MATCHES ROWS IN SAMPLE INFO 
all(colnames(counts_data) %in% rownames(colData))

# MAKE SURE THEY ARE ALSO IN THE SAME ORDER
all(colnames(counts_data) == rownames(colData))

# MAKE A DESEQDATASET OBJECT
dds.object <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds.object

# REMOVE ROWS WITH LOW GENE COUNTS SO IN THIS CASE WE WANT TO KEEP ROWS WHO HAVE A GENE COUNT GREATER THAN 10
keep <- rowSums(counts(dds.object)) >= 10
dds.object <- dds.object[keep,]

dds.object

# SET FACTOR LEVEL
dds.object$dexamethasone <- relevel(dds.object$dexamethasone, ref = "untreated")

# RUN DESEQ
dds <- DESeq(dds.object)
res <- results(dds)

res

# VISUALIZE RESULTS

# MA plot
plotMA(res) +
  title('Human ASM Cells Treated with Dexamethasone vs. Untreated') 
  

# PCA
vst_data <- vst(dds, blind = FALSE)
plotPCA(vst_data, intgroup = 'dexamethasone') 



