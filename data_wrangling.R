if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")
setwd('~/bioinfo_final_project/')
### LOAD LIBRARIES
library(airway)
library(tidyverse)

data('airway') #First lets load the airway dataset in
airway #lets see some information about the dataset

### EXTRACTING DATA AND SAVING TO A CSV
sample_information <- as.data.frame(colData(airway)) #extracts info from colData list in airway
sample_information <- sample_info[,c(2,3)] #extracts only the second and third column
sample_information$dex <- gsub('trt', 'treated', sample_information$dex)
sample_infomation$dex <- gsub('untrt', 'untreated', sample_information$dex)
names(sample_information) <- c('cellLine', 'dexamethasone')
write.table(sample_information, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)
