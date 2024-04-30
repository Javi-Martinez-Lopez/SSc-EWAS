library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

## ----"load data"------------------------------------------------------------------------------------------------------
inputs_dir <- system.file("extdata", package = "decoupleR")
data <- readRDS(file.path(inputs_dir, "bk_data.rds"))

## ----"counts"---------------------------------------------------------------------------------------------------------
# Remove NAs and set row names
counts <- data$counts %>%
  dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x)) %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix()
head(counts)

## ----"design"---------------------------------------------------------------------------------------------------------
design <- data$design
design

## ----"deg"------------------------------------------------------------------------------------------------------------
DEGs <- read.delim("~/PRECISESADS/Meth/JML_methylation/eQTMs/DEGs.tsv", header=FALSE)
DEGs$t <- -log10(DEGs$V4)*DEGs$V5 
colnames(DEGs) <- c('Number','ID','P.Value.NotAdjusted','P.Value','logFC','F','t')
DEGs$ID <- make.names(DEGs$ID)
deg <- DEGs[,c(2,5,7,4)] %>% 
  filter(!ID%in%c('AC012652.1', 'AC138035.1', 'AL592284.1', 'X5S_rRNA', 'X7SK', 'Y_RNA'))%>%
  filter(!is.na(t)) %>% 
  column_to_rownames(var = 'ID') %>%
  as.matrix()
head(deg)

## ----"collectri"------------------------------------------------------------------------------------------------------
net <- get_collectri(organism='human', split_complexes=FALSE)
net

## ----"contrast_umean", message=FALSE----------------------------------------------------------------------------------
# Run umean
contrast_acts <- run_ulm(mat=deg[, 't', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)
                         
SIGNIF.TFs <- contrast_acts[contrast_acts$p_value<0.05&abs(contrast_acts$score)>1,]


