#load the libraries 
library(SeuratObject)
library(SeuratDisk)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

---
title: "data filtering"
output: html_document
date: '2024-04-08'
---

# To save the entire workspace:
save.image(file = "my_workspace.RData")
load("my_workspace.RData")

#clear env:
# rm(list = ls())


library(SeuratObject)
library(SeuratDisk)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)


#loom_data <- SeuratDisk::Connect(filename = "/data/gpfs/projects/punim2183/data_processed/l5_All.loom", mode = 'r')
#mouse_cortex_data <- as.Seurat(loom_data)
#saveRDS(mouse_cortex_data, file = '/data/gpfs/projects/punim2183/data_processed/mouse_cortex_data.rds')

mouse_cortex_data <- readRDS("/data/gpfs/projects/punim2183/data_processed/mouse_cortex_scRNAseq.rds")

mouse_cortex_data <- NormalizeData(mouse_cortex_data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

mouse_cortex_data <- RunUMAP(mouse_cortex_data, reduction = 'pca', dims = 1:25, assay = 'RNA', 
                             reduction.name = "rna_umap", reduction.key = "RNA_umap_")





