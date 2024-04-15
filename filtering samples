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
```{r}
# To save the entire workspace:
save.image(file = "my_workspace.RData")
load("my_workspace.RData")

#clear env:
# rm(list = ls())

```

```{r}
library(SeuratObject)
library(SeuratDisk)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

```

```{r}
#loom_data <- SeuratDisk::Connect(filename = "/data/gpfs/projects/punim2183/data_processed/l5_All.loom", mode = 'r')
#mouse_cortex_data <- as.Seurat(loom_data)
#saveRDS(mouse_cortex_data, file = '/data/gpfs/projects/punim2183/data_processed/mouse_cortex_data.rds')

mouse_cortex_data <- readRDS("/data/gpfs/projects/punim2183/data_processed/mouse_cortex_scRNAseq.rds")

mouse_cortex_data <- NormalizeData(mouse_cortex_data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

mouse_cortex_data <- RunUMAP(mouse_cortex_data, reduction = 'pca', dims = 1:25, assay = 'RNA', 
                             reduction.name = "rna_umap", reduction.key = "RNA_umap_")

```

```{r}
# keep samples with a certain number of Astrocytes, Neurons and Oligos as well
min_astrocyte_count <- 200 
min_neuron_count <- 200  # define the minimum for Neurons
min_oligo_count <- 200   # define the minimum for Oligos

# Summarize cell counts by SampleID and Class
cell_counts <- mouse_cortex_data@meta.data %>%
  group_by(SampleID, Class) %>%
  summarise(CellCount = n(), .groups = 'drop')
as.data.frame(cell_counts)

# Filter the samples with the minimum number of Cells of interest
astrocyte_counts <- cell_counts %>%
  filter(Class == 'Astrocytes' & CellCount >= min_astrocyte_count) %>%
  select(SampleID)  

sum(table(astrocyte_counts))

neuron_counts <- cell_counts %>%
  filter(Class == 'Neurons' & CellCount >= min_neuron_count) %>%
  select(SampleID)

sum(table(neuron_counts))

oligo_counts <- cell_counts %>%
  filter(Class == 'Oligos' & CellCount >= min_oligo_count) %>%
  select(SampleID)

sum(table(oligo_counts))

# Combine the filters to get a list of SampleIDs to keep
samples_to_keep <- Reduce(intersect, list(astrocyte_counts$SampleID, neuron_counts$SampleID, oligo_counts$SampleID))

# Subset the Seurat object to include only cells from the samples that meet the criteria
seurat_subset <- subset(mouse_cortex_data, subset = SampleID %in% samples_to_keep)

# Cell Type per sample
target_cell_counts_by_sample <- table(seurat_subset@meta.data$SampleID, seurat_subset@meta.data$Class)
print(target_cell_counts_by_sample)

# Convert the table to a data frame for plotting
target_counts_df <- as.data.frame(target_cell_counts_by_sample)
names(target_counts_df) <- c("Sample", "CellType", "Count")

# Create a more readable bar plot

ggplot(target_counts_df, aes(x = Sample, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5),
    legend.position = "bottom"
  ) +
  labs(x = "Sample", y = "Cell Count", fill = "Cell Type", title = "Cell Counts by Sample and Type")


# Plot the data in a stacked bar chart
ggplot(target_counts_df, aes(x = Sample, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Sample", y = "Cell Count", fill = "Cell Type", title = "Cell Type Frequencies by Sample")





level of gene expression in those samples/ cells of interest (violin plot), filter out some samples on that basis.



Idents(seurat_subset) <- seurat_subset$Class

# List all the identities in the Seurat object
levels(Idents(seurat_subset))

Violn_plot <- VlnPlot(seurat_subset, features = "Ntrk2", group.by = 'Class', idents = c("Neurons","Astrocytes", "Oligos"), sort = TRUE, pt.size = 0.1) + ylim(c(0,5))+ NoLegend()
Violn_plot


# Ntrk2 expression

ntrk2_expression <- FetchData(seurat_subset, vars = "Ntrk2")
ntrk2_expression


# First, check if the gene is in the RNA assay data
if (!("Ntrk2" %in% rownames(seurat_subset[["RNA"]]@counts))) {
  stop("Gene not found in the dataset")
}

# The subset function expects a logical vector, so create one first:
cells_to_keep <- seurat_subset@assays$RNA@counts["Ntrk2", ] > 0.8
cells_to_keep_df <- as.data.frame(cells_to_keep) #this created a column of True/False

cell_names_to_keep <- names(cells_to_keep)[cells_to_keep] #fetching the cell-barcodes of those cells

# filter the metadata to get only the SampleIDs corresponding to those cell names
metadata_to_keep <- seurat_subset@meta.data[cell_names_to_keep, ]

# kept samples:
sample_ids_to_keep <- unique(metadata_to_keep$SampleID)
# sample_ids_to_keep


# Use the %in% operator to filter cells from the desired SampleIDs
cells_to_subset <- rownames(metadata_to_keep[metadata_to_keep$SampleID %in% samples_to_keep, ])

# or
## sample_ids_to_keep <- unique(metadata_to_keep$SampleID)
## seurat_filtered_by_sample <- subset(seurat_subset, subset = SampleID %in% sample_ids_to_keep)

# Subset Seurat object to include only cells with high Ntrk2 expression
# seurat_filtered_by_ntrk2 <- subset(seurat_subset, cells = cell_names_to_keep)

# subset the Seurat object to keep only the cells from the desired SampleIDs with high Ntrk2 expression
seurat_filtered_by_ntrk2 <- subset(seurat_subset, cells = cells_to_subset)
ntrk2_high_expression <- FetchData(seurat_filtered_by_ntrk2, vars = "Ntrk2")
ntrk2_high_expression


high_ntrk_Violn_plot <- VlnPlot(seurat_filtered_by_ntrk2, features = "Ntrk2", group.by = 'Class', idents = c("Neurons","Astrocytes", "Oligos"), sort = TRUE, pt.size = 0.1) + ylim(c(0,5))+ NoLegend()
high_ntrk_Violn_plot



ntrk2_high_expression <- FetchData(seurat_filtered_by_ntrk2, vars = "Ntrk2")
ntrk2_high_expression


Ntrk2_in_targeted <- Seurat::FeaturePlot(seurat_filtered_by_ntrk2, features = 'Ntrk2', reduction = 'rna_umap', max.cutoff = 2) 
Ntrk2_in_targeted

