#load the libraries 
library(SeuratObject)
library(SeuratDisk)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

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




---
title: "arranged"
output: html_document
date: '2024-05-25'
---

```{r setup, ech=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
source("/data/gpfs/projects/punim2183/data_processed/my_functions.R")

```

Summarizing cell counts for each cell type and calculating the proportion of each cell type for scaling

```{r}
cell_types_of_interest <- c("Neurons", "Astrocytes", "Oligos")
gene_of_interest <- "Ntrk2"

```


```{r}
cell_counts <- mouse_cortex_data@meta.data %>%
    group_by(SampleID) %>%
    summarise(TotalCellCount = n(), .groups = 'drop') %>%
    as.data.frame()
summary(cell_counts)

neurons_cell_counts <- mouse_cortex_data@meta.data %>%
    dplyr::filter(Class == "Neurons") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(NeuronsCellCount = n(), .groups = 'drop') %>%
    as.data.frame()
head(neurons_cell_counts)
summary(neurons_cell_counts)


astrocytes_cell_counts <- mouse_cortex_data@meta.data %>%
    dplyr::filter(Class == "Astrocytes") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(AstrocytesCellCount = n(), .groups = 'drop') %>%
    as.data.frame()
head(astrocytes_cell_counts)
summary(astrocytes_cell_counts)

oligos_cell_counts <- mouse_cortex_data@meta.data %>%
    dplyr::filter(Class == "Oligos") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::summarise(OligosCellCount = n(), .groups = 'drop') %>%
    as.data.frame()
head(oligos_cell_counts)
summary(oligos_cell_counts)


merged_counts <- cell_counts %>%
    left_join(neurons_cell_counts, by = "SampleID") %>%
    left_join(astrocytes_cell_counts, by = "SampleID") %>%
    left_join(oligos_cell_counts, by = "SampleID")
merged_counts[is.na(merged_counts)] <- 0

merged_counts <- merged_counts %>%
    mutate(
        NeuronsProportion = NeuronsCellCount / TotalCellCount,
        AstrocytesProportion = AstrocytesCellCount / TotalCellCount,
        OligosProportion = OligosCellCount / TotalCellCount
    )


merged_counts <- merged_counts %>%
    rowwise() %>%
    mutate(
        SumProportions = NeuronsProportion + AstrocytesProportion + OligosProportion,
        NeuronsProportion = NeuronsProportion / SumProportions,
        AstrocytesProportion = AstrocytesProportion / SumProportions,
        OligosProportion = OligosProportion / SumProportions
    ) %>%
    select(-SumProportions)

head(merged_counts)

summary(merged_counts)

```

selecting samples with minimum cell count of each type:

```{r}

# Function to calculate overall proportions of each cell type of interest in the entire dataset

calculate_overall_proportions <- function(data) {
  
  total_cell_count <- data@meta.data %>%
    summarise(TotalCellCount = n()) %>%
    pull(TotalCellCount)
  
  cell_type_counts <- data@meta.data %>%           # Calculate cell counts for each cell type
    group_by(Class) %>%
    filter(Class == c("Neurons", "Astrocytes", "Oligos")) %>%
    summarise(CellCount = n(), .groups = 'drop')
  
  cell_type_counts <- cell_type_counts %>%          # Calculate proportions for each cell type
    mutate(Proportion = CellCount / total_cell_count)
  
  return(as.data.frame(cell_type_counts))
}


prop <- calculate_overall_proportions(mouse_cortex_data)

prop


```

```{r}
library(dplyr)

# Function to get the proportion for a specific cell type
get_proportion <- function(proportions_df, cell_type) {
  proportion <- proportions_df %>%
    filter(Class == cell_type) %>%
    pull(Proportion)
  
  return(proportion)
}


neuron_proportion <- get_proportion(prop, "Neurons")
print(neuron_proportion)

astrocyte_proportion <- get_proportion(prop, "Astrocytes")
print(astrocyte_proportion)

oligo_proportion <- get_proportion(prop, "Oligos")
print(oligo_proportion)


```

```{r}


# Function to calculate overall proportions of each cell type in the entire dataset

calculate_overall_proportions <- function(data) {
  filtered_data <- data@meta.data %>%
    filter(Class %in% c("Astrocytes", "Neurons", "Oligos"))
  
  total_cell_count <- filtered_data %>%
    summarise(TotalCellCount = n()) %>%
    pull(TotalCellCount)
  
  cell_type_counts <- filtered_data %>%   # Calculate cell counts for each cell type
    group_by(Class) %>%
    summarise(CellCount = n(), .groups = 'drop')
  
  cell_type_counts <- cell_type_counts %>%   # Calculate proportions for each cell type
    mutate(Proportion = CellCount / total_cell_count)
  
  return(cell_type_counts)
}

# Function to get the proportion for a specific cell type
get_proportion <- function(proportions_df, cell_type) {
  proportion <- proportions_df %>%
    filter(Class == cell_type) %>%
    pull(Proportion)
  
  return(proportion)
}


# Function to filter cells by minimum count according to overall proportions
filter_samples_by_min_count <- function(data, min_total_count) {
  
  overall_proportions <- calculate_overall_proportions(data)

  # Scale minimum counts for each cell type based on overall proportions
  min_astrocyte_count <- min_total_count * get_proportion(overall_proportions, "Astrocytes")
  min_neuron_count <- min_total_count * get_proportion(overall_proportions, "Neurons")
  min_oligo_count <- min_total_count * get_proportion(overall_proportions, "Oligos")
  
  # Initialize a list to store SampleIDs that meet the criteria for each cell type
  valid_sample_ids <- list()
  
  # Check each cell type
  for (cell_type in c("Astrocytes", "Neurons", "Oligos")) {
    min_count <- switch(cell_type,
                        "Astrocytes" = min_astrocyte_count,
                        "Neurons" = min_neuron_count,
                        "Oligos" = min_oligo_count)
    
    ids <- data@meta.data %>%
      filter(Class == cell_type) %>%
      group_by(SampleID) %>%
      summarise(CellCount = n(), .groups = 'drop') %>%
      filter(CellCount >= min_count) %>%
      pull(SampleID)
    
    if (length(ids) == 0) {
      cat(sprintf("No samples have the defined minimum count for %s.\n", cell_type))
    } else {
      valid_sample_ids[[cell_type]] <- ids
    }
  }
  
  # Check if there are any SampleIDs that meet the criteria for all specified cell types
  if (length(valid_sample_ids) < 3) {
    cat("Some cell types did not meet the minimum cell count criteria.\n")
    return(NULL)
  }
  
  # Find common SampleIDs across all cell types
  samples_to_keep <- Reduce(intersect, valid_sample_ids)
  
  # Check if any samples meet the criteria
  if (length(samples_to_keep) == 0) {
    cat("No samples meet the minimum cell count across all specified cell types.\n")
    return(NULL)
  }
  
  # Subset the Seurat object to include only samples that meet the criteria
  seurat_subset <- subset(data, subset = SampleID %in% samples_to_keep & Class %in% c("Astrocytes", "Neurons", "Oligos"))
  
  # Verify if subset is empty
  if (seurat_subset@meta.data %>% nrow() == 0) {
    cat("The resulting subset of the data has no cells. Check the filtering criteria.\n")
    return(NULL)
  }
  
  # Cell Type per sample
  target_cell_counts_by_sample <- table(seurat_subset@meta.data$SampleID, seurat_subset@meta.data$Class)
  print(target_cell_counts_by_sample)
  
  # Convert the table to a data frame for plotting
  target_counts_df <- as.data.frame(target_cell_counts_by_sample)
  names(target_counts_df) <- c("Sample", "CellType", "Count")
  target_counts_df <- target_counts_df %>%
    filter(CellType %in% c("Astrocytes", "Oligos", "Neurons"))
  
  return(target_counts_df)
}

# Example usage
# Assuming `mouse_cortex_data` is your Seurat object
min_total_count <- 1000  # Example total minimum count to scale by proportions
filtered_data <- filter_samples_by_min_count(mouse_cortex_data, min_total_count)


```


```{r}
plot_cell_type_counts <- function(counts_df, df_name) {
  title_text <- sprintf("Cell Type Counts per Sample\n(Data Frame: %s)", df_name)

  ggplot(counts_df, aes(x = Sample, y = Count, fill = CellType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = title_text,
         x = "Sample",
         y = "Cell Count",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}


```

```{r}
min_500 <- filter_samples_by_min_count(mouse_cortex_data, 500)

min_750 <- filter_samples_by_min_count(mouse_cortex_data, 750)

min_1000 <- filter_samples_by_min_count(mouse_cortex_data, 1000)

min_1200 <- filter_samples_by_min_count(mouse_cortex_data, 1200)

min_1500 <- filter_samples_by_min_count(mouse_cortex_data, 1500)


```


```{r}
p_500 <- plot_cell_type_counts(min_500, "minimum 500")
p_750 <- plot_cell_type_counts(min_750,"minimum 750")
p_1000 <- plot_cell_type_counts(min_1000, "minimum 1000")
p_1200 <- plot_cell_type_counts(min_1200, "minimum 1200")

p_500
p_750
p_1000
p_1200
```

```{r plot_chunck, fig.width=15, fig.height=8}

# Add a new column 'Threshold' to each dataframe
min_500$Threshold <- "500"
min_750$Threshold <- "750"
min_1000$Threshold <- "1000"
min_1200$Threshold <- "1200"


# Combine the dataframes
combined_df <- rbind(min_500, min_750, min_1000, min_1200)
combined_df$Threshold <- factor(combined_df$Threshold)

combined_df <- combined_df %>%
  mutate(Threshold = factor(Threshold, levels = c(500, 750, 1000, 1200)))

plot <- ggplot(combined_df, aes(x = Sample, y = Count)) +
  geom_boxplot(aes(fill = as.factor(Threshold))) +
  facet_wrap(~ Threshold, scales = "free_y") +
  labs(title = "Cell Counts by Threshold Across Samples",
       x = "Sample",
       y = "Cell Count",
       fill = "Threshold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)
```


```{r}
cell_types_of_interest <- c("Neurons", "Astrocytes", "Oligos")
by_gene_1_50 <- filter_samples_by_gene_expression(mouse_cortex_data, "Ntrk2", 1.0, 50, cell_types_of_interest)
by_gene_2_50 <- filter_samples_by_gene_expression(mouse_cortex_data, "Ntrk2", 2.0, 50, cell_types_of_interest)
by_gene_1_80 <- filter_samples_by_gene_expression(mouse_cortex_data, "Ntrk2", 1.0, 80, cell_types_of_interest)
by_gene_2_80 <- filter_samples_by_gene_expression(mouse_cortex_data, "Ntrk2", 2.0, 80, cell_types_of_interest)

by_gene_1_50
by_gene_2_50
by_gene_1_50
by_gene_1_80
# Print all unique combinations of SampleID and Class
unique_combinations <- unique(mouse_cortex_data@meta.data %>% dplyr::select(SampleID, Class))
#print(unique_combinations)

```

```{r, fig.width=15, fig.height=8}
# Calculate intersections for different conditions
intersected_samples_500a <- intersect(min_500$Sample, by_gene_1_50$SampleID)
intersected_samples_500b <- intersect(min_500$Sample, by_gene_2_50$SampleID)
intersected_samples_500c <- intersect(min_500$Sample, by_gene_1_80$SampleID)
intersected_samples_500d <- intersect(min_500$Sample, by_gene_2_80$SampleID)

intersected_samples_750a <- intersect(min_750$Sample, by_gene_1_50$SampleID)
intersected_samples_750b <- intersect(min_750$Sample, by_gene_2_50$SampleID)
intersected_samples_750c <- intersect(min_750$Sample, by_gene_1_80$SampleID)
intersected_samples_750d <- intersect(min_750$Sample, by_gene_2_80$SampleID)

intersected_samples_1000a <- intersect(min_1000$Sample, by_gene_1_50$SampleID)
intersected_samples_1000b <- intersect(min_1000$Sample, by_gene_2_50$SampleID)
intersected_samples_1000c <- intersect(min_1000$Sample, by_gene_1_80$SampleID)
intersected_samples_1000d <- intersect(min_1000$Sample, by_gene_2_80$SampleID)

intersected_samples_1200a <- intersect(min_1200$Sample, by_gene_1_50$SampleID)
intersected_samples_1200b <- intersect(min_1200$Sample, by_gene_2_50$SampleID)
intersected_samples_1200c <- intersect(min_1200$Sample, by_gene_1_80$SampleID)
intersected_samples_1200d <- intersect(min_1200$Sample, by_gene_2_80$SampleID)

# Create a data frame for plotting
plot_data <- data.frame(
  CellCountThreshold = rep(c("500", "750", "1000", "1200"), each = 4),
  GeneExpressionCondition = rep(c("expression_thresholds_1.0_50%", "expression_thresholds_2.0_50%", "expression_thresholds_1.0_80%", "expression_thresholds_2.0_80%"), times = 4),
  SampleCount = c(
    length(intersected_samples_500a), length(intersected_samples_500b), length(intersected_samples_500c), length(intersected_samples_500d),
    length(intersected_samples_750a), length(intersected_samples_750b), length(intersected_samples_750c), length(intersected_samples_750d),
    length(intersected_samples_1000a), length(intersected_samples_1000b), length(intersected_samples_1000c), length(intersected_samples_1000d),
    length(intersected_samples_1200a), length(intersected_samples_1200b), length(intersected_samples_1200c), length(intersected_samples_1200d)
  ),
  SampleIDs = I(list(
    intersected_samples_500a, intersected_samples_500b, intersected_samples_500c, intersected_samples_500d,
    intersected_samples_750a, intersected_samples_750b, intersected_samples_750c, intersected_samples_750d,
    intersected_samples_1000a, intersected_samples_1000b, intersected_samples_1000c, intersected_samples_1000d,
    intersected_samples_1200a, intersected_samples_1200b, intersected_samples_1200c, intersected_samples_1200d
  ))
)

# Ensure the GeneExpressionCondition is a factor with all levels specified
plot_data$GeneExpressionCondition <- factor(plot_data$GeneExpressionCondition, levels = c("expression_thresholds_1.0_50%", "expression_thresholds_2.0_50%", "expression_thresholds_1.0_80%", "expression_thresholds_2.0_80%"))

# Ensure the CellCountThreshold is a factor with the correct order
plot_data$CellCountThreshold <- factor(plot_data$CellCountThreshold, levels = c("500", "750", "1000", "1200"))
print(plot_data)

line_plot <- ggplot(plot_data, aes(x = CellCountThreshold, y = SampleCount, group = GeneExpressionCondition, color = GeneExpressionCondition)) +
  geom_line() +
  geom_point() +
  labs(title = "Number of Samples by Intersection Conditions", x = "Cell Count Threshold", y = "Number of Samples") +
  theme_minimal()

print(line_plot)

line_plot + facet_wrap(~ GeneExpressionCondition)


```
```{r, include =FALSE}
plot_non_zero <- ggplot(non_zero_table, aes(x = CellType, y = NonZeroCount, color = Dataset, group = Dataset)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Non-Zero Count of Ntrk2 Expression", x = "Cell Type", y = "Non-Zero Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Dataset)

# Print the Non-Zero Count line plot with facets
print(plot_non_zero)

# Plot Percentage as a line plot with facets
plot_percentage <- ggplot(non_zero_table, aes(x = CellType, y = Percentage, color = Dataset, group = Dataset)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Percentage of Non-Zero Ntrk2 Expression", x = "Cell Type", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Dataset)

# Print the Percentage line plot with facets
print(plot_percentage)

```


```{r}

# Store the intersected samples in a list
intersected_samples_list <- list(
  intersected_samples_500a, intersected_samples_500b, intersected_samples_500c, intersected_samples_500d,
  intersected_samples_750a, intersected_samples_750b, intersected_samples_750c, intersected_samples_750d,
  intersected_samples_1000a, intersected_samples_1000b, intersected_samples_1000c, intersected_samples_1000d,
  intersected_samples_1200a, intersected_samples_1200b, intersected_samples_1200c, intersected_samples_1200d
)

# Find the maximum length of the intersected sample lists
max_length <- max(sapply(intersected_samples_list, length))

# Pad each list with NA values to match the maximum length
padded_samples_list <- lapply(intersected_samples_list, function(x) {
  length(x) <- max_length
  return(x)
})

# Combine the padded lists into a data frame
combined_intersected <- as.data.frame(do.call(cbind, padded_samples_list))

# Assign column names based on the conditions
colnames(combined_intersected) <- c(
  "500_expression_thresholds_1.0_50%", "500_expression_thresholds_2.0_50%", "500_expression_thresholds_1.0_80%", "500_expression_thresholds_2.0_80%",
  "750_expression_thresholds_1.0_50%", "750_expression_thresholds_2.0_50%", "750_expression_thresholds_1.0_80%", "750_expression_thresholds_2.0_80%",
  "1000_expression_thresholds_1.0_50%", "1000_expression_thresholds_2.0_50%", "1000_expression_thresholds_1.0_80%", "1000_expression_thresholds_2.0_80%",
  "1200_expression_thresholds_1.0_50%", "1200_expression_thresholds_2.0_50%", "1200_expression_thresholds_1.0_80%", "1200_expression_thresholds_2.0_80%"
)

# Print the combined_intersected data frame to check its contents
print(combined_intersected)


```

```{r}
intersected_samples_list_1 <- list(
  intersected_samples_750b,intersected_samples_750c, intersected_samples_1000a, intersected_samples_1000b,
  intersected_samples_1200a, intersected_samples_1200b
)

# Find the maximum length of the intersected sample lists
max_length <- max(sapply(intersected_samples_list_1, length))

# Pad each list with NA values to match the maximum length
padded_samples_list_1 <- lapply(intersected_samples_list_1, function(x) {
  length(x) <- max_length
  return(x)
})

# Combine the padded lists into a data frame
combined_intersected_1 <- as.data.frame(do.call(cbind, padded_samples_list_1))

# Assign column names based on the conditions
colnames(combined_intersected_1) <- c("750b_expression_thresholds_2.0_50%",
  "750_expression_thresholds_1.0_80%", "1000_expression_thresholds_1.0_50%", 
  "1000_expression_thresholds_2.0_50%", "1200_expression_thresholds_1.0_50%", 
  "1200_expression_thresholds_2.0_50%"
)

# Print the combined_intersected data frame to check its contents
print(combined_intersected_1)


#write.csv(combined_intersected_1, "combined_intersected_1.csv", row.names = FALSE)

#install.packages("writexl")

#write_xlsx(combined_intersected_1, "combined_intersected_1.xlsx")


```

```{r, Seurat subsetting}
chosen_samples_11 <- c("10X26_3", "10X26_4", "10X52_2","10X22_3", "10X22_4", "10X24_1", 
                    "10X52_1","10X52_3", "10X52_4", "10X57_2", "10X57_3")
chosen_samples_3 <- c("10X26_3", "10X26_4", "10X52_2")

eleven_samples_seurat <- subset(mouse_cortex_data, subset = SampleID %in% chosen_samples_11 & Class %in% cell_types_of_interest)

three_samples_seurat <- subset(mouse_cortex_data, subset = SampleID %in% chosen_samples_3 & Class %in% cell_types_of_interest)

ntrk2_expression_all <- FetchData(mouse_cortex_data, vars = "Ntrk2")
non_zero_all <- sum(ntrk2_expression_all$Ntrk2 > 0)


# mouse_cortex data to calculate expression in each cell type
subset_Neurons_all <- subset(mouse_cortex_data, subset = Class == "Neurons")
neurons_expression_all <- FetchData(subset_Neurons_all, vars = "Ntrk2")
non_zero_N_all <- sum(neurons_expression_all[["Ntrk2"]] > 0)

subset_Astrocytes_all <- subset(mouse_cortex_data, subset = Class == "Astrocytes")
astro_expression_all <- FetchData(subset_Astrocytes_all, vars = "Ntrk2")
non_zero_A_al <- sum(astro_expression_all[["Ntrk2"]] > 0)

subset_Oligos_all <- subset(mouse_cortex_data, subset = Class == "Oligos")
oligo_expression_all <- FetchData(subset_Oligos_all, vars = "Ntrk2")
non_zero_O_all <- sum(oligo_expression_all[["Ntrk2"]] > 0)



#subset of 11 samples to calculate expression in each cell type

ntrk2_expression_11 <- FetchData(samples_seurat, vars = "Ntrk2")
num_non_zero_11 <- sum(ntrk2_expression_11$Ntrk2 > 0)

subset_Neurons_11 <- subset(samples_seurat, subset = Class == "Neurons")
neurons_expression_11 <- FetchData(subset_Neurons_11, vars = "Ntrk2")
non_zero_N_11 <- sum(neurons_expression_11[["Ntrk2"]] > 0)

subset_Astrocytes_11 <- subset(eleven_samples_seurat, subset = Class == "Astrocytes")
astro_expression_11 <- FetchData(subset_Astrocytes_11, vars = "Ntrk2")
non_zero_A_11 <- sum(astro_expression_11[["Ntrk2"]] > 0)

subset_Oligos_11 <- subset(eleven_samples_seurat, subset = Class == "Oligos")
oligo_expression_11 <- FetchData(subset_Oligos_11, vars = "Ntrk2")
non_zero_O_11 <- sum(oligo_expression_11[["Ntrk2"]] > 0)


#subset of 3 samples to calculate expression in each cell type

ntrk2_expression_3 <- FetchData(three_samples_seurat, vars = "Ntrk2")
num_non_zero_3 <- sum(ntrk2_expression_3$Ntrk2 > 0)

subset_Neurons_3 <- subset(three_samples_seurat, subset = Class == "Neurons")
neurons_expression_3 <- FetchData(subset_Neurons_3, vars = "Ntrk2")
non_zero_N_3 <- sum(neurons_expression_3[["Ntrk2"]] > 0)

subset_Astrocytes_3 <- subset(three_samples_seurat, subset = Class == "Astrocytes")
astro_expression_3 <- FetchData(subset_Astrocytes_3, vars = "Ntrk2")
non_zero_A_3 <- sum(astro_expression_3[["Ntrk2"]] > 0)

subset_Oligos_3 <- subset(three_samples_seurat, subset = Class == "Oligos")
oligo_expression_3 <- FetchData(subset_Oligos_3, vars = "Ntrk2")
non_zero_O_3 <- sum(oligo_expression_3[["Ntrk2"]] > 0)

```

```{r}
ages <- three_samples_seurat@meta.data$Age[]
ages

# Remove the 'p' character and convert to numeric
numeric_ages <- as.numeric(sub("p", "", ages))

age_summary <- table(numeric_ages)

# Convert the result to a data frame
age_summary_df <- as.data.frame(age_summary)
colnames(age_summary_df) <- c("Age", "Count")

# Print the summary
print(age_summary_df)


```

```{r}
tissues <- three_samples_seurat@meta.data$Tissue
table(tissues)

tis <- mouse_cortex_data@meta.data$Tissue
table(tis)
```

```{r}
sex <- three_samples_seurat@meta.data$Sex
table(sex)


sex_all <- mouse_cortex_data@meta.data$Sex
table(sex_all)

```

```{r}

# Initialize a data frame to store the results
non_zero_table <- data.frame(
  Dataset = character(), 
  CellType = character(), 
  NonZeroCount = numeric(), 
  stringsAsFactors = FALSE
)

count_non_zero_expressions <- function(dataset, dataset_name, cell_type) {
  subset_seurat_object <- subset(dataset, subset = Class == cell_type)
  
  gene_expression <- FetchData(subset_seurat_object, vars = "Ntrk2")
  
  num_non_zero <- sum(gene_expression[["Ntrk2"]] > 0)
  
  # Add the result to the table
  non_zero_table <<- rbind(
    non_zero_table, 
    data.frame(
      Dataset = dataset_name, 
      CellType = cell_type, 
      NonZeroCount = num_non_zero, 
      stringsAsFactors = FALSE
    )
  )
}

# Define datasets, names, and cell types
datasets <- list(
  mouse_cortex_data = "mouse_cortex_data", 
  
  samples_seurat = "eleven_samples_seurat", 
 
  three_samples_seurat = "three_samples_seurat"
)

cell_types_of_interest <- c("Neurons", "Astrocytes", "Oligos")

# Loop through each dataset and cell type to count non-zero expressions
for (dataset_name in names(datasets)) {
  dataset <- get(dataset_name)
  
  # Total non-zero expressions for the entire dataset
  ntrk2_expression <- FetchData(dataset, vars = "Ntrk2")
  num_non_zero <- sum(ntrk2_expression$Ntrk2 > 0)
  non_zero_table <- rbind(
    non_zero_table, 
    data.frame(
      Dataset = dataset_name, 
      CellType = "Total", 
      NonZeroCount = num_non_zero, 
      stringsAsFactors = FALSE
    )
  )
  
  for (cell_type in cell_types_of_interest) {
    count_non_zero_expressions(dataset, dataset_name, cell_type)
  }
}

# Calculate percentage for each cell type within each dataset
non_zero_table <- non_zero_table %>%
  group_by(Dataset) %>%
  mutate(Percentage = (NonZeroCount / NonZeroCount[CellType == "Total"]) * 100)

# Print the resulting table
print(non_zero_table)

# Plot Non-Zero Count
plot_non_zero <- ggplot(non_zero_table, aes(x = CellType, y = NonZeroCount, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Non-Zero Count of Ntrk2 Expression", x = "Cell Type", y = "Non-Zero Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the Non-Zero Count plot
print(plot_non_zero)

# Plot Percentage
plot_percentage <- ggplot(non_zero_table, aes(x = CellType, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(title = "Percentage of Non-Zero Ntrk2 Expression", x = "Cell Type", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the Percentage plot
print(plot_percentage)
# write_xlsx(non_zero_table, "non_zero_expression.xlsx")

```



```{r, fig.width=15, fig.height=8}
## write all that into a function 

# Initialize a data frame to store the results
non_zero_table <- data.frame(
  Dataset = character(), 
  CellType = character(), 
  NonZeroCount = numeric(), 
  stringsAsFactors = FALSE
)

count_non_zero_expressions <- function(dataset, dataset_name, cell_type) {
  subset_seurat_object <- subset(dataset, subset = Class == cell_type)
  
  gene_expression <- FetchData(subset_seurat_object, vars = "Ntrk2")
  
  num_non_zero <- sum(gene_expression[["Ntrk2"]] > 0)
  
  # Add the result to the table
  non_zero_table <<- rbind(
    non_zero_table, 
    data.frame(
      Dataset = dataset_name, 
      CellType = cell_type, 
      NonZeroCount = num_non_zero, 
      stringsAsFactors = FALSE
    )
  )
}

# Define datasets, names, and cell types
datasets <- list(
  mouse_cortex_data = mouse_cortex_data, 
  eleven_samples_seurat = eleven_samples_seurat, 
  three_samples_seurat = three_samples_seurat
)

cell_types_of_interest <- c("Neurons", "Astrocytes", "Oligos")

# Loop through each dataset and cell type to count non-zero expressions
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Total non-zero expressions for the entire dataset
  ntrk2_expression <- FetchData(dataset, vars = "Ntrk2")
  num_non_zero <- sum(ntrk2_expression$Ntrk2 > 0)
  non_zero_table <- rbind(
    non_zero_table, 
    data.frame(
      Dataset = dataset_name, 
      CellType = "Total", 
      NonZeroCount = num_non_zero, 
      stringsAsFactors = FALSE
    )
  )
  
  for (cell_type in cell_types_of_interest) {
    count_non_zero_expressions(dataset, dataset_name, cell_type)
  }
}

# Calculate percentage for each cell type within each dataset
non_zero_table <- non_zero_table %>%
  group_by(Dataset) %>%
  mutate(Percentage = (NonZeroCount / NonZeroCount[CellType == "Total"]) * 100)

# Convert Dataset and CellType to factors for better plotting
non_zero_table$Dataset <- factor(non_zero_table$Dataset, levels = c("mouse_cortex_data", "eleven_samples_seurat", "three_samples_seurat"))
non_zero_table$CellType <- factor(non_zero_table$CellType, levels = c("Total", "Neurons", "Astrocytes", "Oligos"))

# Print the resulting table
print(non_zero_table)

# Plot Non-Zero Count as a line plot
plot_non_zero <- ggplot(non_zero_table, aes(x = CellType, y = NonZeroCount, color = Dataset, group = Dataset)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Non-Zero Count of Ntrk2 Expression", x = "Cell Type", y = "Non-Zero Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the Non-Zero Count line plot
print(plot_non_zero)

# Plot Percentage as a line plot
plot_percentage <- ggplot(non_zero_table, aes(x = CellType, y = Percentage, color = Dataset, group = Dataset)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Percentage of Non-Zero Ntrk2 Expression", x = "Cell Type", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the Percentage line plot
print(plot_percentage)


```


```{r, Ntrk2 expression, fig.height= 8, fig.width= 15}
featureplot_all <- FeaturePlot(mouse_cortex_data, features = "Ntrk2", cols = c("lightgrey", "blue"), raster=FALSE)
                    #scale_color_gradientn(colors = c("lightgrey", "skyblue", "blue"),
                        # breaks = c(0, 0.5, 1, 2, 3, 4),
                        # labels = c("0", "0.5", "1", "2", "3", "4"),
                        # limits = c(0, 4),
                        # oob = scales::squish,
                        #name = "Expression Level" 


featureplot_11 <- FeaturePlot(samples_seurat, features = "Ntrk2", cols = c("lightgrey", "blue"), raster = FALSE)
                    

VlnPlot(samples_seurat, features = "Ntrk2", group.by = "Class")

featureplot_3 <- FeaturePlot(three_samples_seurat, features = "Ntrk2", cols = c("lightgrey", "blue"), raster = FALSE)
                    

vln_3 <- VlnPlot(three_samples_seurat, features = "Ntrk2", group.by = "Class")
vln_11 <- VlnPlot(eleven_samples_seurat, features = "Ntrk2", group.by = "Class")
vln_all <- VlnPlot(mouse_cortex_data, features = "Ntrk2", group.by = "Class")
Dimplot_ <- DimPlot(mouse_cortex_data, group.by = "Class", label = TRUE, repel = TRUE, raster = FALSE)

featureplot_all
featureplot_11
featureplot_3
Dimplot_

vln_3
vln_11
vln_all
```

```{r, Scaling cell count}

# Function to scale the cell counts
scale_cell_counts <- function(cell_counts, neuron_scale, astrocyte_scale, oligo_scale) {
  cell_counts %>%
    dplyr::mutate(
      ScaledCount = dplyr::case_when(
        Class == "Neurons" ~ CellCount * neuron_scale,
        Class == "Astrocytes" ~ CellCount * astrocyte_scale,
        Class == "Oligos" ~ CellCount * oligo_scale,
        TRUE ~ CellCount
      )
    )
}

# Main filtering function with scaling included
filter_cells_by_min_count <- function(data, min_astrocyte_count, min_neuron_count, min_oligo_count) {
  # Summarize cell counts by SampleID and Class
  cell_counts <- data@meta.data %>%
    dplyr::group_by(SampleID, Class) %>%
    dplyr::summarise(CellCount = n(), .groups = 'drop') %>%
    as.data.frame()
  
  # Calculate scaling factors (example values, adjust based on actual proportions)
  neuron_scale <- 1
  astrocyte_scale <- 0.6
  oligo_scale <- 0.4
  
  # Apply scaling to the cell counts
  scaled_cell_counts <- scale_cell_counts(cell_counts, neuron_scale, astrocyte_scale, oligo_scale)
  
  # Define cell-type-specific minimum counts after scaling
  min_counts <- list(Astrocytes = min_astrocyte_count, Neurons = min_neuron_count, Oligos = min_oligo_count)
  
  # Initialize a list to store SampleIDs that meet the criteria for each cell type
  valid_sample_ids <- list()
  
  # Check each cell type with scaled counts
  for (cell_type in names(min_counts)) {
    ids <- scaled_cell_counts %>%
      dplyr::filter(Class == cell_type & ScaledCount >= min_counts[[cell_type]]) %>%
      dplyr::select(SampleID) %>%
      dplyr::distinct() %>%
      pull(SampleID)
    
    if (length(ids) == 0) {
      cat(sprintf("No samples have the defined minimum count of %d for %s.\n", min_counts[[cell_type]], cell_type))
    } else {
      valid_sample_ids[[cell_type]] <- ids
    }
  }
  
  # Check if there are any SampleIDs that meet the criteria for all specified cell types
  if (length(valid_sample_ids) < length(min_counts)) {
    cat("Some cell types did not meet the minimum cell count criteria.\n")
    return(NULL)
  }
  
  # Find common SampleIDs across all cell types
  samples_to_keep <- Reduce(intersect, valid_sample_ids)
  
  # Check if any samples meet the criteria
  if (length(samples_to_keep) == 0) {
    cat("No samples meet the minimum cell count across all specified cell types.\n")
    return(NULL)
  }
  
  # Subset the Seurat object to include only cells from the samples that meet the criteria
  seurat_subset <- subset(data, subset = SampleID %in% samples_to_keep)
  
  # Verify if subset is empty
  if (seurat_subset@meta.data %>% nrow() == 0) {
    cat("The resulting subset of the data has no cells. Check the filtering criteria.\n")
    return(NULL)
  }
  
  # Cell Type per sample
  target_cell_counts_by_sample <- table(seurat_subset@meta.data$SampleID, seurat_subset@meta.data$Class)
  print(target_cell_counts_by_sample)
  
  # Convert the table to a data frame for plotting
  target_counts_df <- as.data.frame(target_cell_counts_by_sample)
  names(target_counts_df) <- c("Sample", "CellType", "Count")
  target_counts_df <- target_counts_df %>%
    dplyr::filter(CellType %in% c("Astrocytes", "Oligos", "Neurons"))
  
  return(target_counts_df)
}



```






