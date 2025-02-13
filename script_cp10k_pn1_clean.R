################################################################################
############### R Script to Analyze scRNA seq from PN1 mice limbs ##############
################################################################################
#Article: https://doi.org/10.1111/imcb.12718

#Remove all:
rm(list = ls())

#Set workind directory:
setwd("C:/Users/jmper/OneDrive/Escritorio/Atlas/sc_mouse")

#Load libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)


#Read data original data and transform it into a Seurat object:
data.orig <- readRDS("E18_PN1_counts.rds") # Data has only counts.

seur.obj <- CreateSeuratObject(counts = data.orig)

#Eliminate rds matrix for more space:
rm(data.orig)


#Mitochondrial percent genes:
seur.obj[["percent.mt"]] <- PercentageFeatureSet(seur.obj, pattern = "^mt-")
View(seur.obj@meta.data)

#Check quality of the data:
qc.vln1 <- VlnPlot(seur.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, layer = "counts")
qc.sc1 <- FeatureScatter(seur.obj, feature1= "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = "lm")

#Filter low quality data:
filtered.obj <- subset(seur.obj, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 10) #According to the original paper.
qc.sc2 <- FeatureScatter(filtered.obj, feature1= "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = "lm")
qc.vln2 <- VlnPlot(filtered.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, layer = "counts")

#Normalize to CP10K:
filtered.obj <- NormalizeData(filtered.obj, normalization.method = "RC", scale.factor = 10000)

#Once normalized, we add the cell type from meta data to the data since they have been kindly provided by our collaborators:
metadata <- read.csv("E18_PN1_metadata.csv", header = TRUE, row.names = 1)
cell_type <- metadata[,6]
filtered.obj$cell_type <- cell_type
View(filtered.obj@meta.data) #Cell types correctly added.

#Now, we are going to split the data between E18 
seurat.pn1 <- subset(filtered.obj, subset = orig.ident == "PN1")
seurat.e18 <- subset(filtered.obj, subset = orig.ident == "E18")

#We are going to change the cell names (barcodes) in the matrix for the cell types as cell names.
#For PN1:
View(seurat.pn1@meta.data)
head(seurat.pn1@meta.data$cell_type) #Check we have the cell types.
head(colnames(seurat.pn1@assays$RNA$data))  # Check first 10 column names in matrix
head(rownames(seurat.pn1@meta.data))        # Check first 10 row names in metadata
rownames(seurat.pn1@meta.data) <- colnames(seurat.pn1) #reset original names (barcodes)
all(rownames(seurat.pn1@meta.data) == colnames(seurat.pn1@assays$RNA$data)) # Check we do not miss any cell #Result is TRUE
colnames(seurat.pn1@assays$RNA) <- seurat.pn1$cell_type #Done changing names

#For E18:
View(seurat.e18@meta.data)
head(seurat.e18@meta.data$cell_type) #Check we have the cell types.
head(colnames(seurat.e18@assays$RNA$data))  # Check first 10 column names in matrix
head(rownames(seurat.e18@meta.data))        # Check first 10 row names in metadata
rownames(seurat.e18@meta.data) <- colnames(seurat.e18) #reset original names (barcodes)
all(rownames(seurat.e18@meta.data) == colnames(seurat.e18@assays$RNA$data)) # Check we do not miss any cell #Result is TRUE
colnames(seurat.e18@assays$RNA) <- seurat.e18$cell_type #Done changing names

#Selecting columns of interest:
pn1.matrix <- as.matrix(seurat.pn1@assays$RNA$data)
pn1.df <- as.data.frame(pn1.matrix)
rm(pn1.matrix)
pn1.int <- select(prueba1, contains(c("Oc","Myo"))) #Select cell types of interest.
any(is.na(pn1.int)) #Check conversion went well and there is no NA

#Stimating replicates:
# Get the column names
col_names <- colnames(pn1.int)

# Identify unique patterns (e.g., B_cells, Basophils, Endo) in column names
# Assuming patterns are unique names without trailing numbers
patterns <- unique(gsub(".\\d+", "", col_names))  # removes trailing numbers after underscore if present

# Initialize an empty list to store the resulting mean columns
mean_columns <- list()

for (pattern in patterns) {
  # Find columns matching the current pattern
  matching_cols <- grep(paste0("^", pattern), col_names, value = TRUE)
  
  # Check if there are columns that match the pattern
  if (length(matching_cols) > 0) {
    # Calculate group size to obtain exactly 5 mean columns for each pattern
    group_size <- ceiling(length(matching_cols) / 10)
    
    # Split the columns into 5 groups (or as close as possible) by calculated group size
    groupings <- split(matching_cols, ceiling(seq_along(matching_cols) / group_size))
    
    # Calculate mean for each group and add to mean_columns list
    for (i in seq_len(min(10, length(groupings)))) {  # limit to 5 groups
      group <- groupings[[i]]
      mean_columns[[paste0(pattern, "_mean_", i)]] <- rowMeans(pn1.int[group], na.rm = TRUE)
    }
  }
}

# Combine mean columns into a new data frame
pn1_mean_df <- as.data.frame(mean_columns)

pn1_mean_names <- colnames(pn1_mean_df)
pn1_mean_names <- gsub("_mean_[0-9]*$","",pn1_mean_names)
colnames(pn1_mean_df) <- pn1_mean_names

#Reorder columns:
pn1_mean_ord_df <-  pn1_mean_df[, order(names(pn1_mean_df))]

#Make some final changes in names and digits:
pn1_final_names <- colnames(pn1_mean_ord_df)
pn1_final_names <- gsub("Oc", "Osteoclasts", pn1_final_names)
pn1_final_names <- gsub("Myo", "Muscle_cells", pn1_final_names)
pn1_final_names <- gsub("\\.[0-9]*$", "", pn1_final_names)
colnames(pn1_mean_ord_df) <- pn1_final_names

#Round to 2 decimal numbers:
pn1_mean_ord_df <- round(pn1_mean_ord_df, digits=2)

#Save results:
write.table(pn1_mean_ord_df, "pn1_data_atlas_CP10K.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
