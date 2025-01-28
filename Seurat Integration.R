# Load required libraries
library(Seurat)
library(patchwork)
library(ggplot2)

# Load individual Seurat objects into the environment
E12.5_1 <- readRDS(".../E12.5_1.RDS")
E12.5_2 <- readRDS(".../E12.5_2.RDS")
E12.5_3 <- readRDS(".../E12.5_3.RDS")
E13.5_1 <- readRDS(".../E13.5_1.RDS")
E13.5_2 <- readRDS(".../E13.5_2.RDS")
E13.5_3 <- readRDS(".../E13.5_3.RDS")
E14.5_1 <- readRDS(".../E14.5_1.RDS")
E14.5_2 <- readRDS(".../E14.5_2.RDS")
E14.5_3 <- readRDS(".../E14.5_3.RDS")
E15.5_1 <- readRDS(".../E15.5_1.RDS")
E15.5_2 <- readRDS(".../E15.5_2.RDS")
E15.5_3 <- readRDS(".../E15.5_3.RDS")
E18.5_1 <- readRDS(".../E18.5_1.RDS")
E18.5_2 <- readRDS(".../E18.5_2.RDS")
E18.5_3 <- readRDS(".../E18.5_3.RDS")


# Assign individual batch information to all datasets
E12.5_1$Batch <- "E12.5_1"
E12.5_2$Batch <- "E12.5_2"
E12.5_3$Batch <- "E12.5_3"
E13.5_1$Batch <- "E13.5_1"
E13.5_2$Batch <- "E13.5_2"
E13.5_3$Batch <- "E13.5_3"
E14.5_1$Batch <- "E14.5_1"
E14.5_2$Batch <- "E14.5_2"
E14.5_3$Batch <- "E14.5_3"
E15.5_1$Batch <- "E15.5_1"
E15.5_2$Batch <- "E15.5_2"
E15.5_3$Batch <- "E15.5_3"
E18.5_1$Batch <- "E18.5_1"
E18.5_2$Batch <- "E18.5_2"
E18.5_3$Batch <- "E18.5_3"

# Assign grouped batch labels

E12.5_1$GroupBatch <- "E12.5"
E12.5_2$GroupBatch <- "E12.5"
E12.5_3$GroupBatch <- "E12.5"
E13.5_1$GroupBatch <- "E13.5"
E13.5_2$GroupBatch <- "E13.5"
E13.5_3$GroupBatch <- "E13.5"
E14.5_1$GroupBatch <- "E14.5"
E14.5_2$GroupBatch <- "E14.5"
E14.5_3$GroupBatch <- "E14.5"
E15.5_1$GroupBatch <- "E15.5"
E15.5_2$GroupBatch <- "E15.5"
E15.5_3$GroupBatch <- "E15.5"
E18.5_1$GroupBatch <- "E18.5"
E18.5_2$GroupBatch <- "E18.5"
E18.5_3$GroupBatch <- "E18.5"

# Merge all datasets into one object
obj <- merge(
  E12.5_1, 
  y = list(
    E12.5_2, E12.5_3, 
    E13.5_1, E13.5_2, E13.5_3, 
    E14.5_1, E14.5_2, E14.5_3, 
    E15.5_1, E15.5_2, E15.5_3, 
    E18.5_1,E18.5_2, E18.5_3
  ),
  add.cell.ids = c(
    "E12.5_1","E12.5_2", "E12.5_3", 
    "E13.5_1", "E13.5_2", "E13.5_3", 
    "E14.5_1", "E14.5_2", "E14.5_3", 
    "E15.5_1", "E15.5_2", "E15.5_3", 
    "E18.5_1","E18.5_2", "E18.5_3"
  )
)


# Preprocess data (normalize, find variable features, scale)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Perform integrations with multiple methods (e.g., CCA Integration)
obj <- IntegrateLayers(
  object = obj, method =  CCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# Perform clustering and visualization for each method
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")


# Visualization: UMAP with batch, group-batch, and clustering information
# Create UMAP plots for multiple grouping variables in one step
p <- DimPlot(object = obj,reduction = "umap.cca",group.by = c("Batch", "GroupBatch", "cca_clusters"),label = TRUE,            
  label.size = 2)
p

# Single plot, grouped by GroupBatch
p1 <- DimPlot(object = obj, reduction = "umap.cca", group.by = "GroupBatch")
p1

# Split plot by GroupBatch, one plot per group
p2 <- DimPlot(object = obj, reduction = "umap.cca", split.by = "GroupBatch")
p2

saveRDS(obj, file = "~/integrated.rds")

# Join layers
joined_obj <- JoinLayers(obj)
saveRDS(joined_obj, file = "joined_obj.rds")

# Split layers
#split_obj <- SplitLayers(joined_obj)
#saveRDS(split_obj, file = "split_obj.rds")

#Differntial Gene Expression analysis using FindAllMarkers
Integrated.markers <- FindAllMarkers(joined_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- Integrated.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(joined_obj, features = top10$gene, label = TRUE, size = 3) + NoLegend()

# Identify markers for a specific cluster compared to all other cells
cluster_0_markers <- FindMarkers(object = joined_obj,ident.1 = 0, ident.2 = NULL,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)

# View the top markers
head(cluster_0_markers, n = 20)

# Evaluate gene exprssion patterns of the Differntial Gene Expression analysis
FeaturePlot(object = obj, features = c("Sp7","Msx1","Barx1",...), label=T, order=T)
DotPlot(obj, features = c("Sp7","Msx1","Barx1",...) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
# Define new cluster identities
new.cluster.ids <- c(
"Mesenchyme",         # Cluster 0
"Mesenchyme",         # Cluster 1
"Mesenchyme",         # Cluster 2
"Mesenchyme",         # Cluster 3
"Mesenchyme",         # Cluster 4
"Mesenchyme",         # Cluster 5
"Mesenchyme",         # Cluster 6
"Mesenchyme",         # Cluster 7
"Epithelial cells",   # Cluster 8
"Epithelial cells",   # Cluster 9
"Mesenchyme",         # Cluster 10
"Endothelial cells",  # Cluster 11
"Epithelial cells",   # Cluster 12
"Erythroid cells",    # Cluster 13
"Myeloid cells",      # Cluster 14
"Mesenchyme",         # Cluster 15
"Glial cells",        # Cluster 16
"Neuronal cells",     # Cluster 17
"Myogenic cells",     # Cluster 18
"Myeloid cells",      # Cluster 19
"Erythroid cells",    # Cluster 20
"Mesenchyme",         # Cluster 21
"Mesenchyme",         # Cluster 22
"Mesenchyme",         # Cluster 23
"Epithelial cells",   # Cluster 24
"Myogenic cells"      # Cluster 25
)
        
# Map new names to the existing cluster levels
names(new.cluster.ids) <- levels(obj)  # Ensure this matches your Seurat object levels
RenamedObj <- RenameIdents(obj, new.cluster.ids)
        
# Verify the updated cluster identities
Idents(RenamedObj)
        
DimPlot(RenamedObj, reduction = "umap.cca", label = TRUE, label.size = 6) + 
theme(legend.text = element_text(size = 10))  # Customize legend font size
        
subset_obj <- subset(obj, idents=c(0, 1, 2, 3, 4, 5, 6, 7, 10, 15, 21, 22, 23))

p3 <- DimPlot(subset_obj, 
              reduction = "umap.cca", 
              group.by = c("Batch", "GroupBatch", "cca_clusters"), 
              label.size = 2)  # Label size for clusters
p3

saveRDS(subset_obj, file = ".../subset.rds")

# Extract metadata
metadata <- seurat_obj@meta.data

# Summarize cell counts by GroupBatch and clusters
cell_counts <- metadata %>%
  group_by(GroupBatch = metadata$grouping_variable, 
           renamed_clusters = metadata$cluster_variable) %>%
  summarize(cell_count = n(), .groups = "drop")

# Calculate percentages
cell_counts <- cell_counts %>%
  group_by(GroupBatch) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100)

# Plot
ggplot(cell_counts, aes(x = GroupBatch, y = percent, fill = renamed_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  labs(
    x = "GroupBatch",
    y = "Percentage",
    fill = "Renamed Clusters"
  )
