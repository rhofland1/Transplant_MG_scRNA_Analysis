library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)

#Read Filtered Matrix
data <- Read10X(data.dir = "D:/ChIP seq and RNA seq analysis/Arlotta_lab scRNA/E17.5.filtered_feature_bc_matrix")
data <- CreatemicroObject(counts = data, project='E17.5', min.cells = 3, min.features = 200)

#Filter Data
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA >200 & nFeature_RNA <2500 & percent.mt <5)



#Process Data
micro<- NormalizeData(micro)
micro<- FindVariableFeatures(micro, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(micro), 10)
top10
plot1 <- VariableFeaturePlot(micro)
LabelPoints(plot = plot1, points = top10, repel = T)
all.genes <- rownames(micro)
micro<- ScaleData(micro, features = all.genes)
micro<- RunPCA(micro, features = VariableFeatures(object = micro))
DimHeatmap(micro, dims = 1:6, cells = 500, balanced = T)
ElbowPlot(micro)
micro<- FindNeighbors(micro, dims = 1:30)
micro<- FindClusters(micro, resolution = c(0.3, 0.5, 0.7, 1), graph.name = "RNA_snn")
head(micro@meta.data)
Idents(micro) <- "RNA_snn_res.0.7"
micro<- RunUMAP(micro, dims = 1:30)
#Pre-load data

micro<- FindClusters(micro, resolution =0.7)
micro<- RunUMAP(micro, dims = 1:30)
DimPlot(micro, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(micro, features = c("Apoe", "Ccr1", "Ms4a6c", "Ms4a7", "Tmem176a", "Ms4a6b", "Gpx3", "Tmem176b", "Ms4a6d"))
#pseudotime
cds <- as.cell_data_set(micro)
head(colData(cds))
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
head(counts(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- micro@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- micro@reductions$umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
micro.pseudo <- as.data.frame(colData(cds))

ggplot(micro.pseudo, aes(monocle3_pseudotime, micro$column, fill = micro$column)) + geom_boxplot()
ggplot(micro.pseudo, aes(monocle3_pseudotime, reorder(micro$column, monocle3_pseudotime), fill = micro$column)) + geom_boxplot()
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
micro$pseudotime <- pseudotime(cds)
FeaturePlot(micro, features = "pseudotime")
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Sox2", "Eomes", "Neurod2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )

#Save Open Plots
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="D:/Rohan/Cluster_Data/20042025_RPlots_Heatmaps_Pseudotime")
plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
sorted.png.names <- gsub(plots.dir.path, "D:/Rohan/Cluster_Data/20042025_RPlots_Heatmaps_Pseudotime", row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0("D:/Rohan/Cluster_Data/20042025_RPlots_Heatmaps_Pseudotime", 1:length(sorted.png.names), ".png")
graphics.off()

DimPlot(micro, reduction = "umap", label = TRUE, repel = TRUE, group.by = "microlia_Type")
DimPlot(micro_and_immune, reduction = "umap", label = TRUE, repel = TRUE, group.by = "microlia_Type")
DimPlot(micro, reduction = "umap", label = TRUE, repel = TRUE, group.by = "microlia_Type")