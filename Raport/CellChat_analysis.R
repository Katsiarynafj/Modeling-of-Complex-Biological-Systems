library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(Seurat)
library(tidyverse)  
library(patchwork)
library(dplyr)
library(SingleR)
library(celldex)

options(future.globals.maxSize = 15 * 1024^3)
future::plan("multisession", workers = 1)

ref = MouseRNAseqData()

parent_dir1 = "C:/Users/Kasia/Downloads/modeling/atlas_data/atlas"
parent_dir2 =  "C:/Users/Kasia/Downloads/modeling"

seurat_atlas = readRDS(file.path(parent_dir1, "atlas_seurat_E75.rds"))
seurat_vasa = readRDS(file.path(parent_dir2, "VASA_seq_E7.5.rds"))

colnames(seurat_atlas@meta.data)
colnames(seurat_vasa@meta.data)

#removing all unnecessary column from seurat_atlas
keep_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample", "celltype")
seurat_atlas@meta.data <- seurat_atlas@meta.data[, keep_cols]
seurat_vasa$celltype <- "unknown"

rownames(seurat_vasa)[is.na(rownames(seurat_vasa))] <- paste0("UnknownGene_", seq_len(sum(is.na(rownames(seurat_vasa)))))
rownames(seurat_vasa) <- make.unique(rownames(seurat_vasa))

DefaultAssay(seurat_vasa) <- "RNA"
DefaultAssay(seurat_atlas) <- "RNA"

seurat_vasa <- NormalizeData(seurat_vasa, verbose = FALSE)
seurat_atlas <- NormalizeData(seurat_atlas, verbose = FALSE)

seurat_vasa <- FindVariableFeatures(seurat_vasa, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat_atlas <- FindVariableFeatures(seurat_atlas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = list(seurat_vasa, seurat_atlas), nfeatures = 2000)


seurat_vasa <- ScaleData(seurat_vasa, features = features, verbose = FALSE)
seurat_atlas <- ScaleData(seurat_atlas, features = features, verbose = FALSE)


target_cells_atlas <- 7000 #subseting the seurat_atlas data
if (ncol(seurat_atlas) > target_cells_atlas) {
  message(paste("Podpróbkowanie seurat_atlas z", ncol(seurat_atlas), "do", target_cells_atlas, "komórek."))
  set.seed(123) # Dla powtarzalności wyników
  seurat_atlas_subsampled <- subset(seurat_atlas, cells = sample(colnames(seurat_atlas), target_cells_atlas))
} else {
  message("seurat_atlas jest już mniejszy lub równy docelowej liczbie komórek. Brak podpróbkowania.")
  seurat_atlas_subsampled <- seurat_atlas
}

seurat_vasa <- JoinLayers(seurat_vasa)

seurat_vasa <- RunPCA(seurat_vasa, features = features, verbose = FALSE, npcs = 50) 
seurat_atlas_subsampled <- RunPCA(seurat_atlas_subsampled, features = features, verbose = FALSE, npcs = 50) 

anchors <- FindIntegrationAnchors(object.list = list(seurat_vasa, seurat_atlas_subsampled), 
                                  anchor.features = features,
                                  reduction = "rpca") 

combined_seurat <- IntegrateData(anchorset = anchors)

combined_seurat

DefaultAssay(combined_seurat) <- "integrated"


combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)

combined_seurat <- RunPCA(combined_seurat, verbose = FALSE)
ElbowPlot(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:15, verbose = FALSE)
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:15, verbose = FALSE)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5, verbose = FALSE)
DimPlot(combined_seurat, reduction = "umap", label = TRUE)



singler_predictions <- SingleR(test = GetAssayData(combined_seurat, assay = "integrated", layer = "data"), 
                               ref = ref,
                               labels = ref$label.fine, 
                               BPPARAM = BiocParallel::bpparam() 
)

combined_seurat$SingleR_cell_type <- singler_predictions$pruned.labels 

DimPlot(combined_seurat, reduction = "umap", label = TRUE, group.by = "seurat_clusters", repel = TRUE) 

DimPlot(combined_seurat, reduction = "umap", label = TRUE, group.by = "SingleR_cell_type", repel = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 2))) # Zwiększ rozmiar punktów w legendzie

saveRDS(combined_seurat, file = "C:/Users/Kasia/Downloads/modeling/combined_data.rds")

combined_seurat <- readRDS(file.path("C:/Users/Kasia/Downloads/modeling/combined_data.rds"))


#-----------CellCHAT Analysis---------------------------------------


#I had to remove 13_cluster, cause it's not represented in altas object, which results in having different amount of layers, so most of the 
#anlaysis on merge CellChat object wont work, even if u force carring the empty layer 
vasa_seurat_subset <- subset(combined_seurat, subset = orig.ident == "SeuratProject")
vasa_seurat_subset <- subset(vasa_seurat_subset, subset = seurat_clusters != "13")
vasa_seurat_subset$seurat_clusters <- droplevels(vasa_seurat_subset$seurat_clusters)
vasa_seurat_subset$seurat_clusters <- factor(paste0("Cluster_", as.character(vasa_seurat_subset$seurat_clusters)))

atlas_seurat_subset <- subset(combined_seurat, subset = orig.ident == "cell") 
atlas_seurat_subset$seurat_clusters <- droplevels(atlas_seurat_subset$seurat_clusters)
atlas_seurat_subset$seurat_clusters <- factor(paste0("Cluster_", as.character(atlas_seurat_subset$seurat_clusters)))


cellchat_vasa <- createCellChat(object = vasa_seurat_subset,
                                group.by = "seurat_clusters", 
                                assay = "RNA")
# For atlas data
cellchat_atlas <- createCellChat(object = atlas_seurat_subset,
                                 group.by = "seurat_clusters", 
                                 assay = "RNA")

cellchat_atlas@idents <- factor(cellchat_atlas@idents, levels = all_cluster_levels)
cellchat_vasa@idents <- factor(cellchat_vasa@idents, levels = all_cluster_levels)


print(table(cellchat_atlas@idents))
print(table(cellchat_vasa@idents))


CellChatDB.use = CellChatDB.mouse
#error fixing
CellChatDB.use[["interaction"]] <- subset(CellChatDB.use[["interaction"]],
                                          !(ligand %in% c("H2-BI", "H2-Ea-ps")))

cellchat_vasa@DB <- CellChatDB.use
cellchat_atlas@DB <- CellChatDB.use


# For vasa
cellchat_vasa <- subsetData(cellchat_vasa) 
cellchat_vasa <- identifyOverExpressedGenes(cellchat_vasa)
cellchat_vasa <- identifyOverExpressedInteractions(cellchat_vasa)
cellchat_vasa <- projectData(cellchat_vasa, PPI.mouse) 

# For atlas
cellchat_atlas <- subsetData(cellchat_atlas)
cellchat_atlas <- identifyOverExpressedGenes(cellchat_atlas)
cellchat_atlas <- identifyOverExpressedInteractions(cellchat_atlas)
cellchat_atlas <- projectData(cellchat_atlas, PPI.mouse)

print(table(cellchat_atlas@idents))
print(table(cellchat_vasa@idents))


# For vasa

options(future.globals.maxSize = 15 * 1024^3)
future::plan("multisession", workers = 1)

cellchat_vasa <- computeCommunProb(cellchat_vasa, type = "triMean") 
cellchat_vasa <- filterCommunication(cellchat_vasa, min.cells = 3) 
cellchat_vasa <- computeCommunProbPathway(cellchat_vasa)
cellchat_vasa <- aggregateNet(cellchat_vasa) 
cellchat_vasa <- netAnalysis_computeCentrality(cellchat_vasa, slot.name = "netP") 

options(future.globals.maxSize = 15 * 1024^3)
future::plan("multisession", workers = 1)

# For atlas
cellchat_atlas <- computeCommunProb(cellchat_atlas, type = "triMean")
cellchat_atlas <- filterCommunication(cellchat_atlas, min.cells = 3)
cellchat_atlas <- computeCommunProbPathway(cellchat_atlas)
cellchat_atlas <- aggregateNet(cellchat_atlas)
cellchat_atlas <- netAnalysis_computeCentrality(cellchat_atlas, slot.name = "netP")

print(table(cellchat_atlas@idents))
print(table(cellchat_vasa@idents))
unique(cellchat_atlas@idents)

object.list <- list(Vasa = cellchat_vasa, Atlas = cellchat_atlas)
cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))

cellchat_merge

saveRDS(cellchat_merge, file = "C:/Users/Kasia/Downloads/modeling/cellchat_merge.rds")

cellchat_merge <- readRDS(file.path("C:/Users/Kasia/Downloads/modeling/cellchat_merge.rds"))


library(ggplot2)
library(patchwork)

# A global comparison of the number of interactions and communication strength (figure 2)
gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "weight")
gg1_enhanced <- gg1 + theme(axis.text.x = element_text(size = 12), 
                            axis.text.y = element_text(size = 12)) 
gg2_enhanced <- gg2 + theme(axis.text.x = element_text(size = 12),
                            axis.text.y = element_text(size = 12)) 

gg1_enhanced + gg2_enhanced


#Circle plots (figure 3)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merge, weight.scale = T)
netVisual_diffInteraction(cellchat_merge, weight.scale = T, measure = "weight")


#Signaling Role Scatter Plots(figure 4)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

#Pathway-Specific Changes in Cell Signaling (figure 5)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "Cluster_12")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = "Cluster_15")
patchwork::wrap_plots(plots = list(gg1,gg2))

# Ranking of Differential Signaling Pathways (figure 6)   
rankNet(cellchat_merge, mode = "comparison", comparison = c(1,2))