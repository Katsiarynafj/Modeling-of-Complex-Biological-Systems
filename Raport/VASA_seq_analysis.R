library(Seurat)

#geting the data
parent_dir <- "C:/Users/Kasia/Downloads/modeling/VASA_seq_data_E7.5"
save_dir =  "C:/Users/Kasia/Downloads/modeling"

file_paths <- list.files(path = parent_dir, pattern = "expression_matrix.tsv", recursive = TRUE, full.names = TRUE)
sample_names <- basename(dirname(file_paths))

#cause the 7.5E is representted by multiple sample i make the seuret object from each, and merging them after 
# Loading each file into a Seurat object, naming as sample ID
seurat_list <- lapply(seq_along(file_paths), function(i) {
  expr <- read.delim(file_paths[i], row.names = 1, check.names = FALSE)
  seu <- CreateSeuratObject(counts = expr)
  seu$sample <- sample_names[i]
  seu
})

dim()

# Merge all Seurat objects into one
combined_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)

dim(combined_seurat)

current_gene_names <- rownames(combined_seurat[["RNA"]])
# Extract gene symbols from gene names(to unite the format)
new_gene_names <- sapply(current_gene_names, function(x) {
  parts <- strsplit(x, split = "-")[[1]]
  return(parts[2])
})
if (any(duplicated(new_gene_names))) {
  new_gene_names <- make.unique(new_gene_names)
}
rownames(combined_seurat[["RNA"]]) <- new_gene_names

save_path <- file.path(save_dir, "VASA_seq_E7.5.rds")
saveRDS(combined_seurat, file = save_path)