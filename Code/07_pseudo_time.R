suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ArchR)
  library(tidyverse)
  library(monocle3)
  library(SeuratWrappers)
})

# Step 1: Load RNA and ATAC reference objects --------------------------------------
rna_path <- "/home/Data_Drive_8TB_4/ghlqls/Thyroid_scRNA/thyroid.all.102845cell_ReAnno.rds"
atac_path <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC/Archr_Result/Archr_TotalOutput_RemoveTwo_Global_240304_Epi/Signac_Obj.rds"

thyroid_rna <- readRDS(rna_path)
thyroid_atac <- readRDS(atac_path)
DefaultAssay(thyroid_atac) <- "ACTIVITY"

# Step 2: Build transfer anchors between modalities -------------------------------
thyroid_rna <- FindVariableFeatures(thyroid_rna, nfeatures = 2000)
thyroid_atac <- FindVariableFeatures(thyroid_atac, nfeatures = 2000)

integration_dir <- "/home/Data_Drive_8TB_4/ghlqls/Seurat_Obj/Integration"
dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)

anchors <- FindTransferAnchors(
  reference = thyroid_atac,
  query = thyroid_rna,
  features = VariableFeatures(thyroid_rna),
  reference.assay = "ACTIVITY",
  query.assay = "RNA",
  reduction = "cca"
)
saveRDS(anchors, file = file.path(integration_dir, "transfer_anchors_epi.rds"))

celltype_predictions <- TransferData(
  anchorset = anchors,
  refdata = thyroid_atac$sub.annotation,
  weight.reduction = thyroid_rna[["harmony"]],
  dims = 2:30
)
write_csv(as.data.frame(celltype_predictions), file.path(integration_dir, "celltype_predictions_epi.csv"))

# Step 3: Co-embed RNA and ATAC profiles ------------------------------------------
genes_use <- intersect(VariableFeatures(thyroid_rna), VariableFeatures(thyroid_atac))
rna_reference <- GetAssayData(thyroid_rna, assay = "RNA", slot = "data")[genes_use, ]
imputed_rna <- TransferData(
  anchorset = anchors,
  refdata = rna_reference,
  weight.reduction = thyroid_atac[["IterativeLSI"]],
  dims = 2:30
)
thyroid_atac[["RNA"]] <- imputed_rna

coembed <- merge(thyroid_rna, thyroid_atac)
coembed <- ScaleData(coembed, features = genes_use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes_use)
coembed <- RunUMAP(coembed, dims = 1:30)

coembed$dataset_label <- ifelse(coembed$orig.ident == "RNA", coembed$dataset, coembed$Sample)

# Step 4: Add module scores for key gene programs ----------------------------------
project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))
pseudotime_dir <- file.path(epi_project_dir, "Pseudotime")
dir.create(pseudotime_dir, recursive = TRUE, showWarnings = FALSE)
score_geneset <- read_csv(file.path(project_root, "RNAseq_Analysis_Data", "score_geneset.csv"))
module_sets <- list(
  TDS = na.omit(score_geneset$TDS_score),
  ERK = na.omit(score_geneset$ERK_score),
  BRS = na.omit(score_geneset$BRS_score)
)

coembed <- AddModuleScore(coembed, features = list(module_sets$TDS), name = "TDS_score")
coembed <- AddModuleScore(coembed, features = list(module_sets$ERK), name = "ERK_score")
coembed <- AddModuleScore(coembed, features = list(module_sets$BRS), name = "BRS_score")

# Step 5: Convert to Monocle object and compute pseudotime -------------------------
cds <- as.cell_data_set(coembed)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)
root_cell <- colnames(coembed)[which.max(coembed$TDS_score1)]
cds <- order_cells(cds, root_cells = root_cell)

# Step 6: Plot pseudotime results --------------------------------------------------
pseudo_min <- min(pseudotime(cds))
pseudo_max <- max(pseudotime(cds))
pseudo_plot <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  cell_size = 1,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  theme_void() +
  scale_color_gradientn(
    colors = c("#000436", "#0B27CF", "#6F33FB", "#FD6E8F", "#FFF95A"),
    limits = c(pseudo_min, pseudo_max)
  ) +
  labs(title = "Pseudo-time", color = "pseudotime") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "bottom"
  )

pdf(file.path(pseudotime_dir, "UMAP_pseudotime.pdf"), width = 6, height = 6)
print(pseudo_plot)
dev.off()

saveRDS(cds, file.path(pseudotime_dir, "cds_coembed.rds"))
