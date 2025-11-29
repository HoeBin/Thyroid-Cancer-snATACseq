suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggpubr)
})

# Step 1: Configure directories and load ArchR project -----------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))
output_dir <- file.path(epi_project_dir, "Markers")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
addArchRThreads(threads = 16)

# Step 2: Prepare gene score matrix for downstream heatmaps -----------------------
gene_score <- getMatrixFromProject(proj_epi, useMatrix = "GeneScoreMatrix")
gene_scores <- as.matrix(assays(gene_score)$GeneScoreMatrix)
rownames(gene_scores) <- rowData(gene_score)$name

# Step 3: Load RNA DEG references --------------------------------------------------
deg_external <- read_csv(file.path(project_root, "RNAseq_Analysis_Data", "epi.deg.0.01.csv")) %>%
  mutate(cluster = str_replace_all(cluster, "_", "."))
deg_internal <- readRDS(file.path(project_root, "RNAseq_Analysis_Data", "Epi_sub_markers.rds")) %>%
  mutate(cluster = str_replace_all(cluster, "_", "."))

get_top_genes <- function(deg_tbl, top_n = 50, min_fc = 0.25) {
  deg_tbl %>%
    filter(avg_log2FC > min_fc) %>%
    group_by(cluster) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

top_external <- get_top_genes(deg_external, top_n = 20)
top_internal <- get_top_genes(deg_internal, top_n = 20)

# Step 4: Identify marker genes from ATAC gene scores -----------------------------
markers_sub <- getMarkerFeatures(
  ArchRProj = proj_epi,
  useMatrix = "GeneScoreMatrix",
  groupBy = "sub.annotation",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markers_sub, file.path(output_dir, "markersGS_sub_annotation.rds"))
marker_list <- getMarkers(markers_sub, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

# Step 5: Plot heatmaps for marker genes ------------------------------------------
heatmap_sub <- plotMarkerHeatmap(
  seMarker = markers_sub,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = FALSE,
  nLabel = 5
)
plotPDF(heatmap_sub, name = "GeneScore_MarkerHeatmap_SubAnno.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 10, height = 12)

# Step 6: Compare marker genes with RNA references ---------------------------------
overlap_summary <- function(marker_genes, deg_tbl) {
  expand_grid(marker_cluster = names(marker_genes), rna_cluster = unique(deg_tbl$cluster)) %>%
    mutate(
      overlap = map2_dbl(marker_cluster, rna_cluster, ~ length(intersect(marker_genes[[.x]]$name, deg_tbl$gene[deg_tbl$cluster == .y]))),
      marker_cluster = factor(marker_cluster, levels = names(marker_genes))
    )
}

overlap_ext <- overlap_summary(marker_list, top_external)
overlap_int <- overlap_summary(marker_list, top_internal)

p_external <- ggplot(overlap_ext, aes(x = marker_cluster, y = rna_cluster, fill = overlap)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Overlap with RNA DEG (External)", x = "ATAC cluster", y = "RNA cluster", fill = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_internal <- ggplot(overlap_int, aes(x = marker_cluster, y = rna_cluster, fill = overlap)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Overlap with RNA DEG (Internal)", x = "ATAC cluster", y = "RNA cluster", fill = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "RNA_ATAC_overlap_external.pdf"), p_external, width = 6, height = 4)
ggsave(file.path(output_dir, "RNA_ATAC_overlap_internal.pdf"), p_internal, width = 6, height = 4)

# Step 7: Build simple Venn diagram-ready lists -----------------------------------
selected_clusters <- names(marker_list)
marker_gene_sets <- map(selected_clusters, ~ head(marker_list[[.x]]$name, 100))
names(marker_gene_sets) <- selected_clusters
saveRDS(marker_gene_sets, file.path(output_dir, "marker_gene_sets.rds"))

# Step 8: Optionally compute cluster-level heatmaps -------------------------------
annotation_data <- as.data.frame(proj_epi@cellColData[, c("sub.annotation", "group", "Sample")])
sample_levels <- unique(annotation_data$Sample)
sample_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(sample_levels)), sample_levels)
annotation_colors <- list(
  sub.annotation = setNames(brewer.pal(6, "Set1"), c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high")),
  group = c(BRAF_PTC = "darkgreen", BRAF_ATC = "red", RAS_FTC = "orange", RAS_ATC = "purple", NBNR_PDTC = "grey50"),
  Sample = sample_colors
)

write_rds(
  list(gene_scores = gene_scores, marker_list = marker_list, annotation = annotation_data, colors = annotation_colors),
  file.path(output_dir, "gene_score_heatmap_inputs.rds")
)
