suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
})

# Step 1: Configure directories and load data --------------------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))
output_dir <- file.path(epi_project_dir, "Motif")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
addArchRThreads(threads = 16)

# Step 2: Annotate motifs and TFBS collections -------------------------------------
proj_epi <- addMotifAnnotations(ArchRProj = proj_epi, motifSet = "cisbp", name = "Motif", force = TRUE)
proj_epi <- addArchRAnnotations(ArchRProj = proj_epi, collection = "EncodeTFBS", force = TRUE)

# Step 3: Load marker peaks for motif enrichment -----------------------------------
marker_peak_path <- file.path(epi_project_dir, "MarkerPeaks", "markersPeaks_sub.annotation.rds")
markers_peaks <- if (file.exists(marker_peak_path)) {
  readRDS(marker_peak_path)
} else {
  getMarkerFeatures(
    ArchRProj = proj_epi,
    useMatrix = "PeakMatrix",
    groupBy = "sub.annotation",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
}

# Step 4: Motif enrichment across marker peaks -------------------------------------
enrich_motifs <- peakAnnoEnrichment(
  seMarker = markers_peaks,
  ArchRProj = proj_epi,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)
heatmap_matrix <- plotEnrichHeatmap(
  enrich_motifs,
  n = 50,
  transpose = FALSE,
  pal = paletteContinuous(set = "solarExtra", n = 100),
  returnMatrix = TRUE
)
pdf(file.path(output_dir, "Motif_enrichment_heatmap.pdf"), width = 8, height = 10)
col_fun <- colorRamp2(c(min(heatmap_matrix), median(heatmap_matrix), max(heatmap_matrix)),
                      c("white", "#FDBE85", "#D94701"))
Heatmap(heatmap_matrix, name = "Enrichment", col = col_fun, cluster_columns = FALSE)
dev.off()

# Step 5: Correlate gene scores with motif deviations ------------------------------
cor_gs_motif <- correlateMatrices(
  ArchRProj = proj_epi,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)
cor_df <- as.data.frame(cor_gs_motif)
write_csv(cor_df, file.path(output_dir, "gene_motif_correlations.csv"))

# Step 6: Cluster motifs based on enrichment profiles ------------------------------
dist_mat <- dist(heatmap_matrix)
hc <- hclust(dist_mat, method = "ward.D2")
clusters <- cutree(hc, k = 9)
cluster_summary <- tibble(motif = rownames(heatmap_matrix), cluster = clusters)
write_csv(cluster_summary, file.path(output_dir, "motif_clusters.csv"))
