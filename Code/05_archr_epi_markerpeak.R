suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(cowplot)
})

# Step 1: Configure directories and load the project -------------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))
output_dir <- file.path(epi_project_dir, "MarkerPeaks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
proj_epi <- addImputeWeights(proj_epi)

orders_celltype <- c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high")
orders_group <- c("BRAF-PTC", "BRAF-ATC", "RAS-FTC", "RAS-ATC")
orders_sample <- c("ATAC3", "ATAC4", "ATAC1", "ATAC2", "ATAC8", "ATAC9", "ATAC10", "ATAC6", "ATAC7", "PDTC")

# Step 2: Visualize UMAP embedding summaries ---------------------------------------
embedding_fields <- c("Sample", "group", "mut", "sub.annotation")
embedding_plots <- lapply(embedding_fields, function(field) {
  plotEmbedding(proj_epi, colorBy = "cellColData", name = field, embedding = "UMAP") +
    theme_ArchR(legendTextSize = 8, legendPosition = "right") +
    ggtitle(field)
})
plotPDF(plotlist = embedding_plots, name = "UMAP_overview.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 5, height = 5)

# Step 3: Generate browser tracks for mutation-centric genes -----------------------
plot_tracks <- function(group_field, groups, genes, filename) {
  plots <- map(genes, function(gene) {
    plotBrowserTrack(
      ArchRProj = proj_epi,
      groupBy = group_field,
      useGroups = groups,
      geneSymbol = gene,
      plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
      upstream = 50000,
      downstream = 50000
    )
  })
  plotPDF(plotList = plots, name = filename, ArchRProj = proj_epi, addDOC = FALSE, width = 8, height = 8)
}

driver_genes <- c("BRAF", "HRAS", "NRAS", "TERT")
plot_tracks("sub.annotation", orders_celltype, driver_genes, "Tracks_subannotation.pdf")
plot_tracks("Sample", orders_sample, driver_genes, "Tracks_sample.pdf")
plot_tracks("group", orders_group, driver_genes, "Tracks_group.pdf")

# Step 4: Identify marker peaks for multiple groupings -----------------------------
marker_configs <- list(
  sub_annotation = list(groupBy = "sub.annotation", cut = "FDR <= 0.1 & Log2FC >= 0.5"),
  sample = list(groupBy = "Sample", cut = "FDR <= 0.1 & Log2FC >= 0.5"),
  clusters = list(groupBy = "Clusters", cut = "FDR <= 0.1 & Log2FC >= 0.5")
)

marker_results <- map(marker_configs, function(cfg) {
  markers <- getMarkerFeatures(
    ArchRProj = proj_epi,
    useMatrix = "PeakMatrix",
    groupBy = cfg$groupBy,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  saveRDS(markers, file.path(output_dir, sprintf("markersPeaks_%s.rds", cfg$groupBy)))
  cut <- cfg$cut
  heatmap <- plotMarkerHeatmap(
    seMarker = markers,
    cutOff = cut,
    transpose = FALSE,
    nLabel = 2,
    clusterCols = TRUE
  )
  plotPDF(heatmap, name = sprintf("MarkerPeaks_%s.pdf", cfg$groupBy), ArchRProj = proj_epi, addDOC = FALSE, width = 8, height = 10)
  list(se = markers, heatmap = heatmap)
})

# Step 5: Explore MA/Volcano plots for selected comparisons ------------------------
markers_sub <- marker_results$sub_annotation$se
groups_to_plot <- c("ERK.high", "ERK.low")
for (grp in groups_to_plot) {
  p_ma <- plotMarkers(seMarker = markers_sub, name = grp, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
  p_volcano <- plotMarkers(seMarker = markers_sub, name = grp, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
  plotPDF(p_ma, p_volcano, name = sprintf("Markers_%s.pdf", grp), ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 6)
}

# Step 6: Optional marker peak browser for high-confidence loci --------------------
features <- getMarkers(markers_sub, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)["ERK.high"]
feature_tracks <- map(driver_genes, function(gene) {
  plotBrowserTrack(
    ArchRProj = proj_epi,
    groupBy = "sub.annotation",
    useGroups = orders_celltype,
    geneSymbol = gene,
    features = features,
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    upstream = 50000,
    downstream = 50000
  )
})
plotPDF(feature_tracks, name = "Tracks_with_marker_peaks.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 8, height = 8)

# Step 7: Pairwise peak tests for contrasts of interest ---------------------------
pairwise_tests <- list(
  ERK_high_vs_low = list(use = "ERK.high", bgd = "ERK.low"),
  TDS_high_vs_low = list(use = "TDS.high", bgd = "TDS.low")
)

walk(names(pairwise_tests), function(name) {
  cfg <- pairwise_tests[[name]]
  marker_test <- getMarkerFeatures(
    ArchRProj = proj_epi,
    useMatrix = "PeakMatrix",
    groupBy = "sub.annotation",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = cfg$use,
    bgdGroups = cfg$bgd
  )
  saveRDS(marker_test, file.path(output_dir, sprintf("pairwise_%s.rds", name)))
})

