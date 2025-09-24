#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggVennDiagram)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(parallel)
  library(TeachingDemos)
})

# ---- Configuration ------------------------------------------------------
config <- list(
  project_root = getwd(),
  renv_activate = file.path("renv", "activate.R"),
  epi_project_dir = "/path/to/Archr_TotalOutput_Global_Epi",
  output_name = file.path("Archr_Result", "Archr_Epi_Markers_Peaks"),
  genome = "hg38",
  threads = min(16L, max(1L, parallel::detectCores(logical = TRUE))),
  log_subdir = "logs",
  plots_subdir = "plots",
  data_subdir = "data",
  function_script = "/path/to/Rfile/Function_list.R",
  rna_deg_csv = "/path/to/epi.deg.0.01.csv",
  rna_deg_rds = "/path/to/Epi_sub_markers.rds",
  rna_cluster_order = c("TC.BRS.high1", "TC.BRS.high2", "TC.BRS.mid1", "TC.BRS.mid2", "TC.BRS.low", "ATC.ERK.low", "ATC.ERK.high"),
  atac_annotation_levels = c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high"),
  marker_cutoff = list(fdr = 0.01, log2fc = 1.25),
  venn_groups = c("ERK.high", "ERK.low", "TDS.high", "TDS.mid2", "TDS.low"),
  macs2_path = "/usr/local/bin/macs2",
  group_by = "sub.annotation"
)

config$output_dir <- normalizePath(
  file.path(config$project_root, config$output_name),
  mustWork = FALSE
)
config$plots_dir <- normalizePath(
  file.path(config$output_dir, config$plots_subdir),
  mustWork = FALSE
)
config$data_dir <- normalizePath(
  file.path(config$output_dir, config$data_subdir),
  mustWork = FALSE
)
config$log_path <- normalizePath(
  file.path(config$output_dir, config$log_subdir, "archr_epi_markers_peaks.log"),
  mustWork = FALSE
)
config$renv_activate <- file.path(config$project_root, config$renv_activate)

if (!dir.exists(config$epi_project_dir)) {
  stop(
    "Epithelial ArchR project not found. Update config$epi_project_dir. Current value: ",
    config$epi_project_dir
  )
}
config$epi_project_dir <- normalizePath(config$epi_project_dir, mustWork = TRUE)

message("Markers and peak calling outputs will be written to: ", config$output_dir)

if (file.exists(config$renv_activate)) {
  message("Activating renv from ", config$renv_activate)
  source(config$renv_activate)
}

addArchRGenome(config$genome)
addArchRThreads(threads = config$threads)


dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(config$log_path), recursive = TRUE, showWarnings = FALSE)

original_wd <- getwd()
setwd(config$output_dir)
on.exit(setwd(original_wd), add = TRUE)

txtStart(config$log_path)
on.exit(try(txtStop(), silent = TRUE), add = TRUE)

if (file.exists(config$function_script)) {
  source(config$function_script)
} else {
  warning("Custom helper script not found: ", config$function_script, "\nContinuing without external helper functions.")
}

# ---- Helper utilities ---------------------------------------------------
save_plot <- function(plot, filename, width = 6, height = 5) {
  filepath <- file.path(config$plots_dir, filename)
  ggsave(filepath, plot = plot, width = width, height = height, limitsize = FALSE)
  message("Saved plot: ", filepath)
}

scale_matrix <- function(mat, limit = 4) {
  scaled <- t(scale(t(mat)))
  scaled[scaled > limit] <- limit
  scaled[scaled < -limit] <- -limit
  scaled[is.na(scaled)] <- 0
  scaled
}

get_marker_cutoff <- function() {
  paste0("FDR <= ", config$marker_cutoff$fdr, " & Log2FC >= ", config$marker_cutoff$log2fc)
}

# ---- Load ArchR project -------------------------------------------------
message("Loading epithelial ArchR project from ", config$epi_project_dir)
proj <- loadArchRProject(
  path = config$epi_project_dir,
  force = FALSE,
  showLogo = FALSE
)

if (!config$group_by %in% colnames(proj@cellColData)) {
  stop("Project is missing column '", config$group_by, "'. Ensure annotation pipeline has been run.")
}

proj <- addImputeWeights(proj)
weights <- getImputeWeights(proj)

# ---- UMAP overviews -----------------------------------------------------
p_umap_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 9, legendPosition = "right")

p_umap_group <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "group",
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 9, legendPosition = "right")

p_umap_annotation <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = config$group_by,
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 8, legendPosition = "right")

save_plot((p_umap_sample | p_umap_group | p_umap_annotation), "UMAP-Sample-Group-Annotation.pdf", width = 15, height = 5)

# ---- Gene score matrix --------------------------------------------------
message("Retrieving GeneScore matrix")
gene_score_se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gene_score <- assays(gene_score_se)$GeneScoreMatrix
rownames(gene_score) <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
gene_score <- as.matrix(gene_score)

# ---- RNA DEG import -----------------------------------------------------
load_deg_table <- function(path, source_label, cluster_order) {
  if (!file.exists(path)) {
    warning(source_label, " DEG file not found: ", path)
    return(NULL)
  }
  tbl <- if (stringr::str_ends(path, ".rds")) {
    readRDS(path)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
  tbl <- tbl %>% dplyr::mutate(cluster = gsub("_", ".", cluster))
  if (!is.null(cluster_order)) {
    tbl$cluster <- factor(tbl$cluster, levels = cluster_order)
  }
  tbl$source <- source_label
  tbl
}

deg_external <- load_deg_table(config$rna_deg_csv, "Seoul", config$rna_cluster_order)
deg_internal <- load_deg_table(config$rna_deg_rds, "Internal", config$rna_cluster_order)

combined_deg <- purrr::compact(list(deg_external, deg_internal))
if (length(combined_deg) == 0) {
  warning("No RNA DEG tables were loaded. Skipping DEG comparisons.")
}

# ---- ArchR marker features ----------------------------------------------
message("Computing GeneScore markers for ", config$group_by)
markers_sub <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = config$group_by,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markers_sub, file.path(config$data_dir, "markersGS_sub.RDS"))
cutoff_str <- get_marker_cutoff()
marker_list_sub <- getMarkers(markers_sub, cutOff = cutoff_str)
marker_list_sub <- purrr::compact(marker_list_sub)

heatmap_sub <- markerHeatmap(
  seMarker = markers_sub,
  cutOff = cutoff_str,
  transpose = FALSE,
  nLabel = 5
)
ComplexHeatmap::draw(heatmap_sub, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(
  heatmap_sub,
  name = "GeneScores-Marker-Heatmap_sub.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 10,
  height = 13
)

if (file.exists(file.path(config$data_dir, "markersGS_Clusters.RDS"))) {
  markers_clusters <- readRDS(file.path(config$data_dir, "markersGS_Clusters.RDS"))
} else {
  message("Cluster markers not found on disk. Computing using current project...")
  markers_clusters <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  saveRDS(markers_clusters, file.path(config$data_dir, "markersGS_Clusters.RDS"))
}
markerHeat_clusters <- markerHeatmap(
  seMarker = markers_clusters,
  cutOff = cutoff_str,
  transpose = FALSE,
  nLabel = 2,
  clusterCols = FALSE,
  binaryClusterRows = TRUE
)
ComplexHeatmap::draw(markerHeat_clusters, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(
  markerHeat_clusters,
  name = "GeneScores-Marker-Heatmap_clusters.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 10,
  height = 13
)

# ---- RNA/ATAC overlap summaries -----------------------------------------
if (length(combined_deg) > 0) {
  marker_gene_lists <- lapply(marker_list_sub, function(df) unique(df$name))
  deg_lists <- lapply(combined_deg, function(tbl) split(tbl$gene, tbl$cluster))
  deg_lists <- purrr::map(deg_lists, purrr::compact)

  overlap_dir <- file.path(config$data_dir, "overlap")
  dir.create(overlap_dir, showWarnings = FALSE)

  purrr::walk2(deg_lists, names(deg_lists), function(cluster_map, label) {
    overlap_mat <- matrix(0, nrow = length(marker_gene_lists), ncol = length(cluster_map), dimnames = list(names(marker_gene_lists), names(cluster_map)))
    for (i in seq_along(marker_gene_lists)) {
      for (j in seq_along(cluster_map)) {
        overlap_mat[i, j] <- length(intersect(marker_gene_lists[[i]], cluster_map[[j]]))
      }
    }
    prop_by_rna <- prop.table(overlap_mat, margin = 2)
    prop_df <- as.data.frame(as.table(prop_by_rna))
    p_overlap <- ggplot(prop_df, aes(Var1, Var2, fill = Freq)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "white", high = "red") +
      labs(title = paste("DEG overlap (", label, ")"), x = "ATAC markers", y = "RNA clusters", fill = "Proportion") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
    save_plot(p_overlap, sprintf("Overlap_%s.pdf", label), width = 7, height = 5)
    write.csv(overlap_mat, file.path(overlap_dir, sprintf("overlap_counts_%s.csv", label)))
  })

  venn_groups <- intersect(config$venn_groups, names(marker_gene_lists))
  if (length(venn_groups) >= 2) {
    venn_data <- marker_gene_lists[venn_groups]
    venn_plot <- ggVennDiagram(venn_data, label_alpha = 0) +
      coord_flip() +
      scale_fill_gradient(low = "grey95", high = "lightsalmon1") +
      scale_color_manual(values = rep("black", length(venn_groups))) +
      ggtitle("Marker gene overlap") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
    save_plot(venn_plot, "Venn-Marker-Gene-Overlap.pdf", width = 5, height = 5)
  }
}

# ---- Gene-based heatmaps -------------------------------------------------
if (exists("CellFeatureHeatmap")) {
  message("Creating gene-by-cell heatmap using helper function")
  annotation_df <- proj@cellColData %>%
    as.data.frame() %>%
    dplyr::select(group, Sample)
  annotation_df$sub.annotation <- proj@cellColData[[config$group_by]]
  annotation_df <- annotation_df %>%
    dplyr::mutate(
      group = gsub("-", "_", group),
      sub.annotation = factor(sub.annotation, levels = config$atac_annotation_levels),
      Sample = factor(Sample, levels = sort(unique(Sample)))
    )

  sample_levels <- levels(annotation_df$Sample)
  if (length(sample_levels) <= 12) {
    sample_colors <- brewer.pal(max(8, length(sample_levels)), "Paired")[seq_along(sample_levels)]
  } else {
    sample_colors <- grDevices::colorRampPalette(brewer.pal(12, "Paired"))(length(sample_levels))
  }
  annotation_colors <- list(
    sub.annotation = setNames(brewer.pal(6, "Set1"), config$atac_annotation_levels),
    group = c(BRAF_PTC = "darkgreen", BRAF_ATC = "red", RAS_FTC = "orange", RAS_ATC = "purple", NBNR_PDTC = "cornflowerblue", Normal = "black"),
    Sample = setNames(sample_colors, sample_levels)
  )

  feature_tbl <- purrr::imap(marker_list_sub, function(df, label) {
    head(df$name, 50) %>%
      tibble::tibble(gene = ., cluster = label)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(cluster = factor(cluster, levels = intersect(config$atac_annotation_levels, unique(cluster))))

  if (nrow(feature_tbl) > 0) {
    feature_df <- data.frame(cluster = feature_tbl$cluster, row.names = feature_tbl$gene)
    CellFeatureHeatmap(
      matrix = gene_score,
      feature = feature_df,
      cell_annotation = annotation_df,
      annotation_color = annotation_colors,
      normalize = TRUE,
      scaling = TRUE,
      legend_title = "Gene score",
      name = "Heatmap_CellGene"
    )
  }
} else {
  message("CellFeatureHeatmap helper not available; skipping cell-by-gene heatmap.")
}

if (exists("ClusterFeatureHeatmap")) {
  message("Creating pseudobulk cluster heatmap using helper function")
  feature_tbl <- purrr::imap(marker_list_sub, function(df, label) {
    head(df$name, 50) %>%
      tibble::tibble(gene = ., cluster = label)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(cluster = factor(cluster, levels = config$atac_annotation_levels))

  feature_tbl <- dplyr::filter(feature_tbl, !is.na(cluster))

  if (nrow(feature_tbl) > 0) {
    annotation_df <- data.frame(sub.annotation = proj@cellColData[[config$group_by]])
    annotation_df$sub.annotation <- factor(annotation_df$sub.annotation, levels = config$atac_annotation_levels)

    ClusterFeatureHeatmap(
      matrix = gene_score,
      feature = data.frame(cluster = feature_tbl$cluster, row.names = feature_tbl$gene),
      feature_anno = head(feature_tbl$gene, 200),
      cell_annotation = annotation_df,
      annotation_color = list(sub.annotation = setNames(brewer.pal(6, "Set1"), config$atac_annotation_levels)),
      scaling = TRUE,
      legend_title = "Gene score",
      name = "Heatmap_CellTypeGene"
    )
  } else {
    message("No markers available for cluster-level heatmap.")
  }
} else {
  message("ClusterFeatureHeatmap helper not available; skipping cluster-level heatmap.")
}

# ---- Pseudobulk matrix export -------------------------------------------
annotation_vec <- proj@cellColData[[config$group_by]]
annotation_levels <- intersect(config$atac_annotation_levels, unique(annotation_vec))

pseudo_bulk <- vapply(annotation_levels, function(label) {
  cells <- proj$cellNames[annotation_vec == label]
  if (length(cells) == 0) {
    return(rep(NA_real_, nrow(gene_score)))
  }
  rowMeans(gene_score[, cells, drop = FALSE])
}, numeric(nrow(gene_score)))
rownames(pseudo_bulk) <- rownames(gene_score)

pseudo_bulk_scaled <- scale_matrix(pseudo_bulk, limit = 4)

write.csv(pseudo_bulk, file.path(config$data_dir, "CellTypeGene_raw.csv"), quote = FALSE)
write.csv(pseudo_bulk_scaled, file.path(config$data_dir, "CellTypeGene_scaled.csv"), quote = FALSE)

# ---- Peak calling -------------------------------------------------------
message("Running pseudo-bulk coverage and peak calling")
proj_with_coverage <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = config$group_by,
  threads = config$threads,
  force = TRUE
)

if (!file.exists(config$macs2_path)) {
  stop("macs2 executable not found at ", config$macs2_path)
}

proj_with_peaks <- addReproduciblePeakSet(
  ArchRProj = proj_with_coverage,
  groupBy = config$group_by,
  pathToMacs2 = config$macs2_path
)

proj_with_peaks <- addPeakMatrix(proj_with_peaks)

# ---- Save project -------------------------------------------------------
message("Saving updated ArchR project (with peak matrix)")
saveArchRProject(ArchRProj = proj_with_peaks, load = FALSE)

message("Markers and peak calling pipeline completed at ", Sys.time())
