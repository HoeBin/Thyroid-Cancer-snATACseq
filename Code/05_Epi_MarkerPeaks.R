#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(grid)
  library(parallel)
  library(TeachingDemos)
})

# ---- Configuration ------------------------------------------------------
config <- list(
  project_root = getwd(),
  renv_activate = file.path("renv", "activate.R"),
  epi_project_dir = "/path/to/Archr_TotalOutput_Global_Epi",
  output_name = file.path("Archr_Result", "Archr_Epi_MarkerPeaks"),
  genome = "hg38",
  threads = min(16L, max(1L, parallel::detectCores(logical = TRUE))),
  log_subdir = "logs",
  plots_subdir = "plots",
  data_subdir = "data",
  marker_cutoff = list(fdr = 0.1, log2fc = 0.5),
  track_genes = tibble::tribble(
    ~gene, ~upstream, ~downstream,
    "BRAF", 300000, 50000,
    "HRAS", 10000, 5000,
    "NRAS", 20000, 20000,
    "TERT", 80000, 20000
  ),
  track_groupings = c("sub.annotation", "Sample", "group"),
  sub_annotation_levels = c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high"),
  group_levels = c("BRAF-PTC", "BRAF-ATC", "RAS-FTC", "RAS-ATC", "NBNR-PDTC"),
  sample_levels = c("ATAC3", "ATAC4", "ATAC1", "ATAC2", "ATAC8", "ATAC9", "ATAC10", "ATAC5", "ATAC6", "ATAC7", "PDTC"),
  pairwise_tests = list(
    list(name = "ERKhigh_vs_ERKlow", use = "ERK.high", background = "ERK.low"),
    list(name = "TDShigh_vs_TDSlow", use = "TDS.high", background = "TDS.low")
  )
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
  file.path(config$output_dir, config$log_subdir, "archr_epi_markerpeaks.log"),
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

message("Marker-peak outputs will be written to: ", config$output_dir)

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

# ---- Helpers ------------------------------------------------------------
save_plot <- function(plot, filename, width = 6, height = 5) {
  filepath <- file.path(config$plots_dir, filename)
  ggsave(filepath, plot = plot, width = width, height = height, limitsize = FALSE)
  message("Saved plot: ", filepath)
}

marker_cutoff_string <- function(fdr = config$marker_cutoff$fdr, log2fc = config$marker_cutoff$log2fc) {
  paste0("FDR <= ", fdr, " & Log2FC >= ", log2fc)
}

scale_matrix <- function(mat, limit = 4) {
  scaled <- t(scale(t(mat)))
  scaled[scaled > limit] <- limit
  scaled[scaled < -limit] <- -limit
  scaled[is.na(scaled)] <- 0
  scaled
}

# ---- Load project -------------------------------------------------------
message("Loading epithelial ArchR project from ", config$epi_project_dir)
proj <- loadArchRProject(
  path = config$epi_project_dir,
  force = FALSE,
  showLogo = FALSE
)

proj <- addImputeWeights(proj)
weights <- getImputeWeights(proj)

cell_meta <- proj@cellColData %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    sub.annotation = factor(sub.annotation, levels = config$sub_annotation_levels),
    group = factor(group, levels = config$group_levels),
    Sample = factor(Sample, levels = config$sample_levels)
  )

# ---- UMAP overview ------------------------------------------------------
p_umap_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 8, legendPosition = "right")

p_umap_group <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "group",
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 8, legendPosition = "right")

p_umap_sub <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "sub.annotation",
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 8, legendPosition = "right")

save_plot((p_umap_sample | p_umap_group | p_umap_sub), "UMAP-Sample-Group-Sub.pdf", width = 15, height = 5)

# ---- Track plotting -----------------------------------------------------
plot_browser_tracks <- function(group_by, use_levels) {
  message("Generating browser tracks for grouping: ", group_by)
  plot_dir <- file.path(config$plots_dir, paste0("tracks_", group_by))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  valid_levels <- stats::na.omit(use_levels)
  present_levels <- unique(proj@cellColData[[group_by]])
  present_levels <- present_levels[!is.na(present_levels)]
  valid_levels <- valid_levels[valid_levels %in% present_levels]
  if (length(valid_levels) == 0) {
    warning("No matching levels found for grouping ", group_by, ". Skipping track plots.")
    return(invisible(NULL))
  }

  plot_list <- list()
  for (gene_info in split(config$track_genes, seq_len(nrow(config$track_genes)))) {
    gene <- gene_info$gene
    upstream <- gene_info$upstream
    downstream <- gene_info$downstream

    track <- plotBrowserTrack(
      ArchRProj = proj,
      groupBy = group_by,
      useGroups = valid_levels,
      geneSymbol = gene,
      plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
      upstream = upstream,
      downstream = downstream
    )
    plot_list[[paste(group_by, gene, sep = "_")]] <- track[[1]]
  }

  pdf_path <- file.path(plot_dir, paste0("Plot-Tracks-", group_by, ".pdf"))
  grDevices::pdf(pdf_path, width = 8, height = 8)
  purrr::walk(plot_list, grid::grid.draw)
  grDevices::dev.off()
  message("Saved browser tracks: ", pdf_path)
}

plot_browser_tracks("sub.annotation", config$sub_annotation_levels)
plot_browser_tracks("Sample", config$sample_levels)
plot_browser_tracks("group", config$group_levels)

# ---- Marker peaks -------------------------------------------------------
message("Computing marker peaks across groupings")
marker_groupings <- list(
  sub.annotation = config$sub_annotation_levels,
  Sample = config$sample_levels,
  Clusters = NULL
)

marker_results <- purrr::imap(marker_groupings, function(levels, group_by) {
  message("  - groupBy = ", group_by)
  markers <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  output_path <- file.path(config$data_dir, sprintf("markersPeaks_%s.rds", group_by))
  saveRDS(markers, output_path)
  list(markers = markers, path = output_path)
})

markers_sub <- marker_results$sub.annotation$markers
cutoff_str <- marker_cutoff_string()
marker_list_sub <- getMarkers(markers_sub, cutOff = cutoff_str)
marker_list_sub <- purrr::compact(marker_list_sub)

peak_counts <- purrr::map_int(marker_list_sub, nrow)
if (length(peak_counts) > 0) {
  count_df <- tibble::tibble(type = names(peak_counts), NumPeaks = peak_counts) %>%
    dplyr::mutate(type = factor(type, levels = config$sub_annotation_levels)) %>%
    tidyr::drop_na(type)
  p_counts <- ggplot(count_df, aes(x = type, y = NumPeaks)) +
    geom_col(fill = "steelblue") +
    theme_cowplot() +
    labs(x = NULL, y = "Marker peaks", title = "Marker peaks by sub annotation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_counts, "MarkerPeaks_Counts.pdf", width = 6, height = 4)
}

heatmap_peaks <- plotMarkerHeatmap(
  seMarker = markers_sub,
  cutOff = cutoff_str,
  transpose = FALSE,
  clusterCols = TRUE,
  nLabel = 2
)
plotPDF(
  heatmap_peaks,
  name = "MarkerPeaks-Heatmap-sub.annotation.pdf",
  ArchRProj = proj,
  width = 10,
  height = 12,
  addDOC = FALSE
)

# ---- MA / Volcano examples ---------------------------------------------
if (length(marker_list_sub) > 0) {
  example_group <- names(marker_list_sub)[1]
  ma_plot <- plotMarkers(
    seMarker = markers_sub,
    name = example_group,
    cutOff = marker_cutoff_string(fdr = 0.1, log2fc = 1),
    plotAs = "MA"
  ) + theme_ArchR(legendTextSize = 8, legendPosition = "right")
  vol_plot <- plotMarkers(
    seMarker = markers_sub,
    name = example_group,
    cutOff = marker_cutoff_string(fdr = 0.1, log2fc = 1),
    plotAs = "Volcano"
  ) + theme_ArchR(legendTextSize = 8, legendPosition = "right")
  save_plot(ma_plot | vol_plot, sprintf("Markers-MA-Volcano_%s.pdf", example_group), width = 10, height = 4)
}

# ---- Pairwise marker tests ----------------------------------------------
pairwise_dir <- file.path(config$data_dir, "pairwise")
dir.create(pairwise_dir, showWarnings = FALSE)

for (comp in config$pairwise_tests) {
  message("Running pairwise comparison: ", comp$name)
  marker_test <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "sub.annotation",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = comp$use,
    bgdGroups = comp$background
  )

  test_path <- file.path(pairwise_dir, paste0("markerTest_", comp$name, ".rds"))
  saveRDS(marker_test, test_path)

  marker_list <- getMarkers(marker_test, cutOff = marker_cutoff_string(fdr = 0.1, log2fc = 1))
  if (length(marker_list) == 0) {
    next
  }

  export_dir <- file.path(pairwise_dir, comp$name)
  dir.create(export_dir, showWarnings = FALSE)

  purrr::iwalk(marker_list, function(df, group_name) {
    if (nrow(df) == 0) {
      return()
    }
    bed_df <- df %>%
      dplyr::filter(FDR <= 0.01) %>%
      dplyr::mutate(
        peak = paste0(seqnames, ":", start, "-", end)
      )
    if (nrow(bed_df) == 0) {
      return()
    }
    up_file <- file.path(export_dir, paste0(group_name, "_up.bed"))
    down_file <- file.path(export_dir, paste0(group_name, "_down.bed"))

    bed_df %>%
      dplyr::filter(Log2FC >= 1) %>%
      dplyr::select(seqnames, start, end) %>%
      readr::write_tsv(up_file, col_names = FALSE)

    bed_df %>%
      dplyr::filter(Log2FC <= -1) %>%
      dplyr::select(seqnames, start, end) %>%
      readr::write_tsv(down_file, col_names = FALSE)
  })
}

# ---- Heatmap preparation utilities --------------------------------------
prepare_peak_matrix <- function(proj, peak_list, top_n = 100) {
  peak_df <- purrr::imap_dfr(peak_list, function(df, group_name) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    df <- as.data.frame(df)
    df$group_name <- group_name
    df
  })
  if (nrow(peak_df) == 0) {
    return(list(matrix = NULL, feature = NULL))
  }
  peak_df <- peak_df %>%
    dplyr::mutate(peak = paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::group_by(group_name) %>%
    dplyr::arrange(dplyr::desc(Log2FC)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  peak_mat <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
  peak_access <- assays(peak_mat)$PeakMatrix

  all_peaks <- paste0(seqnames(proj@peakSet), ":", start(proj@peakSet), "-", end(proj@peakSet))
  rownames(peak_access) <- all_peaks

  selected <- peak_access[peak_df$peak, , drop = FALSE]
  feature_df <- data.frame(cluster = peak_df$group_name, row.names = peak_df$peak)

  list(
    matrix = selected,
    feature = feature_df
  )
}

markers_peaks_top <- getMarkers(markers_sub, cutOff = marker_cutoff_string(fdr = 0.001, log2fc = 0.5))
peak_heatmap_data <- prepare_peak_matrix(proj, markers_peaks_top)

if (!is.null(peak_heatmap_data$matrix)) {
  annotation_df <- cell_meta %>%
    dplyr::select(cell, sub.annotation, group, Sample) %>%
    dplyr::mutate(group = gsub("-", "_", group))
  rownames(annotation_df) <- annotation_df$cell
  annotation_df <- annotation_df[colnames(peak_heatmap_data$matrix), , drop = FALSE]

  sample_levels <- config$sample_levels
  if (length(sample_levels) <= 12) {
    sample_colors <- brewer.pal(max(8, length(sample_levels)), "Paired")[seq_along(sample_levels)]
  } else {
    sample_colors <- grDevices::colorRampPalette(brewer.pal(12, "Paired"))(length(sample_levels))
  }
  annotation_colors <- list(
    sub.annotation = setNames(brewer.pal(6, "Set1"), config$sub_annotation_levels),
    group = c(BRAF_PTC = "darkgreen", BRAF_ATC = "red", RAS_FTC = "orange", RAS_ATC = "purple", NBNR_PDTC = "cornflowerblue", Normal = "black"),
    Sample = setNames(sample_colors, sample_levels)
  )

  if (exists("CellFeatureHeatmap", mode = "function")) {
    CellFeatureHeatmap_fn <- get("CellFeatureHeatmap", mode = "function")
    CellFeatureHeatmap_fn(
      matrix = peak_heatmap_data$matrix,
      feature = peak_heatmap_data$feature,
      cell_annotation = annotation_df,
      annotation_color = annotation_colors,
      normalize = TRUE,
      scaling = TRUE,
      legend_title = "Accessibility",
      name = "Heatmap_PeaksCells"
    )
  } else {
    message("CellFeatureHeatmap function not available; skipping cell-level peak heatmap.")
  }
}

# ---- Save project -------------------------------------------------------
message("Saving ArchR project (no modification made to peak matrices)")
saveArchRProject(ArchRProj = proj, load = FALSE)

message("Marker peak visualisation pipeline completed at ", Sys.time())
