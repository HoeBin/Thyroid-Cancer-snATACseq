#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(data.table)
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
  output_name = file.path("Archr_Result", "Archr_Epi_Motif"),
  genome = "hg38",
  threads = min(16L, max(1L, parallel::detectCores(logical = TRUE))),
  log_subdir = "logs",
  plots_subdir = "plots",
  data_subdir = "data",
  sub_annotation_levels = c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high"),
  group_levels = c("BRAF-PTC", "BRAF-ATC", "RAS-FTC", "RAS-ATC", "NBNR-PDTC"),
  sample_levels = c("ATAC3", "ATAC4", "ATAC1", "ATAC2", "ATAC8", "ATAC9", "ATAC10", "ATAC5", "ATAC6", "ATAC7", "PDTC"),
  pairwise_tests = list(
    list(name = "ERKhigh_vs_ERKlow", file = "markerTest_ErkHighLow.rds", up_cutoff = list(fdr = 0.01, log2fc = 1), down_cutoff = list(fdr = 0.01, log2fc = -1)),
    list(name = "TDShigh_vs_TDSlow", file = "markerTest_TdsHighLow.rds", up_cutoff = list(fdr = 0.01, log2fc = 2), down_cutoff = list(fdr = 0.01, log2fc = -2))
  ),
  motif_cutoff = list(fdr = 0.05, log2fc = 1),
  top_motifs = 50,
  homer_dir = NULL
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
  file.path(config$output_dir, config$log_subdir, "archr_epi_motif.log"),
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

message("Motif analysis outputs will be written to: ", config$output_dir)

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

cutoff_string <- function(cutoff) {
  if (cutoff$log2fc >= 0) {
    paste0("FDR <= ", cutoff$fdr, " & Log2FC >= ", cutoff$log2fc)
  } else {
    paste0("FDR <= ", cutoff$fdr, " & Log2FC <= ", cutoff$log2fc)
  }
}

scale_matrix <- function(mat, limit = 4) {
  scaled <- t(scale(t(mat)))
  scaled[scaled > limit] <- limit
  scaled[scaled < -limit] <- -limit
  scaled[is.na(scaled)] <- 0
  scaled
}

prepare_enrichment_df <- function(se_obj) {
  df <- data.frame(TF = rownames(se_obj), mlog10Padj = assay(se_obj)[, 1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE), ]
  df$rank <- seq_len(nrow(df))
  df
}

plot_enrichment_scatter <- function(df, title = "Motif enrichment") {
  ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[seq_len(min(30, nrow(df))), ],
      aes(label = TF),
      size = 3,
      nudge_x = 2,
      color = "black"
    ) +
    theme_ArchR() +
    ylab("-log10(FDR) motif enrichment") +
    xlab("Rank sorted TFs") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    ggtitle(title)
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

# ---- QC UMAP ------------------------------------------------------------
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

# ---- Add motif annotations ----------------------------------------------
message("Adding motif annotations (cisbp, EncodeTFBS)")
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = FALSE)
proj <- addArchRAnnotations(ArchRProj = proj, collection = "EncodeTFBS", force = FALSE)

# ---- Pairwise motif enrichment -----------------------------------------
pairwise_dir <- file.path(config$data_dir, "pairwise")
dir.create(pairwise_dir, showWarnings = FALSE)

for (comp in config$pairwise_tests) {
  test_path <- file.path(config$epi_project_dir, "Output", comp$file)
  if (!file.exists(test_path)) {
    warning("Pairwise marker test not found: ", test_path)
    next
  }
  marker_test <- readRDS(test_path)

  up_res <- peakAnnoEnrichment(
    seMarker = marker_test,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = cutoff_string(comp$up_cutoff)
  )
  down_res <- peakAnnoEnrichment(
    seMarker = marker_test,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = cutoff_string(comp$down_cutoff)
  )

  up_df <- prepare_enrichment_df(up_res)
  down_df <- prepare_enrichment_df(down_res)

  saveRDS(list(up = up_res, down = down_res), file.path(pairwise_dir, paste0("motifEnrichment_", comp$name, ".rds")))
  readr::write_csv(up_df, file.path(pairwise_dir, paste0("motifs_up_", comp$name, ".csv")))
  readr::write_csv(down_df, file.path(pairwise_dir, paste0("motifs_down_", comp$name, ".csv")))

  p_up <- plot_enrichment_scatter(up_df, paste0(comp$name, " up motifs"))
  p_down <- plot_enrichment_scatter(down_df, paste0(comp$name, " down motifs"))
  save_plot(p_up | p_down, paste0("Motif-Enrichment-Scatter_", comp$name, ".pdf"), width = 10, height = 4)
}

# ---- Marker peak motif enrichment ---------------------------------------
markers_path <- file.path(config$epi_project_dir, "Output", "markersPeaks_sub.annotation.rds")
if (!file.exists(markers_path)) {
  stop("Marker peaks file not found: ", markers_path)
}
markers_peaks <- readRDS(markers_path)
marker_list <- getMarkers(markers_peaks, cutOff = cutoff_string(config$motif_cutoff))
marker_list <- purrr::compact(marker_list)

peak_counts <- purrr::map_int(marker_list, nrow)
if (length(peak_counts) > 0) {
  count_df <- tibble::tibble(type = names(peak_counts), NumPeaks = peak_counts) %>%
    dplyr::mutate(type = factor(type, levels = config$sub_annotation_levels)) %>%
    tidyr::drop_na(type)
  p_counts <- ggplot(count_df, aes(x = type, y = NumPeaks)) +
    geom_col(fill = "steelblue") +
    theme_cowplot() +
    labs(x = NULL, y = "Marker peaks", title = "Marker peaks per sub annotation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_counts, "MarkerPeaks_Counts.pdf", width = 6, height = 4)
}

motif_enrich <- peakAnnoEnrichment(
  seMarker = markers_peaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = cutoff_string(config$motif_cutoff)
)

motif_matrix <- plotEnrichHeatmap(
  motif_enrich,
  n = config$top_motifs,
  transpose = FALSE,
  pal = paletteContinuous(set = "whiteBlue", n = 100),
  returnMatrix = TRUE
)

motif_matrix <- motif_matrix[, intersect(config$sub_annotation_levels, colnames(motif_matrix)), drop = FALSE]

motif_matrix_scaled <- scale_matrix(motif_matrix, limit = 4)

annotate_heatmap <- function(matrix_data, highlight_genes = NULL, title = "Motif enrichment") {
  col_fun <- colorRamp2(
    breaks = c(min(matrix_data), median(matrix_data), max(matrix_data)),
    colors = c("antiquewhite", "#FDBE85", "#D94701")
  )
  highlight_idx <- integer(0)
  if (!is.null(highlight_genes)) {
    highlight_idx <- which(rownames(matrix_data) %in% highlight_genes)
  }
  row_anno <- NULL
  if (length(highlight_idx) > 0) {
    row_anno <- rowAnnotation(
      highlight = anno_mark(
        at = highlight_idx,
        labels = rownames(matrix_data)[highlight_idx],
        labels_gp = gpar(fontsize = 9)
      )
    )
  }
  ht <- Heatmap(
    matrix_data,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    name = "Norm.Enrichment",
    column_title = title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    col = col_fun,
    right_annotation = row_anno
  )
  print(ht)
}

highlight_motifs <- c("NKX2-1", "NFIC", "EGR1", "STAT1", "IRF8", "FOSL1", "RUNX1")
annotate_heatmap(motif_matrix_scaled, highlight_motifs, title = "Motif enrichment in marker peaks")

if (exists("CellFeatureHeatmap", mode = "function")) {
  CellFeatureHeatmap_fn <- get("CellFeatureHeatmap", mode = "function")
  CellFeatureHeatmap_fn(
    matrix = motif_matrix,
    feature = data.frame(cluster = rep("motif", nrow(motif_matrix)), row.names = rownames(motif_matrix)),
    cell_annotation = cell_meta %>% dplyr::select(cell, sub.annotation),
    annotation_color = list(sub.annotation = setNames(brewer.pal(6, "Set1"), config$sub_annotation_levels)),
    normalize = TRUE,
    scaling = TRUE,
    legend_title = "Enrichment",
    name = "Heatmap_Motifs_Cells"
  )
}

# ---- Homer motif summaries ----------------------------------------------
if (!is.null(config$homer_dir) && dir.exists(config$homer_dir)) {
  message("Summarising HOMER motif results from ", config$homer_dir)
  homer_types <- c("tds_high", "tds_mid2", "tds_low", "erk_low", "erk_high")
  homer_summary <- purrr::map_dfr(homer_types, function(type) {
    homer_file <- file.path(config$homer_dir, type, "HOMER-Result", "knownResults.txt")
    if (!file.exists(homer_file)) {
      warning("HOMER result not found for ", type)
      return(NULL)
    }
    tbl <- data.table::fread(homer_file, sep = "\t", check.names = FALSE) %>%
      as.data.frame()
    colnames(tbl) <- gsub("[ -]", "_", colnames(tbl))
    if (!"q-value_(Benjamini)" %in% colnames(tbl)) {
      warning("HOMER result missing q-value column for ", type)
      return(NULL)
    }
    tbl$q_value <- as.numeric(tbl$`q-value_(Benjamini)`)
    tbl$motif_name <- stringr::str_split(tbl$Motif_Name, "/", simplify = TRUE)[, 1]
    tbl <- tbl %>%
      dplyr::filter(!is.na(q_value) & q_value < 0.001) %>%
      dplyr::mutate(type = type) %>%
      dplyr::select(type, motif_name, P_value, q_value)
    tbl
  })
  if (nrow(homer_summary) > 0) {
    readr::write_csv(homer_summary, file.path(config$data_dir, "HOMER_motif_summary.csv"))
  }
}

# ---- Save project -------------------------------------------------------
message("Saving ArchR project (motif annotations cached)")
saveArchRProject(ArchRProj = proj, load = FALSE)

message("Motif analysis pipeline completed at ", Sys.time())
