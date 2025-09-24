#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(mclust)
  library(parallel)
  library(TeachingDemos)
})

# ---- Configuration ------------------------------------------------------
config <- list(
  project_root = getwd(),
  renv_activate = file.path("renv", "activate.R"),
  epi_project_dir = "/path/to/Archr_TotalOutput_Global_Epi",
  output_name = file.path("Archr_Result", "Archr_Epi_Annotation"),
  genome = "hg38",
  threads = min(16L, max(1L, parallel::detectCores(logical = TRUE))),
  geneset_csv = "/path/to/RNAseq_Analysis_Data/score_geneset.csv",
  module_scores = list(
    TDS = list(column = "TDS_score1", geneset = "TDS_score", name = "TDS_score"),
    ERK = list(column = "ERK_score1", geneset = "ERK_score", name = "ERK_score"),
    BRS = list(column = "BRS_score1", geneset = "BRS_score", name = "BRS_score")
  ),
  annotation_levels = c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high"),
  annotation_colors = c(
    TDS.high = "green4",
    TDS.mid1 = "darkorange",
    TDS.mid2 = "gold",
    TDS.low = "mediumorchid",
    ERK.low = "royalblue",
    ERK.high = "firebrick"
  ),
  group_colors = c(
    "BRAF-PTC" = "darkgreen",
    "BRAF-ATC" = "red",
    "RAS-FTC" = "orange",
    "RAS-ATC" = "purple",
    "NBNR-PDTC" = "black",
    "Normal" = "grey50"
  ),
  annotation_map = c(
    `1` = "TDS.high",
    `2` = "TDS.mid1",
    `3` = "TDS.mid2",
    `4` = "TDS.low",
    `5` = "ERK.low",
    `6` = "ERK.high"
  ),
  gmm_clusters = 6L,
  random_seed = 1234L,
  run_kmeans = FALSE,
  kmeans_clusters = 6L,
  log_subdir = "logs",
  plots_subdir = "plots"
)

config$output_dir <- normalizePath(
  file.path(config$project_root, config$output_name),
  mustWork = FALSE
)
config$plots_dir <- normalizePath(
  file.path(config$output_dir, config$plots_subdir),
  mustWork = FALSE
)
config$log_path <- normalizePath(
  file.path(config$output_dir, config$log_subdir, "archr_epi_annotation.log"),
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

message("Annotation outputs will be written to: ", config$output_dir)

if (file.exists(config$renv_activate)) {
  message("Activating renv from ", config$renv_activate)
  source(config$renv_activate)
}

addArchRGenome(config$genome)
addArchRThreads(threads = config$threads)


dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(config$log_path), recursive = TRUE, showWarnings = FALSE)
dir.create(config$plots_dir, recursive = TRUE, showWarnings = FALSE)

original_wd <- getwd()
setwd(config$output_dir)
on.exit(setwd(original_wd), add = TRUE)

txtStart(config$log_path)
on.exit(try(txtStop(), silent = TRUE), add = TRUE)

# ---- Helper functions ---------------------------------------------------
save_plot <- function(plot, filename, width = 6, height = 5) {
  filepath <- file.path(config$plots_dir, filename)
  ggsave(filepath, plot = plot, width = width, height = height, limitsize = FALSE)
  message("Saved plot: ", filepath)
}

get_score_table <- function(proj, score_column, weights) {
  raw <- proj@cellColData[[score_column]]
  score_matrix <- matrix(raw, nrow = 1)
  colnames(score_matrix) <- proj$cellNames
  imputed <- as.numeric(imputeMatrix(score_matrix, imputeWeights = weights))
  tibble::tibble(
    cell = proj$cellNames,
    score_raw = raw,
    score_imputed = imputed
  )
}

# ---- Load ArchR project -------------------------------------------------
message("Loading epithelial ArchR project from ", config$epi_project_dir)
proj <- loadArchRProject(
  path = config$epi_project_dir,
  force = FALSE,
  showLogo = FALSE
)

if (!"group" %in% colnames(proj@cellColData)) {
  stop("Project is missing 'group' metadata. Please run epi subset pipeline first.")
}

proj <- addImputeWeights(proj)
weights <- getImputeWeights(proj)

# ---- Ensure module scores ------------------------------------------------
if (!file.exists(config$geneset_csv)) {
  stop("Gene set file not found: ", config$geneset_csv)
}

message("Loading gene sets from ", config$geneset_csv)
geneset_tbl <- readr::read_csv(config$geneset_csv, show_col_types = FALSE)

for (score_name in names(config$module_scores)) {
  score_info <- config$module_scores[[score_name]]
  if (!score_info$column %in% colnames(proj@cellColData)) {
    genes <- geneset_tbl[[score_info$geneset]]
    genes <- genes[!is.na(genes)]
    matched_genes <- intersect(getFeatures(proj), genes)
    if (length(matched_genes) == 0) {
      warning("No genes found for module score ", score_name, ". Skipping.")
      next
    }
    message("Adding module score for ", score_name)
    proj <- addModuleScore(
      ArchRProj = proj,
      features = list(matched_genes),
      useMatrix = "GeneScoreMatrix",
      name = score_info$name
    )
  }
}

proj <- addImputeWeights(proj)
weights <- getImputeWeights(proj)

cell_meta <- proj@cellColData %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    Sample = factor(Sample),
    group = factor(group),
    mut = factor(mut)
  )

score_data <- cell_meta %>%
  dplyr::select(cell, Sample, group, mut)

for (score_name in names(config$module_scores)) {
  score_info <- config$module_scores[[score_name]]
  column <- score_info$column
  if (!column %in% colnames(cell_meta)) {
    next
  }
  score_tbl <- get_score_table(proj, column, weights)
  colnames(score_tbl)[colnames(score_tbl) == "score_raw"] <- column
  colnames(score_tbl)[colnames(score_tbl) == "score_imputed"] <- paste0(column, "_imputed")
  score_data <- score_data %>% dplyr::left_join(score_tbl, by = "cell")
}

score_data <- score_data %>%
  dplyr::mutate(
    Sample = factor(Sample, levels = sort(unique(as.character(Sample)))),
    group = factor(group, levels = names(config$group_colors)),
    mut = factor(mut)
  )

# ---- Sample summaries ----------------------------------------------------
sample_summary <- cell_meta %>%
  dplyr::count(group, Sample, name = "n") %>%
  dplyr::arrange(group, Sample)
message("Sample composition:\n", paste(utils::capture.output(print(sample_summary)), collapse = "\n"))

# ---- UMAP visualisations -------------------------------------------------
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

p_umap_mut <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "mut",
  embedding = "UMAP"
) + theme_ArchR(legendTextSize = 8, legendPosition = "right")

save_plot((p_umap_sample | p_umap_group | p_umap_mut), "UMAP-Sample-Group-Mutation.pdf", width = 15, height = 5)

module_plots <- purrr::map(
  config$module_scores,
  function(score_info) {
    column <- score_info$column
    if (!column %in% colnames(proj@cellColData)) {
      return(NULL)
    }
    plotEmbedding(
      ArchRProj = proj,
      colorBy = "cellColData",
      name = column,
      embedding = "UMAP",
      imputeWeights = weights
    ) + theme_ArchR(legendTextSize = 8, legendPosition = "right")
  }
)
module_plots <- Filter(Negate(is.null), module_plots)

if (length(module_plots) > 0) {
  save_plot(wrap_plots(module_plots, nrow = 1), "UMAP-Module-Scores.pdf", width = 5 * length(module_plots), height = 5)
}

# ---- Score distributions -------------------------------------------------
plot_score_by_sample <- function(data, score_col, title) {
  ggplot(data, aes(x = Sample, y = .data[[score_col]], fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_jitter(aes(color = group), width = 0.2, size = 0.4, alpha = 0.4) +
    scale_fill_manual(values = config$group_colors) +
    scale_color_manual(values = config$group_colors) +
    theme_cowplot() +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    ylab(title) +
    ggtitle(title)
}

score_plot_list <- list()
if (!is.null(config$module_scores$TDS$column) && paste0(config$module_scores$TDS$column, "_imputed") %in% colnames(score_data)) {
  score_plot_list[["TDS"]] <- plot_score_by_sample(score_data, paste0(config$module_scores$TDS$column, "_imputed"), "TDS Score (imputed)")
}
if (!is.null(config$module_scores$ERK$column) && paste0(config$module_scores$ERK$column, "_imputed") %in% colnames(score_data)) {
  score_plot_list[["ERK"]] <- plot_score_by_sample(score_data, paste0(config$module_scores$ERK$column, "_imputed"), "ERK Score (imputed)")
}
if (length(score_plot_list) > 0) {
  save_plot(wrap_plots(score_plot_list, nrow = 1), "Scores-by-Sample.pdf", width = 6 * length(score_plot_list), height = 4)
}

# ---- Gaussian mixture annotation ----------------------------------------
if (!all(c(paste0(config$module_scores$TDS$column, "_imputed"), paste0(config$module_scores$ERK$column, "_imputed")) %in% colnames(score_data))) {
  stop("Imputed module scores not available for GMM analysis. Check module score configuration.")
}

gmm_input <- score_data %>%
  dplyr::select(
    cell,
    Sample,
    group,
    mut,
    TDS = .data[[paste0(config$module_scores$TDS$column, "_imputed")]],
    ERK = .data[[paste0(config$module_scores$ERK$column, "_imputed")]]
  ) %>%
  tidyr::drop_na(TDS, ERK)

set.seed(config$random_seed)
message("Running Gaussian mixture model with ", config$gmm_clusters, " components")
gmm_fit <- Mclust(gmm_input[, c("TDS", "ERK")], G = config$gmm_clusters)
gmm_input$cluster <- factor(gmm_fit$classification)

unmapped_clusters <- setdiff(levels(gmm_input$cluster), names(config$annotation_map))
if (length(unmapped_clusters) > 0) {
  warning("Clusters missing in annotation_map: ", paste(unmapped_clusters, collapse = ", "))
}

gmm_input$annotation <- dplyr::recode(as.character(gmm_input$cluster), !!!config$annotation_map)

gmm_plot <- ggplot(gmm_input, aes(x = TDS, y = ERK, color = annotation)) +
  geom_point(size = 0.4, alpha = 0.6) +
  theme_cowplot() +
  scale_color_manual(values = config$annotation_colors, drop = FALSE) +
  ggtitle("Gaussian Mixture Annotation") +
  labs(color = "Annotation")

centroids <- gmm_input %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(TDS = mean(TDS), ERK = mean(ERK), .groups = "drop")

gmm_plot <- gmm_plot +
  geom_text(
    data = centroids,
    aes(label = annotation),
    fontface = "bold",
    color = "black",
    size = 3,
    vjust = -0.5
  )

save_plot(gmm_plot, "Annotation-GMM-Scatter.pdf", width = 6, height = 5)

faceted_plot <- ggplot(gmm_input, aes(x = TDS, y = ERK, color = annotation)) +
  geom_point(size = 0.3, alpha = 0.6) +
  facet_wrap(~group) +
  theme_cowplot() +
  scale_color_manual(values = config$annotation_colors, drop = FALSE) +
  ggtitle("Annotation by Group") +
  labs(color = "Annotation")

save_plot(faceted_plot, "Annotation-GMM-Scatter-By-Group.pdf", width = 8, height = 6)

score_data <- score_data %>%
  dplyr::left_join(gmm_input %>% dplyr::select(cell, annotation), by = "cell")
score_data$annotation <- factor(score_data$annotation, levels = config$annotation_levels)

proj$sub.annotation <- score_data$annotation[match(proj$cellNames, score_data$cell)]
proj$sub.annotation <- factor(proj$sub.annotation, levels = config$annotation_levels)

# ---- Summary plots ------------------------------------------------------
annotation_counts <- score_data %>%
  tidyr::drop_na(annotation) %>%
  dplyr::count(Sample, annotation) %>%
  dplyr::mutate(Sample = factor(Sample, levels = sort(unique(Sample))))

bar_sample <- ggplot(annotation_counts, aes(x = Sample, y = n, fill = annotation)) +
  geom_col(color = "black") +
  theme_cowplot() +
  scale_fill_manual(values = config$annotation_colors, drop = FALSE) +
  labs(x = "Sample", y = "Cells", fill = "Annotation", title = "Annotation counts per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(bar_sample, "Annotation-Counts-by-Sample.pdf", width = 7, height = 5)

group_proportions <- score_data %>%
  tidyr::drop_na(annotation) %>%
  dplyr::count(group, annotation) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(prop = n / sum(n), .groups = "drop")

bar_group <- ggplot(group_proportions, aes(x = group, y = prop, fill = annotation)) +
  geom_col(color = "black", position = "stack") +
  theme_cowplot() +
  scale_fill_manual(values = config$annotation_colors, drop = FALSE) +
  labs(x = "Group", y = "Proportion", fill = "Annotation", title = "Annotation proportion by group")

save_plot(bar_group, "Annotation-Proportion-by-Group.pdf", width = 6, height = 5)

score_annotation_plot <- ggplot(score_data, aes(x = Sample, y = .data[[paste0(config$module_scores$TDS$column, "_imputed")]], fill = group)) +
  geom_violin(alpha = 0.7) +
  geom_point(aes(color = annotation), size = 0.4, position = position_jitter(width = 0.2), alpha = 0.5) +
  scale_fill_manual(values = config$group_colors) +
  scale_color_manual(values = config$annotation_colors, drop = FALSE) +
  theme_cowplot() +
  labs(x = "Sample", y = "TDS Score", color = "Annotation", fill = "Group", title = "TDS Score by sample and annotation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(score_annotation_plot, "TDS-Score-by-Annotation.pdf", width = 8, height = 5)

if (isTRUE(config$run_kmeans)) {
  message("Running k-means clustering with ", config$kmeans_clusters, " centers")
  k_input <- as.matrix(gmm_input[, c("TDS", "ERK")])
  set.seed(config$random_seed)
  km_fit <- kmeans(k_input, centers = config$kmeans_clusters)
  gmm_input$kmeans_cluster <- factor(km_fit$cluster)
  kmeans_plot <- ggplot(gmm_input, aes(x = TDS, y = ERK, color = kmeans_cluster)) +
    geom_point(size = 0.4, alpha = 0.6) +
    theme_cowplot() +
    ggtitle("k-means clusters (TDS vs ERK)") +
    labs(color = "Cluster")
  save_plot(kmeans_plot, "Annotation-KMeans-Scatter.pdf", width = 6, height = 5)
}

# ---- Save project -------------------------------------------------------
message("Saving updated ArchR project with sub.annotation column")
saveArchRProject(ArchRProj = proj, load = FALSE)

message("Annotation pipeline completed at ", Sys.time())
