#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(harmony)
  library(ggplot2)
  library(parallel)
  library(TeachingDemos)
  library(dplyr)
  library(readr)
})

# ---- Configuration ------------------------------------------------------
config <- list(
  project_root = getwd(),
  renv_activate = file.path("renv", "activate.R"),
  global_project_dir = "/path/to/Archr_TotalOutput_Global",
  epi_output_name = file.path("Archr_Result", "Archr_TotalOutput_Global_Epi"),
  genome = "hg38",
  threads = min(16L, max(1L, parallel::detectCores(logical = TRUE))),
  include_groups = c("Epithelial"),
  exclude_samples = c("PDTC"),
  log_subdir = "logs",
  geneset_csv = "/path/to/RNAseq_Analysis_Data/score_geneset.csv",
  filter_iterations = list(
    default = list(iterations = 2, var_features = 25000, dims = 1:30),
    recompute = list(iterations = 4, var_features = 15000, dims = 1:25)
  ),
  umap_defaults = list(n_neighbors = 30, min_dist = 0.5, metric = "cosine"),
  recompute_embeddings = FALSE
)

config$epi_output_dir <- normalizePath(
  file.path(config$project_root, config$epi_output_name),
  mustWork = FALSE
)
config$log_path <- normalizePath(
  file.path(config$epi_output_dir, config$log_subdir, "archr_epi_subset.log"),
  mustWork = FALSE
)
config$renv_activate <- file.path(config$project_root, config$renv_activate)

if (!dir.exists(config$global_project_dir)) {
  stop(
    "Global ArchR project not found. Update config$global_project_dir. Current value: ",
    config$global_project_dir
  )
}
config$global_project_dir <- normalizePath(config$global_project_dir, mustWork = TRUE)

message("Epithelial ArchR subset will be written to: ", config$epi_output_dir)

if (file.exists(config$renv_activate)) {
  message("Activating renv from ", config$renv_activate)
  source(config$renv_activate)
}

addArchRGenome(config$genome)
addArchRThreads(threads = config$threads)


dir.create(config$epi_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(config$log_path), recursive = TRUE, showWarnings = FALSE)

original_wd <- getwd()
setwd(config$epi_output_dir)
on.exit(setwd(original_wd), add = TRUE)

txtStart(config$log_path)
on.exit(try(txtStop(), silent = TRUE), add = TRUE)

message("Loading global ArchR project from ", config$global_project_dir)
proj_global <- loadArchRProject(
  path = config$global_project_dir,
  force = FALSE,
  showLogo = FALSE
)

required_cols <- c("predictedGroup_global", "Sample")
missing_cols <- setdiff(required_cols, colnames(proj_global@cellColData))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns in global project: ", paste(missing_cols, collapse = ", "))
}

cells_epi <- proj_global$cellNames[
  proj_global$predictedGroup_global %in% config$include_groups &
    !proj_global$Sample %in% config$exclude_samples
]

if (length(cells_epi) == 0) {
  stop("No cells matched the epithelial subset criteria. Check config$include_groups and exclude_samples.")
}

proj_epi <- proj_global[cells_epi, ]
rm(proj_global)

message("Subset project contains ", nCells(proj_epi), " cells across ", length(unique(proj_epi$Sample)), " samples")

saveArchRProject(
  ArchRProj = proj_epi,
  outputDirectory = config$epi_output_dir,
  load = FALSE
)

proj_epi <- loadArchRProject(
  force = FALSE,
  showLogo = FALSE
)

proj_epi <- addImputeWeights(proj_epi)

message("Running dimensionality reduction (Iterative LSI)")
proj_epi <- addIterativeLSI(
  ArchRProj = proj_epi,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = config$filter_iterations$default$iterations,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = config$filter_iterations$default$var_features,
  dimsToUse = config$filter_iterations$default$dims,
  force = TRUE
)

proj_epi <- addHarmony(
  ArchRProj = proj_epi,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

proj_epi <- addClusters(
  input = proj_epi,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)

message("TSNE/UMAP embeddings (IterativeLSI)")
proj_epi <- addTSNE(
  ArchRProj = proj_epi,
  reducedDims = "IterativeLSI",
  name = "TSNE",
  perplexity = 30,
  force = TRUE
)

plotPDF(
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Sample", embedding = "TSNE"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Clusters", embedding = "TSNE"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "predictedGroup_global", embedding = "TSNE"),
  name = "Plot-TSNE-Sample-Clusters.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

proj_epi <- addUMAP(
  ArchRProj = proj_epi,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = config$umap_defaults$n_neighbors,
  minDist = config$umap_defaults$min_dist,
  metric = config$umap_defaults$metric,
  force = TRUE
)

plotPDF(
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Sample", embedding = "UMAP"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Clusters", embedding = "UMAP"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "predictedGroup_global", embedding = "UMAP"),
  name = "Plot-UMAP-Sample-Clusters.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

plotPDF(
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "group", embedding = "UMAP"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "mut", embedding = "UMAP"),
  name = "Plot-UMAP-group-mut.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

message("TSNE/UMAP embeddings (Harmony)")
proj_epi <- addTSNE(
  ArchRProj = proj_epi,
  reducedDims = "Harmony",
  name = "TSNEHarmony",
  perplexity = 30,
  force = TRUE
)

plotPDF(
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "predictedGroup_global", embedding = "TSNEHarmony"),
  name = "Plot-TSNEHarmony-Sample-Clusters.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

proj_epi <- addUMAP(
  ArchRProj = proj_epi,
  reducedDims = "Harmony",
  name = "UMAPHarmony",
  nNeighbors = config$umap_defaults$n_neighbors,
  minDist = config$umap_defaults$min_dist,
  metric = config$umap_defaults$metric,
  force = TRUE
)

plotPDF(
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony"),
  plotEmbedding(proj_epi, colorBy = "cellColData", name = "predictedGroup_global", embedding = "UMAPHarmony"),
  name = "Plot-UMAPHarmony-Sample-Clusters.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveArchRProject(ArchRProj = proj_epi, load = FALSE)

message("Identifying marker genes at default clustering")
markers_clusters <- getMarkerFeatures(
  ArchRProj = proj_epi,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markers_clusters, file.path(config$epi_output_dir, "markersGS_Clusters.RDS"))

heatmap_clusters <- markerHeatmap(
  seMarker = markers_clusters,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = TRUE,
  nLabel = 3,
  clusterCols = FALSE
)

ComplexHeatmap::draw(heatmap_clusters, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(
  heatmap_clusters,
  name = "GeneScores-Marker-Heatmap-Clusters.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 8,
  height = 6
)

message("Identifying marker genes at resolution 0.8")
proj_epi <- addClusters(
  input = proj_epi,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  resolution = 0.8,
  name = "Clusters0.8",
  force = TRUE
)

plotPDF(
  plotEmbedding(proj_epi, name = "Clusters0.8", embedding = "UMAP"),
  name = "Plot-UMAP-Clusters0.8.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

markers_clusters08 <- getMarkerFeatures(
  ArchRProj = proj_epi,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters0.8",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markers_clusters08, file.path(config$epi_output_dir, "markersGS_Clusters0.8.RDS"))

heatmap_clusters08 <- markerHeatmap(
  seMarker = markers_clusters08,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = TRUE,
  nLabel = 3,
  clusterCols = FALSE
)

ComplexHeatmap::draw(heatmap_clusters08, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(
  heatmap_clusters08,
  name = "GeneScores-Marker-Heatmap-Clusters0.8.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 8,
  height = 6
)

saveArchRProject(ArchRProj = proj_epi, load = FALSE)

message("Loading gene-set scores from ", config$geneset_csv)
if (!file.exists(config$geneset_csv)) {
  stop("Gene set file not found: ", config$geneset_csv)
}

geneset_tbl <- readr::read_csv(config$geneset_csv, show_col_types = FALSE)

score_sets <- list(
  TDS_score = stats::na.omit(geneset_tbl$TDS_score),
  ERK_score = stats::na.omit(geneset_tbl$ERK_score),
  BRS_score = stats::na.omit(geneset_tbl$BRS_score)
)

for (score_name in names(score_sets)) {
  genes <- unique(score_sets[[score_name]])
  matched_genes <- intersect(getFeatures(proj_epi), genes)
  if (length(matched_genes) == 0) {
    warning("No genes from ", score_name, " were found in the project features")
    next
  }
  proj_epi <- addModuleScore(
    ArchRProj = proj_epi,
    features = list(matched_genes),
    useMatrix = "GeneScoreMatrix",
    name = score_name
  )
}

plotPDF(
  plotEmbedding(proj_epi, name = "TDS_score1", embedding = "UMAP", imputeWeights = getImputeWeights(proj_epi)),
  plotEmbedding(proj_epi, name = "ERK_score1", embedding = "UMAP", imputeWeights = getImputeWeights(proj_epi)),
  plotEmbedding(proj_epi, name = "BRS_score1", embedding = "UMAP", imputeWeights = getImputeWeights(proj_epi)),
  name = "Plot-UMAP-ModuleScores.pdf",
  ArchRProj = proj_epi,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveArchRProject(ArchRProj = proj_epi, load = FALSE)

if (isTRUE(config$recompute_embeddings)) {
  message("Recomputing embeddings with extended Iterative LSI parameters")
  proj_re <- loadArchRProject(force = FALSE, showLogo = FALSE)
  proj_re <- addIterativeLSI(
    ArchRProj = proj_re,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = config$filter_iterations$recompute$iterations,
    clusterParams = list(
      resolution = c(0.1, 0.2, 0.4),
      sampleCells = 50000,
      n.start = 10
    ),
    varFeatures = config$filter_iterations$recompute$var_features,
    dimsToUse = config$filter_iterations$recompute$dims,
    force = TRUE
  )
  proj_re <- addUMAP(
    ArchRProj = proj_re,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = config$umap_defaults$n_neighbors,
    minDist = config$umap_defaults$min_dist,
    metric = config$umap_defaults$metric,
    force = TRUE
  )
  plotPDF(
    plotEmbedding(proj_re, colorBy = "cellColData", name = "Sample", embedding = "UMAP"),
    plotEmbedding(proj_re, colorBy = "cellColData", name = "Clusters", embedding = "UMAP"),
    plotEmbedding(proj_re, colorBy = "cellColData", name = "group", embedding = "UMAP"),
    plotEmbedding(proj_re, colorBy = "cellColData", name = "mut", embedding = "UMAP"),
    name = "Plot-UMAP-Sample-Clusters-group-mut.pdf",
    ArchRProj = proj_re,
    addDOC = FALSE,
    width = 5,
    height = 5
  )
  plotPDF(
    plotEmbedding(proj_re, name = "TDS_score1", embedding = "UMAP", imputeWeights = getImputeWeights(proj_re)),
    plotEmbedding(proj_re, name = "ERK_score1", embedding = "UMAP", imputeWeights = getImputeWeights(proj_re)),
    plotEmbedding(proj_re, name = "BRS_score1", embedding = "UMAP", imputeWeights = getImputeWeights(proj_re)),
    name = "Plot-UMAP-ModuleScores-Extended.pdf",
    ArchRProj = proj_re,
    addDOC = FALSE,
    width = 5,
    height = 5
  )
  saveArchRProject(ArchRProj = proj_re, load = FALSE)
}

message("ArchR epithelial subset pipeline completed at ", Sys.time())
