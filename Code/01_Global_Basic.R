#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(stringr)
  library(harmony)
  library(ggplot2)
  library(parallel)
  library(TeachingDemos)
  library(dplyr)
})

# ---- Configuration ------------------------------------------------------
config <- list(
  project_root = getwd(),
  input_directory = "/path/to/atac_project",
  output_directory_name = file.path("Archr_Result", "Archr_TotalOutput_Global"),
  seurat_rds = "/path/to/thyroid.all.102845cell.rds",
  genome = "hg38",
  threads = min(16L, max(1L, parallel::detectCores(logical = TRUE))),
  filter_tss = 4,
  filter_frags = 1000,
  log_subdir = "logs",
  renv_activate = file.path("renv", "activate.R")
)

config$output_directory <- normalizePath(
  file.path(config$project_root, config$output_directory_name),
  mustWork = FALSE
)
config$renv_activate <- file.path(config$project_root, config$renv_activate)
config$log_path <- normalizePath(
  file.path(config$output_directory, config$log_subdir, "archr_total_output_global.log"),
  mustWork = FALSE
)

message("ArchR outputs will be written to: ", config$output_directory)

if (file.exists(config$renv_activate)) {
  message("Activating renv from ", config$renv_activate)
  source(config$renv_activate)
}

dir.create(config$output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(config$log_path), recursive = TRUE, showWarnings = FALSE)

original_wd <- getwd()
setwd(config$output_directory)
on.exit(setwd(original_wd), add = TRUE)

txtStart(config$log_path)
on.exit(try(txtStop(), silent = TRUE), add = TRUE)

# ---- Sample metadata ----------------------------------------------------
sample_metadata <- data.frame(
  raw_sample_id = c(
    "ATAC1_TC1_9(CR_23_09370_TS_C_SSQ_1)",
    "ATAC2_TC1_10(CR_23_09371_TS_C_SSQ_1)",
    "ATAC3_TC1_3(CR_23_09372_TS_C_SSQ_1)",
    "ATAC4_TC2_4(CR_23_09373_TS_C_SSQ_1)",
    "ATAC6_TC1_6(CR_23_09375_TS_C_SSQ_1)",
    "ATAC7_TC2_2(CR_23_09376_TS_C_SSQ_1)",
    "ATAC8_TC2_5(CR_23_09377_TS_C_SSQ_1)",
    "ATAC9_TC1_8(CR_23_09378_TS_C_SSQ_1)",
    "ATAC10_TC2_14(CR_23_09379_TS_C_SSQ_1)",
    "CR_22_12890_TS_C_SSQ_1"
  ),
  alias = c(
    "ATAC1", "ATAC2", "ATAC3", "ATAC4", "ATAC6",
    "ATAC7", "ATAC8", "ATAC9", "ATAC10", "PDTC"
  ),
  fragments_subpath = c(
    rep("Result/fragments.tsv.gz", 9),
    "outs/fragments.tsv.gz"
  ),
  stringsAsFactors = FALSE
)

sample_metadata$fragments_path <- file.path(
  config$input_directory,
  sample_metadata$raw_sample_id,
  sample_metadata$fragments_subpath
)

missing_fragments <- sample_metadata$fragments_path[!file.exists(sample_metadata$fragments_path)]
if (length(missing_fragments) > 0) {
  stop(
    "Fragments files not found. Update config$input_directory or sample_metadata. Missing: ",
    paste(missing_fragments, collapse = ", ")
  )
}

input_files <- sample_metadata$fragments_path
names(input_files) <- sample_metadata$alias

# ---- ArchR project setup -------------------------------------------------
addArchRGenome(config$genome)
addArchRThreads(threads = config$threads)

arrow_files <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  cleanTmp = TRUE,
  nChunk = 1,
  threads = config$threads,
  filterTSS = config$filter_tss,
  filterFrags = config$filter_frags
)

doublet_scores <- addDoubletScores(
  input = arrow_files,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 1
)

proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = config$output_directory,
  copyArrows = TRUE
)

message("Created ArchRProject with ", nCells(proj), " cells")

# ---- QC metrics ---------------------------------------------------------
qc_metrics <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

qc_scatter <- ggPoint(
  x = qc_metrics[, 1],
  y = qc_metrics[, 2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(qc_metrics[, 1], probs = 0.99)),
  ylim = c(0, quantile(qc_metrics[, 2], probs = 0.99))
) +
  geom_hline(yintercept = config$filter_tss, lty = "dashed") +
  geom_vline(xintercept = log10(config$filter_frags), lty = "dashed")

plotPDF(qc_scatter, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

sample_tss_ridge <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)

sample_tss_violin <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

sample_frags_ridge <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
)

sample_frags_violin <- plotGroups(
  ArchRProj = proj,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(
  sample_tss_ridge,
  sample_tss_violin,
  sample_frags_ridge,
  sample_frags_violin,
  name = "QC-Sample-Statistics.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 4,
  height = 4
)

fragment_sizes <- plotFragmentSizes(ArchRProj = proj)
tsse_profile <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(
  fragment_sizes,
  tsse_profile,
  name = "QC-Sample-FragSizes-TSSProfile.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveArchRProject(ArchRProj = proj, outputDirectory = config$output_directory, load = FALSE)

# ---- Dimensional reduction & clustering ---------------------------------
proj <- filterDoublets(proj)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force = TRUE
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)

proj <- addTSNE(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "TSNE",
  perplexity = 30
)

tsne_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "TSNE"
)
tsne_clusters <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "TSNE"
)
plotPDF(
  tsne_sample,
  tsne_clusters,
  name = "Plot-TSNE-Sample-Clusters.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

umap_sample <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)
umap_clusters <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)
plotPDF(
  umap_sample,
  umap_clusters,
  name = "Plot-UMAP-Sample-Clusters.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveArchRProject(ArchRProj = proj, outputDirectory = config$output_directory, load = FALSE)

# ---- Marker discovery ---------------------------------------------------
markers_gs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markers_gs, file.path(config$output_directory, "markersGS.RDS"))

marker_heatmap <- markerHeatmap(
  seMarker = markers_gs,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = TRUE
)

ComplexHeatmap::draw(marker_heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(
  marker_heatmap,
  name = "GeneScores-Marker-Heatmap.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 8,
  height = 6
)

# ---- Sample annotations -------------------------------------------------
sample_groups <- data.frame(
  alias = c(
    "ATAC1", "ATAC2", "ATAC3", "ATAC4", "ATAC5", "ATAC6",
    "ATAC7", "ATAC8", "ATAC9", "ATAC10", "PDTC", "Normal"
  ),
  group = c(
    "BRAF-ATC", "BRAF-ATC", "BRAF-PTC", "BRAF-PTC", "RAS-ATC", "RAS-ATC",
    "RAS-ATC", "RAS-FTC", "RAS-FTC", "RAS-FTC", "NBNR-PDTC", "Normal"
  ),
  mutation = c(
    "BRAF", "BRAF", "BRAF", "BRAF", "RAS", "RAS",
    "RAS", "RAS", "RAS", "RAS", "NBNR", "Normal"
  ),
  stringsAsFactors = FALSE
)

proj$group <- sample_groups$group[match(proj$Sample, sample_groups$alias)]
proj$group <- factor(
  proj$group,
  levels = c("BRAF-PTC", "BRAF-ATC", "RAS-FTC", "RAS-ATC", "NBNR-PDTC", "Normal")
)
proj$mut <- sample_groups$mutation[match(proj$Sample, sample_groups$alias)]
proj$mut <- factor(
  proj$mut,
  levels = c("BRAF", "RAS", "NBNR", "Normal")
)

if (any(is.na(proj$group))) {
  warning("Some samples lack group annotations. Check 'sample_groups'.")
}

if (any(is.na(proj$mut))) {
  warning("Some samples lack mutation annotations. Check 'sample_groups'.")
}

cell_summary <- as.data.frame(proj@cellColData) %>%
  dplyr::group_by(group, Sample) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(group) %>%
  dplyr::mutate(
    color = dplyr::case_when(
      group == "BRAF-PTC" ~ "darkgreen",
      group == "BRAF-ATC" ~ "red",
      group == "RAS-FTC" ~ "orange",
      group == "RAS-ATC" ~ "purple",
      group == "NBNR-PDTC" ~ "black",
      TRUE ~ "grey50"
    )
  )
print(cell_summary)

umap_group <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "group",
  embedding = "UMAP"
)
umap_mut <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "mut",
  embedding = "UMAP"
)
plotPDF(
  umap_group,
  umap_mut,
  name = "Plot-UMAP-group-mut.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

proj <- addImputeWeights(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = config$output_directory, load = FALSE)

# ---- scRNA-seq integration ----------------------------------------------
scrna <- readRDS(config$seurat_rds)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = scrna,
  addToArrow = FALSE,
  force = FALSE,
  groupRNA = "global.annotation.final",
  nameCell = "predictedCell_global_Un",
  nameGroup = "predictedGroup_global_Un",
  nameScore = "predictedScore_global_Un",
  threads = config$threads
)

unconstrained_group <- plotEmbedding(
  proj,
  colorBy = "cellColData",
  name = "predictedGroup_global_Un",
  embedding = "UMAP"
)
unconstrained_score <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "predictedScore_global_Un",
  embedding = "UMAP"
)
plotPDF(
  unconstrained_group,
  unconstrained_score,
  name = "RNA-Integration-UnConstrained.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

conf_mat <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_global_Un))
cluster_assignments <- colnames(conf_mat)[apply(conf_mat, 1, which.max)]

clusters_epi <- rownames(conf_mat)[cluster_assignments == "Epithelial"]
clusters_tme <- rownames(conf_mat)[cluster_assignments %in% c("Fibroblast", "Endothelial", "Myeloid", "T_NK")]

rna_epi <- colnames(scrna)[scrna$global.annotation.final == "Epithelial"]
rna_tme <- colnames(scrna)[scrna$global.annotation.final %in% c("Fibroblast", "Endothelial", "Myeloid", "T_NK")]

group_list <- S4Vectors::SimpleList(
  Epi = S4Vectors::SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% clusters_epi],
    RNA = rna_epi
  ),
  TME = S4Vectors::SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% clusters_tme],
    RNA = rna_tme
  )
)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = scrna,
  addToArrow = TRUE,
  force = TRUE,
  groupList = group_list,
  groupRNA = "global.annotation.final",
  nameCell = "predictedCell_global",
  nameGroup = "predictedGroup_global",
  nameScore = "predictedScore_global",
  threads = config$threads
)

constrained_group <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "predictedGroup_global",
  embedding = "UMAP"
)
constrained_score <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "predictedScore_global",
  embedding = "UMAP"
)
plotPDF(
  constrained_group,
  constrained_score,
  name = "RNA-Integration-Constrained.pdf",
  ArchRProj = proj,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveArchRProject(ArchRProj = proj, load = FALSE)

message("ArchR global QC pipeline completed at ", Sys.time())
