#!/usr/bin/env Rscript

# Cleaned ArchR Global Pipeline (preserves original behavior)
# - Focus: code organization, readability, minimal functional changes
# - Original file: 1_ArchR_Total_QcArchR_Global.R
# - Notes:
#   * Centralizes configuration at the top (original absolute paths kept).
#   * Groups logic into labeled sections and helper functions.
#   * Removes unused or duplicate lines and tidy prints.
#   * Keeps thresholds, resolution, and steps identical to original.

suppressPackageStartupMessages({
  library(ArchR)
  library(stringr)
  library(harmony)
  library(ggplot2)
  options(warn = 1)
  library(parallel)
  library(TeachingDemos)
  library(dplyr)
})

rm(list = ls())

# =============================
# Configuration
# =============================
config <- list(
  archr_genome     = "hg38",
  threads          = 16,
  # Paths (kept from original; adjust as needed)
  renv_dir         = "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration/",
  base_dir         = "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC/",
  input_dir        = "/home/Data_Drive_8TB/kim3712/atac_project/",
  output_dir       = "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC/Archr_Result/Archr_TotalOutput_RemoveTwo_Global_240304/",
  # QC thresholds
  filterTSS        = 4,
  filterFrags      = 1000,
  # Dimred / clustering / embedding
  lsi_varFeatures  = 25000,
  lsi_dims         = 1:30,
  seurat_resolution= 0.5,
  umap_nNeighbors  = 30,
  umap_minDist     = 0.5,
  # scRNA reference
  scrna_rds        = "/home/Data_Drive_8TB_2/ghlqls/Thyroid_scRNA/thyroid.all.102845cell.rds",
  scrna_label_col  = "global.annotation.final",
  # Logging
  log_txt          = "./output_txt/Archr_TotalOutput_RemoveTwo_Global_240304.txt"
)

# Sample list (preserved)
data_list <- c(
  'ATAC1_TC1_9(CR_23_09370_TS_C_SSQ_1)',
  'ATAC2_TC1_10(CR_23_09371_TS_C_SSQ_1)',
  'ATAC3_TC1_3(CR_23_09372_TS_C_SSQ_1)',
  'ATAC4_TC2_4(CR_23_09373_TS_C_SSQ_1)',
  # 'ATAC5_TC1_5(CR_23_09374_TS_C_SSQ_1)',
  'ATAC6_TC1_6(CR_23_09375_TS_C_SSQ_1)',
  'ATAC7_TC2_2(CR_23_09376_TS_C_SSQ_1)',
  'ATAC8_TC2_5(CR_23_09377_TS_C_SSQ_1)',
  'ATAC9_TC1_8(CR_23_09378_TS_C_SSQ_1)',
  'ATAC10_TC2_14(CR_23_09379_TS_C_SSQ_1)',
  # 'CR_22_12889_TS_C_SSQ_1',
  'CR_22_12890_TS_C_SSQ_1'
)

# =============================
# Helpers
# =============================
ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

make_input_files <- function(data_list, directory_input) {
  sample_id <- c(); inputFiles <- c()
  for (sample in data_list) {
    if (sample == "CR_22_12889_TS_C_SSQ_1") {
      sample_id_tmp <- "Normal"
      inputFiles_tmp <- file.path(directory_input, sample, 'outs/fragments.tsv.gz')
    } else if (sample == "CR_22_12890_TS_C_SSQ_1") {
      sample_id_tmp <- "PDTC"
      inputFiles_tmp <- file.path(directory_input, sample, 'outs/fragments.tsv.gz')
    } else {
      sample_id_tmp <- gsub("\\(.*?)", "", sample)  # remove parentheses contents (preserved)
      sample_id_tmp <- str_split(sample_id_tmp, "_", simplify = TRUE)[, 1]
      inputFiles_tmp <- file.path(directory_input, sample, 'Result/fragments.tsv.gz')
    }
    sample_id <- c(sample_id, sample_id_tmp)
    inputFiles <- c(inputFiles, inputFiles_tmp)
  }
  names(inputFiles) <- sample_id
  list(sample_id = sample_id, inputFiles = inputFiles)
}

plot_qc <- function(proj) {
  df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
  p <- ggPoint(
    x = df[, 1], y = df[, 2], colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[, 1], probs = 0.99)),
    ylim = c(0, quantile(df[, 2], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 2, lty = "dashed")
  plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

  p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
  p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
  p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges")
  p4 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
  plotPDF(p1, p2, p3, p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

  pfs <- plotFragmentSizes(ArchRProj = proj)
  pts <- plotTSSEnrichment(ArchRProj = proj)
  plotPDF(pfs, pts, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

run_dimred_cluster_embed <- function(proj, cfg) {
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
    varFeatures = cfg$lsi_varFeatures,
    dimsToUse = cfg$lsi_dims,
    force = TRUE
  )

  proj <- addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")

  proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = cfg$seurat_resolution, force = TRUE)

  proj <- addTSNE(ArchRProj = proj, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
  plotPDF(p1, p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = cfg$umap_nNeighbors, minDist = cfg$umap_minDist, metric = "cosine", force = TRUE)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
  plotPDF(p3, p4, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  proj
}

find_markers <- function(proj) {
  markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
  saveRDS(markersGS, file.path(getOutputDirectory(proj), "markersGS.RDS"))
  heatmapGS <- markerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", transpose = TRUE)
  ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)
  invisible(markerList)
}

add_group_labels <- function(proj) {
  # Map samples to groups (preserved mapping)
  proj@cellColData$group <- ifelse(proj$Sample %in% c("ATAC3", "ATAC4"), "BRAF-PTC",
                             ifelse(proj$Sample %in% c("ATAC1", "ATAC2"), "BRAF-ATC",
                               ifelse(proj$Sample %in% c("ATAC8", "ATAC9", "ATAC10"), "RAS-FTC",
                                 ifelse(proj$Sample %in% c("ATAC5", "ATAC6", "ATAC7"), "RAS-ATC",
                                   ifelse(proj$Sample %in% c("PDTC"), "NBNR-PDTC",
                                     ifelse(proj$Sample %in% c("Normal"), "Normal", NA))))))
  )
  proj@cellColData$group <- factor(proj$group, levels = c("BRAF-PTC", "BRAF-ATC", "RAS-FTC", "RAS-ATC", "NBNR-PDTC", "Normal"))
  proj@cellColData$mut <- ifelse(proj$group %in% c("BRAF-ATC", "BRAF-PTC"), "BRAF",
                           ifelse(proj$group %in% c("RAS-ATC", "RAS-FTC"), "RAS",
                             ifelse(proj$group %in% c("NBNR-PDTC"), "NBNR",
                               ifelse(proj$group %in% c("Normal"), "Normal", NA))))
  )

  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "group", embedding = "UMAP")
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "mut", embedding = "UMAP")
  plotPDF(p1, p2, name = "Plot-UMAP-group-mut.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
  proj
}

integrate_scrna <- function(proj, cfg) {
  if (!file.exists(cfg$scrna_rds)) {
    message("scRNA RDS not found; skipping integration: ", cfg$scrna_rds)
    return(proj)
  }
  seRNA <- readRDS(cfg$scrna_rds)

  # Unconstrained integration
  proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    force = FALSE,
    groupRNA = cfg$scrna_label_col,
    nameCell = "predictedCell_global_Un",
    nameGroup = "predictedGroup_global_Un",
    nameScore = "predictedScore_global_Un",
    threads = config$threads
  )
  p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_global_Un", embedding = "UMAP")
  p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedScore_global_Un", embedding = "UMAP")
  plotPDF(p1, p2, name = "RNA-Integration-UnConstrained.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  # Constrained integration (Epithelial vs TME)
  cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_global_Un))
  preClust <- colnames(cM)[apply(cM, 1, which.max)]
  clustEpi <- rownames(cM)[preClust == "Epithelial"]
  clustFibro <- rownames(cM)[preClust == "Fibroblast"]
  clustEndo <- rownames(cM)[preClust == "Endothelial"]
  clustMye <- rownames(cM)[preClust == "Myeloid"]
  clustTNK <- rownames(cM)[preClust == "T_NK"]
  clustTME <- c(clustFibro, clustEndo, clustMye, clustTNK)

  rnaEpi <- colnames(seRNA)[seRNA@meta.data[[cfg$scrna_label_col]] == "Epithelial"]
  rnaTME <- colnames(seRNA)[seRNA@meta.data[[cfg$scrna_label_col]] %in% c("Fibroblast", "Endothelial", "Myeloid", "T_NK")]

  groupList <- S4Vectors::SimpleList(
    Epi = S4Vectors::SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustEpi], RNA = rnaEpi),
    TME = S4Vectors::SimpleList(ATAC = proj$cellNames[proj$Clusters %in% clustTME], RNA = rnaTME)
  )

  proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force = TRUE,
    groupList = groupList,
    groupRNA = cfg$scrna_label_col,
    nameCell = "predictedCell_global",
    nameGroup = "predictedGroup_global",
    nameScore = "predictedScore_global",
    threads = config$threads
  )
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_global", embedding = "UMAP")
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedScore_global", embedding = "UMAP")
  plotPDF(p3, p4, name = "RNA-Integration-Constrained.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  proj
}

# =============================
# Main
# =============================
setwd(config$renv_dir)
if (file.exists("./renv/activate.R")) {
  source("./renv/activate.R")
}
setwd(config$base_dir)

ensure_dir(dirname(config$log_txt))
txtStart(config$log_txt)
on.exit(txtStop(), add = TRUE)

addArchRGenome(config$archr_genome)
addArchRThreads(threads = config$threads)

ensure_dir(config$output_dir)
setwd(config$output_dir)

# Prepare inputs
inp <- make_input_files(data_list, config$input_dir)
inputFiles <- inp$inputFiles
rm(inp)

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  cleanTmp = TRUE,
  nChunk = 1,
  threads = config$threads,
  filterTSS = config$filterTSS,
  filterFrags = config$filterFrags
)

# Doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 1
)

# Project init
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = config$output_dir, copyArrows = TRUE)

# QC plots
plot_qc(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = config$output_dir, load = FALSE)

# Doublet filtering
proj <- filterDoublets(proj)

# Dimred / clustering / embeddings
proj <- run_dimred_cluster_embed(proj, config)
saveArchRProject(ArchRProj = proj, outputDirectory = config$output_dir, load = FALSE)

# Marker discovery
find_markers(proj)

# Add group labels and UMAP
proj <- add_group_labels(proj)

# Impute weights
proj <- addImputeWeights(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = config$output_dir, load = FALSE)

# scRNA integration
proj <- integrate_scrna(proj, config)

# Final save
saveArchRProject(ArchRProj = proj, load = FALSE)

message("Pipeline completed.")
