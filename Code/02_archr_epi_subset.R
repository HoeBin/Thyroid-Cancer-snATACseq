suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
})

# Step 1: Configure paths and load the full project -------------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
full_project_dir <- file.path(project_root, "Archr_Result", project_name)
epi_project_dir <- paste0(full_project_dir, "_Epi")
dir.create(epi_project_dir, recursive = TRUE, showWarnings = FALSE)

proj_full <- loadArchRProject(full_project_dir, showLogo = FALSE)
addArchRThreads(threads = 16)

# Step 2: Subset epithelial cells and drop PDTC sample ----------------------------
epi_cells <- which(proj_full$predictedGroup_global == "Epithelial")
non_pdtc <- which(proj_full$Sample != "PDTC")
cells_to_keep <- intersect(epi_cells, non_pdtc)
proj_epi <- proj_full[cells_to_keep, ]
saveArchRProject(ArchRProj = proj_epi, outputDirectory = epi_project_dir, load = FALSE)
rm(proj_full)

# Step 3: Reload subset project and compute embeddings -----------------------------
proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
proj_epi <- addImputeWeights(proj_epi)

proj_epi <- addIterativeLSI(
  ArchRProj = proj_epi,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
  varFeatures = 25000,
  dimsToUse = 1:30,
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

# Step 4: Generate TSNE/UMAP for raw and Harmony embeddings -----------------------
proj_epi <- addTSNE(proj_epi, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30, force = TRUE)
proj_epi <- addUMAP(proj_epi, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
proj_epi <- addTSNE(proj_epi, reducedDims = "Harmony", name = "TSNEHarmony", perplexity = 30, force = TRUE)
proj_epi <- addUMAP(proj_epi, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)

# Helper to save multiple embedding plots with consistent styling
save_embedding_panel <- function(embedding_name, file_tag) {
  p_sample <- plotEmbedding(proj_epi, colorBy = "cellColData", name = "Sample", embedding = embedding_name)
  p_cluster <- plotEmbedding(proj_epi, colorBy = "cellColData", name = "Clusters", embedding = embedding_name)
  p_annotation <- plotEmbedding(proj_epi, colorBy = "cellColData", name = "predictedGroup_global", embedding = embedding_name)
  plotPDF(p_sample, p_cluster, p_annotation,
          name = sprintf("Embedding_%s_%s.pdf", embedding_name, file_tag),
          ArchRProj = proj_epi,
          addDOC = FALSE,
          width = 5,
          height = 5)
}

save_embedding_panel("TSNE", "IterativeLSI")
save_embedding_panel("UMAP", "IterativeLSI")
save_embedding_panel("TSNEHarmony", "Harmony")
save_embedding_panel("UMAPHarmony", "Harmony")

# Step 5: Save the processed project ----------------------------------------------
saveArchRProject(ArchRProj = proj_epi, outputDirectory = epi_project_dir, load = FALSE)

