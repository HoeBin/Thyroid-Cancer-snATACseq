suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(cowplot)
})

# Step 1: Configure directories and load the epithelial project --------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))
output_dir <- file.path(epi_project_dir, "Annotation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
proj_epi <- addImputeWeights(proj_epi)
weight_mat <- getImputeWeights(proj_epi)

# Step 2: Summarize per-sample cell counts -----------------------------------------
meta <- as.data.frame(proj_epi@cellColData)
orders_group <- c("BRAF-PTC", "BRAF-ATC", "RAS-FTC", "RAS-ATC", "NBNR-PDTC")
orders_sample <- c("ATAC3", "ATAC4", "ATAC1", "ATAC2", "ATAC8", "ATAC9", "ATAC10", "ATAC6", "ATAC7", "PDTC")

group_counts <- meta %>%
  count(group, Sample) %>%
  mutate(
    Sample = factor(Sample, levels = orders_sample),
    group = factor(group, levels = orders_group)
  ) %>%
  arrange(group, Sample)
write_csv(group_counts, file.path(output_dir, "group_counts.csv"))

# Step 3: Plot key UMAPs -----------------------------------------------------------
plot_umap <- function(field, file_tag) {
  p <- plotEmbedding(proj_epi, colorBy = "cellColData", name = field, embedding = "UMAP") +
    theme_ArchR(legendTextSize = 8, legendPosition = "right")
  plotPDF(p, name = sprintf("UMAP_%s.pdf", file_tag), ArchRProj = proj_epi, addDOC = FALSE, width = 5, height = 5)
}

umap_fields <- c(Sample = "Sample", group = "Group", mut = "Mutation")
walk(names(umap_fields), function(field) plot_umap(field, umap_fields[[field]]))

# Step 4: Visualize module scores on embeddings ------------------------------------
module_fields <- c("TDS_score1", "ERK_score1", "BRS_score1")
module_plots <- lapply(module_fields, function(field) {
  plotEmbedding(
    ArchRProj = proj_epi,
    name = field,
    imputeWeights = weight_mat,
    embedding = "UMAP"
  ) + theme_ArchR(legendTextSize = 8, legendPosition = "right") +
    ggtitle(field)
})
plotPDF(plotlist = module_plots, name = "ModuleScore_UMAPs.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 5, height = 5)

# Step 5: Export imputed module scores per sample ----------------------------------
extract_scores <- function(field) {
  scores <- proj_epi@cellColData[[field]]
  mat <- matrix(scores, nrow = 1, dimnames = list(field, proj_epi$cellNames))
  data.frame(
    Sample = proj_epi$Sample,
    group = proj_epi$group,
    score = scores,
    imputed_score = imputeMatrix(mat, imputeWeights = weight_mat)[1, ]
  )
}

score_tables <- lapply(module_fields, extract_scores)
names(score_tables) <- module_fields

walk(names(score_tables), function(field) {
  write_csv(score_tables[[field]], file.path(output_dir, sprintf("%s_by_sample.csv", field)))
})

# Step 6: Save workspace for downstream usage --------------------------------------
saveRDS(score_tables, file.path(output_dir, "module_score_tables.rds"))
