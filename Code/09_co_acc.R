suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
})

# Step 1: Configure directories and load project -----------------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))
output_dir <- file.path(epi_project_dir, "CoAccessibility")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
proj_epi <- addImputeWeights(proj_epi)

# Step 2: Prepare peak annotations --------------------------------------------------
peak_info <- data.frame(getPeakSet(proj_epi))
peak_info$peak <- paste0(peak_info$seqnames, "-", peak_info$start, "-", peak_info$end)
peak_info$idx <- seq_len(nrow(peak_info))

# Step 3: Compute co-accessibility matrix -------------------------------------------
proj_epi <- addCoAccessibility(
  ArchRProj = proj_epi,
  reducedDims = "IterativeLSI"
)
coacc <- getCoAccessibility(
  ArchRProj = proj_epi,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
coacc_df <- as.data.frame(coacc)

# Step 4: Join peak metadata to co-accessibility pairs ------------------------------
peak_meta <- peak_info %>%
  select(idx, peak, peakType, nearestGene, distToGeneStart)
coacc_df <- coacc_df %>%
  left_join(peak_meta, by = c("subjectHits" = "idx")) %>%
  rename(
    subject_peak = peak,
    subject_type = peakType,
    subject_gene = nearestGene,
    subject_dist = distToGeneStart
  ) %>%
  left_join(peak_meta, by = c("queryHits" = "idx")) %>%
  rename(
    query_peak = peak,
    query_type = peakType,
    query_gene = nearestGene,
    query_dist = distToGeneStart
  )
write_csv(coacc_df, file.path(output_dir, "coaccessibility_all_pairs.csv"))

# Step 5: Focus on promoter-distal interactions ------------------------------------
promoter_distal <- coacc_df %>%
  filter(subject_type == "Promoter", query_type == "Distal")
write_csv(promoter_distal, file.path(output_dir, "coaccessibility_promoter_distal.csv"))

# Step 6: Summaries per annotation --------------------------------------------------
summary_counts <- promoter_distal %>%
  count(subject_gene, sort = TRUE) %>%
  head(50)
write_csv(summary_counts, file.path(output_dir, "top_promoter_links.csv"))

