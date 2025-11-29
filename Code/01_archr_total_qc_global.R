suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ggplot2)
})

# Step 1: Configure project- and data-specific paths ---------------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
fragment_root <- "/home/Data_Drive_8TB/kim3712/atac_project"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
output_dir <- file.path(project_root, "Archr_Result", project_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Step 2: Build the fragment metadata table -----------------------------------------
raw_samples <- c(
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
)

sample_table <- tibble(raw_name = raw_samples) %>%
  mutate(
    stripped = str_remove(raw_name, "\\(.*\\)$"),
    sample_id = case_when(
      str_detect(raw_name, "CR_22_12889") ~ "Normal",
      str_detect(raw_name, "CR_22_12890") ~ "PDTC",
      TRUE ~ str_split(stripped, "_", simplify = TRUE)[, 1]
    ),
    fragments = if_else(
      str_detect(raw_name, "CR_22_"),
      file.path(fragment_root, raw_name, "outs", "fragments.tsv.gz"),
      file.path(fragment_root, raw_name, "Result", "fragments.tsv.gz")
    )
  )
stopifnot(!any(duplicated(sample_table$sample_id)))
input_files <- setNames(sample_table$fragments, sample_table$sample_id)

# Step 3: Generate Arrow files with QC-aware filters --------------------------------
addArchRGenome("hg38")
addArchRThreads(threads = 16)
arrow_files <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  filterTSS = 4,
  filterFrags = 1000,
  cleanTmp = TRUE,
  threads = 16
)

# Step 4: Score potential doublets ---------------------------------------------------
addDoubletScores(
  input = arrow_files,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 4
)

# Step 5: Create the ArchR project ---------------------------------------------------
proj_total <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = output_dir,
  copyArrows = TRUE
)

# Step 6: Plot fragment and TSS-based QC metrics -------------------------------------
df_qc <- getCellColData(proj_total, select = c("log10(nFrags)", "TSSEnrichment"))
tss_plot <- ggPoint(
  x = df_qc[, 1],
  y = df_qc[, 2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df_qc[, 1], 0.99)),
  ylim = c(0, quantile(df_qc[, 2], 0.99))
) + geom_hline(yintercept = 4, linetype = "dashed") +
  geom_vline(xintercept = 2, linetype = "dashed")
plotPDF(tss_plot, name = "QC_TSS_vs_Frags.pdf", ArchRProj = proj_total, addDOC = FALSE)

# Step 7: Visualize sample-level QC distributions -----------------------------------
qc_ridge <- plotGroups(
  ArchRProj = proj_total,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)
qc_violin <- plotGroups(
  ArchRProj = proj_total,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(qc_ridge, qc_violin, name = "QC_Sample_Distributions.pdf", ArchRProj = proj_total, addDOC = FALSE)

# Step 8: Save project metadata ------------------------------------------------------
saveArchRProject(ArchRProj = proj_total, outputDirectory = output_dir, load = FALSE)
