suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(minfi)
  library(cowplot)
})

# Step 1: Configure directories and load ArchR resources ---------------------------
renv_dir <- "/home/Data_Drive_8TB/ghlqls/Renv/ArchR_Integration"
source(file.path(renv_dir, "renv/activate.R"))

project_root <- "/home/Data_Drive_8TB/ghlqls/Thyroid_ATAC"
project_name <- "Archr_TotalOutput_RemoveTwo_Global_240304"
epi_project_dir <- file.path(project_root, "Archr_Result", paste0(project_name, "_Epi"))

proj_epi <- loadArchRProject(epi_project_dir, showLogo = FALSE)
orders_celltype <- c("TDS.high", "TDS.mid1", "TDS.mid2", "TDS.low", "ERK.low", "ERK.high")

# Step 2: Load GEO methylation metadata --------------------------------------------
meth_root <- "/home/Data_Drive_8TB_5/ghlqls/Thyroid_Methylation"
dataset <- "GSE97466"
dataset_dir <- file.path(meth_root, dataset)
dir.create(file.path(dataset_dir, "DAR_Methylation"), recursive = TRUE, showWarnings = FALSE)

sample_info_raw <- read.table(file.path(dataset_dir, paste0(dataset, "_series_matrix_sub.txt")),
                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rownames(sample_info_raw) <- c("Sample_GEO", "Sample_type", "Sample_sex", "Sample_race", "Sample_age", "Sample_variant")
sample_info <- as.data.frame(t(sample_info_raw))[-1, ]
sample_info <- mutate(sample_info,
                      Sample_type = str_remove(Sample_type, ".*: "),
                      Sample_sex = str_remove(Sample_sex, ".*: "),
                      Sample_race = str_remove(Sample_race, ".*: "),
                      Sample_age = str_remove(Sample_age, ".*: "),
                      Sample_variant = str_remove(Sample_variant, ".*: "))
rownames(sample_info) <- sample_info$Sample_GEO

type_map <- list(
  ATC = "anaplastic thyroid cancer",
  FTC = c("follicullar thyroid cancer", "minimally invasive follicular carcinomas"),
  PTC = "papillary thyroid cancer",
  PDTC = "poorly differentiated thyroid carcinoma",
  NAT = "non-neoplastic adjacent tissue",
  Hurthle = "Hürthle cell carcinomas",
  benign_lesions = c("follicular adenoma", "follicular adenoma/Hürthle cell", "lymphocytic thyroiditis", "nodular goiter")
)
sample_info$Sample_type <- map_chr(sample_info$Sample_type, function(x) {
  matched <- map_lgl(type_map, ~ x %in% .x)
  hit <- names(which(matched))
  if (length(hit) == 0) NA_character_ else hit[1]
})

# Step 3: Load beta matrix and probe annotations ----------------------------------
matrix_lines <- readLines(file.path(dataset_dir, paste0(dataset, "_series_matrix.txt")))
beta_matrix <- read.table(text = paste(matrix_lines[!grepl("^!", matrix_lines)], collapse = "\n"),
                          header = TRUE, row.names = 1)

probe_info <- read.csv(list.files(dataset_dir, pattern = "\\.csv$", full.names = TRUE), skip = 7) %>%
  filter(!is.na(MAPINFO)) %>%
  mutate(
    CHR = paste0("chr", CHR),
    Start = MAPINFO,
    End = MAPINFO + 50,
    Strand = if_else(Strand == "F", "+", "-")
  ) %>%
  select(IlmnID, CHR, Start, End, Strand)
rownames(probe_info) <- probe_info$IlmnID

probe_gr <- makeGRangesFromDataFrame(probe_info, keep.extra.columns = TRUE)
probe_gr_hg38 <- liftOver(probe_gr, import.chain(file.path(meth_root, "genome", "hg19ToHg38.over.chain"))) %>%
  unlist()
beta_matrix <- beta_matrix[rownames(beta_matrix) %in% names(probe_gr_hg38), ]
probe_gr_hg38 <- probe_gr_hg38[rownames(beta_matrix)]

save(beta_matrix, sample_info, probe_gr_hg38, file = file.path(dataset_dir, "methylation.rda"))

# Step 4: Summarize methylation within ArchR peaks ---------------------------------
peak_info <- data.frame(getPeakSet(proj_epi))
peak_info$peak <- paste0(peak_info$seqnames, "-", peak_info$start, "-", peak_info$end)
peak_gr <- makeGRangesFromDataFrame(peak_info, keep.extra.columns = TRUE)

overlap_counts <- countOverlaps(peak_gr, probe_gr_hg38)
peaks_with_probes <- peak_gr[overlap_counts > 0]

betamat_peak <- matrix(NA_real_, nrow = length(peaks_with_probes), ncol = ncol(beta_matrix),
                       dimnames = list(peaks_with_probes$peak, colnames(beta_matrix)))

message("Collapsing probe-level beta values to peak-level beta estimates...")
for (i in seq_along(peaks_with_probes)) {
  probes <- names(probe_gr_hg38)[countOverlaps(probe_gr_hg38, peaks_with_probes[i]) > 0]
  if (length(probes) == 1) {
    betamat_peak[i, ] <- beta_matrix[probes, ]
  } else if (length(probes) > 1) {
    betamat_peak[i, ] <- colMeans(beta_matrix[probes, , drop = FALSE], na.rm = TRUE)
  }
  if (i %% 500 == 0) message(sprintf("Processed %s / %s peaks", i, length(peaks_with_probes)))
}

betamat_peak <- as.data.frame(t(betamat_peak))
betamat_peak <- betamat_peak %>%
  rownames_to_column("Sample_GEO")
write_csv(betamat_peak, file.path(dataset_dir, "DAR_Methylation", "betamat_peak_all.csv"))

# Step 5: Compare methylation across cancer groups --------------------------------
sample_info <- mutate(sample_info, cancer_type = factor(Sample_type, levels = c("NAT", "benign_lesions", "FTC", "PTC", "ATC", "PDTC", "Hurthle")))
shared_samples <- intersect(colnames(beta_matrix), sample_info$Sample_GEO)
beta_matrix <- beta_matrix[, shared_samples]
sample_info <- sample_info[match(shared_samples, sample_info$Sample_GEO), ]
betamat_peak_filtered <- betamat_peak %>% filter(Sample_GEO %in% shared_samples)

group_means <- betamat_peak_filtered %>%
  pivot_longer(-Sample_GEO, names_to = "peak", values_to = "beta") %>%
  left_join(sample_info[, c("Sample_GEO", "cancer_type")], by = "Sample_GEO") %>%
  group_by(cancer_type, peak) %>%
  summarise(beta = mean(beta, na.rm = TRUE), .groups = "drop")
write_csv(group_means, file.path(dataset_dir, "DAR_Methylation", "group_mean_beta.csv"))

# Step 6: Overlay ATAC marker peaks with methylation --------------------------------
markers_peaks <- readRDS(file.path(epi_project_dir, "MarkerPeaks", "markersPeaks_sub.annotation.rds"))
marker_df <- as.data.frame(getMarkers(markers_peaks, cutOff = "FDR <= 0.001 & Log2FC >= 0.5"))
marker_df$peak <- paste0(marker_df$seqnames, "-", marker_df$start, "-", marker_df$end)

peak_summary <- marker_df %>%
  group_by(group_name, peakType) %>%
  summarise(n = n(), .groups = "drop")
write_csv(peak_summary, file.path(dataset_dir, "DAR_Methylation", "marker_peak_summary.csv"))

# Step 7: Simple visualization helpers --------------------------------------------
plot_peak_proportions <- peak_summary %>%
  ggplot(aes(x = group_name, y = n, fill = peakType)) +
  geom_col(color = "black") +
  scale_x_discrete(limits = orders_celltype) +
  theme_cowplot() +
  labs(x = "", y = "Number of peaks")
ggsave(file.path(dataset_dir, "DAR_Methylation", "marker_peak_counts.pdf"), plot_peak_proportions, width = 7, height = 4)
