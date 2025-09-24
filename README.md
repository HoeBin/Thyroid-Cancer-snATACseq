# Thyroid Cancer snATAC-seq Pipelines

This repository contains a scripted ArchR analysis workflow for thyroid cancer single-nucleus ATAC-seq data. Each module is packaged as an executable `Rscript` with a shared configuration pattern to simplify customization and reproducibility.

![Summary](https://github.com/HoeBin/Thyroid-Cancer-snATACseq/blob/main/Fig/Summary.png)

## Getting Started

1. Clone the repository and install the required R packages (ArchR and its dependencies).
2. For each script, adjust the `config` block near the top to point at your data directories, project output locations, and optional resources (e.g. `renv` activation, macs2 binary, HOMER outputs).
3. Execute individual steps with `Rscript Code/<script-name>.R` to replicate or adapt portions of the analysis pipeline.

## Pipeline Modules

- `Code/01_Global_Basic.R`
  - Builds the global ArchR project, performs QC plotting, dimensional reduction, clustering, RNA integration, and saves intermediate projects and marker results.
  - Accepts configurable sample metadata, output paths, and filtering thresholds via the `config` list.

- `Code/02_Epi_Subset.R`
  - Subsets epithelial cells from the global project, re-clusters with tailored parameters, regenerates embeddings, and computes module scores for predefined gene sets.
  - Persists the subset project to a dedicated directory for downstream analyses.

- `Code/03_Epi_Annotation.R`
  - Derives refined epithelial annotations using Gaussian mixture modelling on TDS/ERK module scores, updates `sub.annotation`, and exports summary plots (counts, proportions, violin plots).
  - Logs outputs and warns when expected metadata or colour mappings are missing.

- `Code/04_Epi_Markers_Peaks.R`
  - Compares epithelial groups to define marker peaks, draws browser tracks, exports pseudobulk gene-score matrices, and triggers macs2-based peak calling.
  - Generates heatmaps and Venn-style overlap summaries for marker genes versus RNA DEG sets.

- `Code/05_Epi_MarkerPeaks.R`
  - Focuses on visualization of marker peaks: browser tracks by grouping, marker heatmaps, pairwise tests, and cell-level peak heatmaps when helper functions are available.
  - Writes marker peak objects, BED files for up/down peaks, and styled ArchR plots under `plots/`.

- `Code/06_Epi_Motifs.R`
  - Adds motif annotations, performs motif enrichment for pairwise marker tests and marker peaks, and produces ranked scatter plots plus enrichment heatmaps.
  - Optionally summarizes HOMER motif results if a HOMER output directory is provided.

## Output Structure

Each script writes logs, figures, and intermediate RDS/CSV exports to a script-specific output directory (configured via `config$output_name`). Re-running scripts with updated inputs will refresh the contents while preserving the organized layout.

## Citation

If you use these scripts, please cite ArchR and any accompanying datasets or publications described in your work.
