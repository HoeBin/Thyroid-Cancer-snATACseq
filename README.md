# Thyroid Cancer snATAC-seq Analysis

Thyroid cancer exhibits substantial heterogeneity that reflects its diverse mutational landscape. To capture the regulatory complexity underlying this diversity, we profiled chromatin accessibility in nine single-nucleus ATAC-seq samples spanning follicular, papillary, and anaplastic thyroid cancer. Malignant epithelial cells were stratified with curated thyroid cancer gene signatures, uncovering subtype-specific compositions across patients. Transcription factor programs highlighted distinct regulatory circuits, enhancer activity varied by cell identity, and promoter–enhancer coupling suggested subtype-dependent 3D genome organization. Integration with external methylation data further linked CpG accessibility to methylation status, revealing epigenetic patterns that may inform subtype-focused therapeutic strategies.

![Summary](https://github.com/HoeBin/Thyroid-Cancer-snATACseq/blob/main/Fig/Summary.png)

## Getting Started

1. Clone the repository and ensure that the necessary R packages for snATAC-seq analysis are installed.
2. Update the `config` list near the top of each script to point to local data directories, output destinations, and optional resources (such as `renv` activation scripts, macs2 binaries, or HOMER outputs).
3. Run individual modules with `Rscript Code/<script-name>.R` to execute the corresponding stage of the pipeline or adapt it to new datasets.

## Pipeline Modules

- `Code/01_Global_Basic.R`
  - Builds the core single-nucleus chromatin accessibility project, performs quality control, dimensionality reduction, clustering, and RNA integration, and stores intermediate results for reuse.
  - Offers configurable sample metadata, output paths, and filtering thresholds through the `config` block.

- `Code/02_Epi_Subset.R`
  - Extracts epithelial cells, re-optimizes embeddings and clustering parameters, and computes module scores for thyroid cancer–related gene sets before saving a focused project for downstream work.

- `Code/03_Epi_Annotation.R`
  - Refines epithelial annotations via Gaussian mixture modelling on module scores, updates `sub.annotation`, and generates plots summarizing subtype distributions, sample proportions, and score variation.

- `Code/04_Epi_Markers_Peaks.R`
  - Identifies marker peaks across epithelial groups, writes browser tracks, exports pseudobulk gene-score matrices, and orchestrates peak calling. It also compares accessibility markers to RNA differential-expression sets to highlight shared drivers.

- `Code/05_Epi_MarkerPeaks.R`
  - Focuses on visualization of marker peaks through browser plots, heatmaps, pairwise comparisons, and cell-level accessibility summaries, producing publication-ready figures and BED exports.

- `Code/06_Epi_Motifs.R`
  - Performs motif enrichment analyses for pairwise contrasts and marker peaks, creates ranked motif scatter plots and enrichment heatmaps, and optionally aggregates HOMER outputs when available.

## Output Structure

Each module records logs, figures, and intermediate RDS/CSV files in a script-specific output directory (controlled via `config$output_name`). Re-running modules with updated inputs refreshes their respective outputs while preserving the organized layout of derived results.

## Citation

If these scripts contribute to your research, please cite the relevant datasets and publications and acknowledge the creators of the analysis pipeline.
