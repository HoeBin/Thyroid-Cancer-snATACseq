# Thyroid Cancer snATAC-seq Analysis

Thyroid cancer is a highly heterogeneous disease whose phenotype shifts with each mutational background. While transcriptomic studies have been plentiful, the chromatin landscape that orchestrates these states remains less explored. Here, single-nucleus ATAC-seq profiles from nine tumors spanning follicular (FTC), papillary (PTC), and anaplastic (ATC) thyroid cancers resolve malignant epithelial heterogeneity by leveraging curated thyroid cancer gene signatures. Cell type–specific transcription factor programs uncover regulatory factors that define lineage and aggressiveness, while enhancer-driven networks highlight state-dependent regulatory outputs. By integrating external DNA methylation resources, we further link CpG accessibility to methylation status, demonstrating concordant epigenetic wiring across modalities. Promoter–enhancer coupling and predicted Hi-C contacts reveal subtype-biased 3D chromatin organization, providing a framework to interpret epigenetic regulation and to prioritize subtype-specific therapeutic strategies.

![Summary](https://github.com/HoeBin/Thyroid-Cancer-snATACseq/blob/main/Fig/Summary.png)

## Organized Pipeline Overview

This folder contains a streamlined ArchR-based workflow for the Thyroid ATAC project. Each script focuses on a single stage, keeps only the core logic, and includes step-by-step comments so you can rerun or adapt individual analyses easily.

## Prerequisites

- R with ArchR ecosystem (ArchR, Signac, Seurat, monocle3, etc.).
- Access to the original renv environments referenced in the scripts (`/home/Data_Drive_8TB/ghlqls/Renv/...`).
- Project input data on the same paths used inside each script (fragment files, GEO methylation tables, etc.).
- Enough CPU and RAM for ArchR to parallelize (most scripts call `addArchRThreads(threads = 16)`).

> **Tip:** Each script hard-codes the project/data roots at the top. Update those paths if your environment differs.

## Script Details

- **01_archr_total_qc_global.R** – Sets up the ArchR project, converts fragment files into Arrow files, scores doublets, and exports QC plots summarizing fragment counts and TSS enrichment. Use this to regenerate the master ArchR project.
- **02_archr_epi_subset.R** – Loads the global project, filters to epithelial cells (excluding PDTC), recomputes Iterative LSI + Harmony embeddings, clusters cells, and saves UMAP/TSNE PDFs. This is the entry point for all epithelial downstream analyses.
- **03_archr_annotation.R** – Focuses on per-sample annotations: exports grouped cell counts, UMAPs colored by metadata, and imputed module scores (TDS/ERK/BRS) per sample and group for quick QC or reporting.
- **04_archr_epi_markers_peakcalling.R** – Identifies gene-score markers for each epithelial subtype, compares them with external/internal RNA DEG lists, and writes convenient overlap summaries and ready-to-plot heatmap inputs.
- **05_archr_epi_markerpeak.R** – Works at the peak level, generating browser tracks for driver genes, computing marker peaks across multiple groupings, plotting heatmaps/MA/volcano plots, and saving pairwise differential results.
- **06_tf_motif.R** – Adds motif and TFBS annotations, runs motif enrichment on the stored marker peaks, exports enrichment heatmaps, and correlates gene scores with motif deviations to suggest candidate regulators.
- **07_pseudo_time.R** – Bridges Seurat RNA and ATAC objects: builds transfer anchors, co-embeds both modalities, adds module scores, converts to a Monocle3 CDS, and plots pseudotime trajectories with the saved trajectory object.
- **08_dna_methylation.R** – Processes the GSE97466 methylation dataset (metadata, beta matrix, liftOver), aggregates beta values per ArchR differential peak, summarizes cancer-type differences, and compares with ATAC marker peaks.
- **09_co_acc.R** – Computes co-accessibility directly from the epithelial ArchR project, annotates promoter/distal pairs, and exports full and filtered tables to support peak–gene linkage analyses.

## Running the Pipeline

1. Launch an R session with access to the same renv environment as the original project.
2. Execute the scripts in numeric order (`Rscript Organized_Pipeline/01_archr_total_qc_global.R`, etc.). Each script handles its own working directories and saves files beneath the project root.
3. Inspect the PDFs/CSVs generated in the subfolders (`Annotation`, `Markers`, `MarkerPeaks`, `Motif`, `Pseudotime`, `DAR_Methylation`, `CoAccessibility`).

With this organization you can re-run the entire workflow or pick individual stages (e.g., re-run motif enrichment only) without wading through large, monolithic notebooks.

## Citation

If you reuse or adapt this pipeline, please cite the Thyroid Cancer snATAC-seq study along with ArchR and any other core packages you rely on. Replace the placeholders below with the final publication or preprint details once available:

```
Author List. Title of the Thyroid Cancer snATAC-seq study. Journal/Preprint Server, Year.
Granja et al. ArchR: An integrative framework for single-cell chromatin accessibility analysis. bioRxiv 2020.
```
