# scClustBench

**scClustBench** is a label-free, multi-metric pipeline for benchmarking
single-cell RNA-seq clustering and testing whether cell-type / sub-state
structure is present and reproducible. It evaluates seven clustering configurations (PCA / scDHA ×
K-means / GMM / Louvain / Leiden) across a range of cluster counts on five
non-redundant quality metrics, and distinguishes two kinds of stability:
**within-embedding** (bootstrap resampling of cells) and **across-embedding**
(re-training the scDHA representation under independent seeds).

Demonstrated on an integrated six-dataset mouse retinal atlas: whole-dataset
cell-type recovery, a within-cell-type positive control (bipolar cells), a
within-cell-type negative case (rod photoreceptors), and synthetic
(Splatter) ground-truth validation.

## Reproducibility

The pipeline runs end-to-end from raw 10x Cell Ranger inputs with fixed random
seeds (`SEED = 1` in `config.R`). There are no one-off scripts — every reported
number and figure is produced by a numbered stage.

```
Rscript run_pipeline.R              # run all stages in order
Rscript run_pipeline.R 06           # rebuild figures from cached tables
Rscript run_pipeline.R 04b 04c      # rerun stability only (uses stage-04 cache)
```

| Stage | Script | Purpose |
|-------|--------|---------|
| 01 | `01_qc_harmony_integration.R` | QC, scDblFinder doublet removal, SCTransform, Harmony integration |
| 02 | `02_cell_type_annotation.R` | Module-score annotation of retinal cell types; rod subset |
| 03 | `03_phase1_full_dataset_benchmark.R` | Whole-dataset benchmark: 7 configs × k, 5 metrics + ARI/NMI |
| 04 | `04_phase2_rod_benchmark.R` | Rod-subset benchmark: 7 configs × k, 5 metrics |
| 04b | `04b_chosen_partition_stability.R` | Within-embedding bootstrap stability (per-cluster Jaccard + ARI) |
| 04c | `04c_seed_robustness.R` | Across-embedding stability: scDHA seed sweep |
| 04d | `04d_stability_figures.R` | Stability summary figure |
| 05 | `05_chosen_partition_analysis.R` | DEGs + downstream analysis of the selected partition |
| 06 | `06_publication_figures.R` | Benchmark + evidence figures |
| 07 / 07b | `07_positive_control_bc.R`, `07b_bc_extended_ksweep.R` | Bipolar positive control + extended k-sweep |
| 08 / 08b | `08_discovery_comparison_figure.R`, `08b_bc_positive_control_figure.R` | Cross-cell-type stability comparison + bipolar detail |
| 09–11 | `09_synthetic_validation.R`, `11_synthetic_umaps.R`, `10_synthetic_figure.R` | Splatter ground-truth validation |

Shared helpers live in `utils.R`; all paths and parameters in `config.R`.

## Stability tests

The pipeline reports two forms of reproducibility separately, because they can
disagree. A partition may be stable to cell resampling within a fixed embedding
(high bootstrap Jaccard) yet not reproducible when the scDHA embedding is
retrained under a different seed. Stages `04b` (within-embedding) and `04c`
(across-embedding) quantify both; scDHA is non-deterministic across runs even at
a fixed seed, so multi-seed sensitivity is reported rather than a single run.

## Environment

R 4.5.3 with: Seurat v5 (5.4.0), scDHA (1.2.3, torch backend), harmony (2.0.2),
scDblFinder (1.24.10), ROGUE (1.0), mclust (6.1.2), aricode (1.0.3), splatter,
igraph, FNN, Matrix, cluster, uwot, dplyr, tidyr, ggplot2, patchwork.

## Data access

Public mouse-retina 10x datasets (GEO): GSM8133387, GSM8133389 (P15),
GSM7734028 (P49), GSM7720305, GSM7720306 (P75), GSM6205478 (P91). Place each
Cell Ranger output folder one level above the pipeline directory (see `SAMPLES`
in `config.R`).

## Layout

```
.
├── *.R                  numbered pipeline stages
├── config.R / utils.R   shared config + helpers
└── output/tables/       CSV result tables (large .rds and latent CSVs gitignored)
```
