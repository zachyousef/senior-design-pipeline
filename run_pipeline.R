# ============================================================
# run_pipeline.R   (pipeline_v2)
#
# Master orchestrator for the two-phase pipeline.
#
# Stages:
#   01 — QC + Harmony integration of full retinal dataset
#   02 — Automated cell-type annotation + rod isolation
#   03 — PHASE 1: positive control benchmark (full dataset, 7×11=77 partitions)
#   04 — PHASE 2: rod heterogeneity benchmark (rod subset, 7×8=56 partitions)
#  04b — Bootstrap stability of the chosen method across k (within-embedding)
#  04c — scDHA seed robustness of the chosen partition (across-embedding)
#   05 — Chosen partition analysis (Phase 2 selection: scDHA + Louvain k=3)
#   06 — Publication figures (Figure 1 + 2 + Phase 1 P1–P4 + Phase 2 P5–P7)
#
# Usage:
#   Rscript run_pipeline.R                 # run all stages in order
#   Rscript run_pipeline.R 01 02 06        # run only specified stages
#   Rscript run_pipeline.R 04b 04c         # rerun stability only (needs stage 04 cache)
#   Rscript run_pipeline.R 06              # rebuild figures from cached data
#
# Wall-time on a workstation (approx):
#   01:  ~30 min (full QC + integration)
#   02:  ~10 min (annotation + rod isolation)
#   03:  ~12 hr  (Phase 1 — long; 77 partitions × 5 metrics + ARI/NMI)
#   04:  ~4.5 hr (Phase 2 — 56 partitions × 5 metrics, bootstrap=50)
#  04b:  ~1 hr   (bootstrap stability of chosen method, k in STABILITY_K_RANGE)
#  04c:  ~20 min (scDHA seed sweep: trains SCDHA_SEED_LIST embeddings)
#   05:  ~10 min (apply chosen partition + DEGs)
#   06:  ~3 min  (figure generation from cached data)
# ============================================================

source("config.R")

args <- commandArgs(trailingOnly = TRUE)
all_stages <- c("01"  = "01_qc_harmony_integration.R",
                "02"  = "02_cell_type_annotation.R",
                "03"  = "03_phase1_full_dataset_benchmark.R",
                "04"  = "04_phase2_rod_benchmark.R",
                "04b" = "04b_chosen_partition_stability.R",
                "04c" = "04c_seed_robustness.R",
                "04d" = "04d_stability_figures.R",
                "05"  = "05_chosen_partition_analysis.R",
                "06"  = "06_publication_figures.R",
                "07"  = "07_positive_control_bc.R",
                "07b" = "07b_bc_extended_ksweep.R",
                "08b" = "08b_bc_positive_control_figure.R",
                "08"  = "08_discovery_comparison_figure.R",
                "09"  = "09_synthetic_validation.R",
                "11"  = "11_synthetic_umaps.R",
                "10"  = "10_synthetic_figure.R")
stages <- if (length(args) == 0) names(all_stages) else args
stages <- intersect(stages, names(all_stages))

T0 <- Sys.time()
for (s in stages) {
  message("\n###################################################")
  message("###  STAGE ", s, ": ", all_stages[s])
  message("###################################################\n")
  st_t0 <- Sys.time()
  source(all_stages[s], chdir = FALSE, echo = FALSE)
  message(sprintf("\n>>> Stage %s done (%s)\n", s,
                   format(Sys.time() - st_t0, digits = 3)))
  gc()
}
message(sprintf("\n>>> Pipeline_v2 total wall-time: %s",
                 format(Sys.time() - T0, digits = 3)))
