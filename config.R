# ============================================================
# config.R — pipeline_v2 (Harmony-only architecture)
#
# Differences from pipeline/config.R:
#   - No CCA-specific parameters; Harmony is the single batch-correction
#     tool, applied at both the integration step (stage 01) and the
#     rod-subset embedding step (stage 03).
#   - CHOSEN_K = 3 (metric-selected rod partition)
#   - Cell-type annotation is auto-derived from module scores in stage 02
#     rather than a hard-coded cluster→type map (cluster IDs from
#     Harmony-integrated FindClusters are not stable across pipeline
#     versions, so a marker-based assignment is more portable).
#
# Source this file at the top of every stage script.
# ============================================================

# ---- Paths (all relative to pipeline_v2/) ----
PIPELINE_DIR   <- "."
RAW_DATA_BASE  <- file.path(PIPELINE_DIR, "..")

RDS_DIR    <- file.path(PIPELINE_DIR, "output", "rds")
TABLE_DIR  <- file.path(PIPELINE_DIR, "output", "tables")
EXPORT_DIR <- file.path(PIPELINE_DIR, "output", "exports")
FIG_QC     <- file.path(PIPELINE_DIR, "output", "figures", "fig1_qc")
FIG_CT     <- file.path(PIPELINE_DIR, "output", "figures", "fig2_celltypes")
FIG_BENCH  <- file.path(PIPELINE_DIR, "output", "figures", "fig3_benchmark")
FIG_SUB    <- file.path(PIPELINE_DIR, "output", "figures", "fig4_subclusters")
FIG_PUB    <- file.path(PIPELINE_DIR, "output", "figures", "publication")

for (d in c(RDS_DIR, TABLE_DIR, EXPORT_DIR,
            FIG_QC, FIG_CT, FIG_BENCH, FIG_SUB, FIG_PUB)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Sample metadata (8 datasets total; S7/S8 excluded) ----
SAMPLES <- tibble::tribble(
  ~sample, ~dir,                                                         ~time_label, ~animal,  ~condition,           ~source,
  "S1",    file.path(RAW_DATA_BASE, "P15_GSM8133387"),                   "P15",       "Mouse",  "P015 B6C3 WT1",      "YuChen Lab",
  "S2",    file.path(RAW_DATA_BASE, "P15_GSM8133389"),                   "P15",       "Mouse",  "P015 B6C3 WT2",      "YuChen Lab",
  "S3",    file.path(RAW_DATA_BASE, "P49_GSM7734028"),                   "P49",       "Mouse",  "P049 C57BL/6J",      "University Oxford",
  "S4",    file.path(RAW_DATA_BASE, "P75_GSM7720305"),                   "P75",       "Mouse",  "P075 C57BL/6J",      "Johns Hopkins University",
  "S5",    file.path(RAW_DATA_BASE, "P75_GSM7720306"),                   "P75",       "Mouse",  "P075 C57BL/6J",      "Johns Hopkins University",
  "S6",    file.path(RAW_DATA_BASE, "P91_GSM6205478_genes_no_features"), "P91",       "Mouse",  "P091 C57BL/6J",      "Peking Univ First Hospital",
  "S7",    file.path(RAW_DATA_BASE, "P75_Cailab"),                       "P75",       "Mouse",  "P075 WT",            "Cai Lab",
  "S8",    file.path(RAW_DATA_BASE, "P225_Cailab"),                      "P225",      "Mouse",  "P225 WT",            "Cai Lab"
)

EXCLUDE_SAMPLES <- c("S7", "S8")   # Unpublished; consistent with pipeline v1.

# ---- QC thresholds (match v1) ----
QC_MIN_COUNTS     <- 400
QC_MAX_COUNTS     <- 80000
QC_MIN_FEATURES   <- 400
QC_MAX_FEATURES   <- 8000
QC_MAX_PCT_MT     <- 18
QC_MAX_PCT_HB     <- 0.30
QC_GENE_MIN_CELLS <- 5

# ---- Integration / Harmony parameters ----
N_PCS_INTEGRATION   <- 100
N_DIMS_NEIGHBORS    <- 15
CLUSTER_RESOLUTIONS <- c(0.2, 0.4, 0.6, 0.8, 1.0)
DEFAULT_RESOLUTION  <- 0.6
HARMONY_THETA       <- 2     # Harmony default
HARMONY_LAMBDA      <- 1     # Harmony default
HARMONY_MAX_ITER    <- 10    # Harmony default

# ---- Clustering benchmarking (rod subset, stage 03) ----
K_RANGE          <- 2:9       # k=1 is degenerate; not benchmarked
N_PCS_CLUSTERING <- 50
SCDHA_N_GENES    <- 3000
KNN_K            <- 20
SNN_PRUNE        <- 1/15
MAX_SIL_CELLS    <- 10000
JACCARD_BOOTSTRAPS <- 50      # bootstrap re-clusterings for the RIGOROUS stability
                              # stages (04b/04c). Definitive stability is reported there.
BENCH_JACCARD_BOOTSTRAPS <- 3 # cheaper Jaccard for the per-(method,k) benchmark heatmap
                              # (stages 03/04). The heatmap Jaccard is illustrative; a high
                              # value here is prohibitively slow because graph methods do a
                              # 60-step resolution search inside every bootstrap.
SEED             <- 1

# ---- Stability / robustness (stages 04b, 04c) ----
STABILITY_K_RANGE <- c(2, 3, 4)        # k values profiled for bootstrap stability (04b)
SCDHA_SEED_LIST   <- c(1, 2, 3, 7, 42) # scDHA training seeds compared in 04c

METHODS <- c("PCA_KMeans", "PCA_GMM", "scDHA_KMeans",
             "PCA_LL",     "scDHA_LL",
             "PCA_Leiden", "scDHA_Leiden")

# Final selected configuration for downstream stages.
# scDHA + Louvain @ k=3 is what the methods-blind metric panel selects on the
# seed=1 embedding (sizes 30,417 / 6,887 / 582). IMPORTANT — stages 04b/04c
# show this 3-way split is NOT reproducible across scDHA training seeds
# (pairwise ARI well below 1; the dominant ~80% rod cluster is the only
# seed-stable structure, while the secondary/minor clusters reshuffle). The
# minor cluster carries a ribosomal/translation signature characteristic of a
# technical artifact. Downstream stages therefore treat k=3 as the metric-
# selected partition for ONE embedding, reported together with its (poor)
# cross-seed reproducibility rather than as a robust subtype.
# See 04b_chosen_partition_stability.R and 04c_seed_robustness.R.
CHOSEN_METHOD <- "scDHA_LL"
CHOSEN_K      <- 3

# ---- Phase 1 (positive control) benchmark settings ----
# Phase 1 sweeps a wider k-range because the full retinal dataset has 8
# known cell types — we want to see the methods recover that structure.
POC_K_RANGE <- 2:12

# ---- Cone contamination filter (stage 03) ----
CONE_PANEL <- c("Opn1sw","Opn1mw","Pde6h","Pde6c","Gngt2","Gnat2",
                "Arr3","Pcp2","Gng13","Cngb3","Cnga3","Crx","Tulp1",
                "Slc24a2","Guca1a","Guca1b")
ROD_PANEL  <- c("Rho","Nrl","Nr2e3","Gnat1","Pde6a","Pde6b","Pde6g",
                "Gngt1","Rcvrn","Sag","Cnga1","Cngb1","Slc24a1","Reep6",
                "Rom1","Prph2","Aipl1")
CONE_QUANTILE_CUT <- 0.95
CONE_DELTA_CUT    <- 0.10

# ---- Cell-type marker panels for stage 02 auto-annotation ----
# Each entry: cell-type name -> vector of canonical marker genes.
# AddModuleScore is run for each panel; for each cluster the panel with
# the highest mean module score is assigned as the cell type, provided
# the score exceeds CELLTYPE_MIN_SCORE. Otherwise "Unassigned".
CELLTYPE_PANELS <- list(
  "Rods"        = c("Rho","Nrl","Nr2e3","Gnat1","Pde6a","Pde6b","Pde6g",
                    "Sag","Cnga1","Cngb1","Reep6"),
  "Cones"       = c("Opn1sw","Opn1mw","Pde6h","Pde6c","Arr3","Gngt2",
                    "Gnat2","Cngb3","Cnga3","Pcp2"),
  "BC"          = c("Vsx2","Otx2","Grik1","Scgn","Prkca","Grm6","Vsx1"),
  "Muller Glia" = c("Vim","Apoe","Glul","Rlbp1","Aqp4","Slc1a3"),
  "RGC/Neuron"  = c("Slc17a6","Pou4f1","Pou4f2","Rbpms","Nefl","Snhg11",
                    "Thy1"),
  "AC"          = c("Slc6a9","Pax6","Gad1","Gad2","Slc32a1","C1ql2"),
  "RPE"         = c("Rpe65","Mlana"),
  "Microglia"   = c("C1qa","Cx3cr1","Tmem119","Csf1r")
)
CELLTYPE_MIN_SCORE <- 0.05    # min module score above which assignment is trusted

# ---- Marker gene lists for figures ----
DOTPLOT_MARKERS <- c("Rho","Nrl","Vim","Apoe","Vsx2","Isl1","Snhg11",
                     "Opn1sw","Pde6h","Rpe65","Pax6","Slc32a1","Tmem119")

# ---- Color palettes ----
CELLTYPE_COLORS <- c(
  "Rods"        = "#E41A1C", "Cones"       = "#377EB8",
  "BC"          = "#4DAF4A", "Muller Glia" = "#984EA3",
  "RGC/Neuron"  = "#FF7F00", "AC"          = "#A65628",
  "RPE"         = "#F781BF", "Microglia"   = "#999999",
  "Unassigned"  = "#CCCCCC"
)

METHOD_COLORS <- c(
  "PCA_KMeans"   = "#1B9E77", "PCA_GMM"      = "#D95F02",
  "scDHA_KMeans" = "#7570B3", "scDHA_LL"     = "#E7298A",
  "PCA_LL"       = "#66A61E", "scDHA_Leiden" = "#A6761D",
  "PCA_Leiden"   = "#666666"
)

CLUSTER_COLORS <- c(
  "1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
  "4" = "#984EA3", "5" = "#FF7F00", "6" = "#A65628",
  "7" = "#F781BF", "8" = "#999999", "9" = "#66C2A5"
)

TIMEPOINT_COLORS <- c(
  "P15" = "#1B9E77", "P49" = "#D95F02",
  "P75" = "#7570B3", "P91" = "#E7298A", "P225" = "#66A61E"
)

# ---- Figure settings ----
FIG_DPI       <- 600
PUB_BASE_SIZE <- 10

# ---- Parallel processing ----
N_CORES <- max(1L, parallel::detectCores() - 1L)
