# ============================================================
# 01_qc_harmony_integration.R   (pipeline_v2)
#
# Harmony-only integration replacing CCA. Pipeline:
#   1. Load 10X data per sample, tag metadata.
#   2. Per-sample doublet removal (scDblFinder).
#   3. Cell-level QC filtering (counts/features/MT/HB).
#   4. SCTransform per sample (joined into one model via SCT v2).
#   5. PCA on the merged SCT-normalized object (N_PCS_INTEGRATION components).
#   6. **RunHarmony** with vars_use = "sample" — produces "harmony" reduction.
#      This is the SINGLE batch-correction step in the pipeline at the
#      full-dataset level. No CCA. No IntegrateLayers.
#   7. FindNeighbors/FindClusters/UMAP all on the harmony reduction.
#
# Input:  Raw 10X directories (config.R SAMPLES)
# Output: output/rds/01_integrated.rds     (Seurat object, harmony reduction)
#         output/figures/fig1_qc/          (QC violins + scatter + UMAPs)
#         output/tables/01_doublet_summary.csv
#         output/tables/01_qc_filter_summary.csv
# ============================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(sctransform)
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(Matrix)
})

options(future.globals.maxSize = 32000 * 1024^2)
set.seed(SEED)

# ============================================================
# PHASE 1 — Load 10X data and tag metadata
# ============================================================
message("=== PHASE 1: Loading 10X data ===")

samples_to_use <- SAMPLES
if (!is.null(EXCLUDE_SAMPLES) && length(EXCLUDE_SAMPLES) > 0) {
  samples_to_use <- SAMPLES[!SAMPLES$sample %in% EXCLUDE_SAMPLES, ]
  message("  Excluding samples: ", paste(EXCLUDE_SAMPLES, collapse = ", "))
}

obj_list <- list()
for (i in seq_len(nrow(samples_to_use))) {
  s <- samples_to_use[i, ]
  message("  Loading ", s$sample, " from ", s$dir)

  counts <- read10x_robust(s$dir)
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
  obj$sample     <- s$sample
  obj$animal     <- s$animal
  obj$condition  <- s$condition
  obj$source     <- s$source
  obj$time_label <- s$time_label
  obj <- RenameCells(obj, add.cell.id = s$sample)
  Idents(obj) <- "condition"
  obj_list[[s$sample]] <- obj
  message("    Cells: ", ncol(obj), " | Genes: ", nrow(obj))
}

# ============================================================
# PHASE 2 — Doublet removal (per-sample scDblFinder)
# ============================================================
message("\n=== PHASE 2: Doublet removal ===")

singlet_list   <- list()
doublet_summary <- data.frame()

for (nm in names(obj_list)) {
  message("  Processing ", nm, " ...")
  obj <- obj_list[[nm]]
  n_before <- ncol(obj)

  set.seed(SEED)
  sce <- scDblFinder(as.SingleCellExperiment(obj), dbr.sd = 0)
  obj <- as.Seurat(sce, data = NULL)
  obj <- subset(obj, subset = scDblFinder.class == "singlet")
  Idents(obj) <- "condition"

  singlet_list[[nm]] <- obj
  doublet_summary <- rbind(doublet_summary, data.frame(
    sample = nm, before = n_before, after = ncol(obj),
    removed = n_before - ncol(obj),
    pct = round(100 * (n_before - ncol(obj)) / n_before, 2)))
  message("    Kept ", ncol(obj), " / ", n_before, " cells")
}
rm(obj_list); gc()

# ---- Merge all singlets into one object ----
All_Data <- Reduce(function(a, b) merge(a, b, project = "Retina_v2",
                                         merge.data = TRUE), singlet_list)
rm(singlet_list); gc()
write.csv(doublet_summary, file.path(TABLE_DIR, "01_doublet_summary.csv"),
          row.names = FALSE)
message("  Merged: ", ncol(All_Data), " cells x ", nrow(All_Data), " genes")

# ============================================================
# PHASE 3 — Quality control filtering
# ============================================================
message("\n=== PHASE 3: QC filtering ===")

DefaultAssay(All_Data) <- "RNA"
All_Data[["percent.mt"]] <- PercentageFeatureSet(All_Data, pattern = "^mt-")
All_Data[["percent.hb"]] <- PercentageFeatureSet(All_Data, pattern = "^Hb")
All_Data[["percent.rb"]] <- PercentageFeatureSet(All_Data, pattern = "^Rp[sl]")
gc()

# ---- QC figures BEFORE filtering ----
message("  Pre-QC figures ...")
meta_before <- All_Data@meta.data
p_vln_before <- meta_before %>%
  select(condition, nCount_RNA, nFeature_RNA, percent.mt, percent.hb) %>%
  pivot_longer(-condition, names_to = "metric", values_to = "value") %>%
  ggplot(aes(condition, value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "QC metrics BEFORE filtering", x = "", y = "")
save_fig(p_vln_before, "qc_violin_before", FIG_QC, width = 14, height = 5)

# ---- Apply filters ----
message("  Applying QC filters ...")
n_before <- ncol(All_Data)
meta <- All_Data@meta.data
keep <- (meta$nCount_RNA   > QC_MIN_COUNTS) &
        (meta$nCount_RNA   < QC_MAX_COUNTS) &
        (meta$nFeature_RNA > QC_MIN_FEATURES) &
        (meta$nFeature_RNA < QC_MAX_FEATURES) &
        (meta$percent.mt   < QC_MAX_PCT_MT) &
        (meta$percent.hb   < QC_MAX_PCT_HB)
All_Data <- subset(All_Data, cells = colnames(All_Data)[keep])
rm(meta, keep); gc()

cts <- GetAssayData(All_Data, layer = "counts")
keep_genes <- Matrix::rowSums(cts > 0) >= QC_GENE_MIN_CELLS
All_Data <- subset(All_Data, features = rownames(All_Data)[keep_genes])
rm(cts, keep_genes); gc()

rm_genes <- grepl("^mt-", rownames(All_Data)) | grepl("^Hb", rownames(All_Data))
All_Data <- subset(All_Data, features = rownames(All_Data)[!rm_genes])
rm(rm_genes); gc()

n_after <- ncol(All_Data)
message("  Cells: ", n_before, " -> ", n_after,
        " (removed ", n_before - n_after, ")")

# Recompute metrics + log-normalize for the RNA assay
All_Data[["percent.mt"]] <- PercentageFeatureSet(All_Data, pattern = "^mt-")
All_Data[["percent.hb"]] <- PercentageFeatureSet(All_Data, pattern = "^Hb")
All_Data <- NormalizeData(All_Data, verbose = FALSE)
gc()

filter_summary <- data.frame(
  metric        = c("nCount_RNA", "nFeature_RNA", "percent.mt",
                     "percent.hb", "gene_min_cells"),
  threshold     = c(paste0(QC_MIN_COUNTS, "-", QC_MAX_COUNTS),
                     paste0(QC_MIN_FEATURES, "-", QC_MAX_FEATURES),
                     paste0("<", QC_MAX_PCT_MT),
                     paste0("<", QC_MAX_PCT_HB),
                     paste0(">=", QC_GENE_MIN_CELLS)),
  cells_before  = n_before,
  cells_after   = n_after)
write.csv(filter_summary, file.path(TABLE_DIR, "01_qc_filter_summary.csv"),
          row.names = FALSE)

# ---- Post-QC violin ----
meta_after <- All_Data@meta.data
p_vln_after <- meta_after %>%
  select(condition, nCount_RNA, nFeature_RNA, percent.mt, percent.hb) %>%
  pivot_longer(-condition, names_to = "metric", values_to = "value") %>%
  ggplot(aes(condition, value, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.5) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "QC metrics AFTER filtering", x = "", y = "")
save_fig(p_vln_after, "qc_violin_after", FIG_QC, width = 14, height = 5)
rm(meta_before, meta_after); gc()

# ============================================================
# PHASE 4 — SCTransform + PCA + Harmony integration
# ============================================================
message("\n=== PHASE 4: SCTransform + PCA + Harmony ===")

# Split RNA assay by sample so SCTransform builds a per-sample model;
# this is the standard SCT-v2 / Harmony workflow.
All_Data[["RNA"]] <- split(All_Data[["RNA"]], f = All_Data$sample)

message("  SCTransform ...")
All_Data <- SCTransform(All_Data, vars.to.regress = "percent.mt",
                         verbose = FALSE)

message("  PCA (npcs = ", N_PCS_INTEGRATION, ") ...")
All_Data <- RunPCA(All_Data, npcs = N_PCS_INTEGRATION, verbose = FALSE)

message("  RunHarmony on PCA (vars_use = 'sample', theta = ",
        HARMONY_THETA, ", lambda = ", HARMONY_LAMBDA,
        ", max_iter = ", HARMONY_MAX_ITER, ") ...")
set.seed(SEED)
All_Data <- harmony::RunHarmony(
  object         = All_Data,
  group.by.vars  = "sample",
  reduction.use  = "pca",
  reduction.save = "harmony",
  dims.use       = 1:N_PCS_INTEGRATION,
  theta          = HARMONY_THETA,
  lambda         = HARMONY_LAMBDA,
  max_iter       = HARMONY_MAX_ITER,
  verbose        = FALSE)

message("  FindNeighbors + FindClusters on harmony reduction ...")
All_Data <- FindNeighbors(All_Data, reduction = "harmony",
                           dims = 1:N_DIMS_NEIGHBORS, verbose = FALSE)
All_Data <- FindClusters(All_Data, resolution = CLUSTER_RESOLUTIONS,
                          verbose = FALSE)

message("  RunUMAP on harmony reduction ...")
All_Data <- RunUMAP(All_Data, reduction = "harmony",
                     dims = 1:N_DIMS_NEIGHBORS, verbose = FALSE)

# Default identity = res 0.6 clusters
Idents(All_Data) <- All_Data[[paste0("SCT_snn_res.", DEFAULT_RESOLUTION)]][, 1]
All_Data$seurat_clusters <- Idents(All_Data)
gc()

# ---- Save integrated object ----
saveRDS(All_Data, file.path(RDS_DIR, "01_integrated.rds"))
message("\n  Saved: output/rds/01_integrated.rds")
message("  Cells: ", ncol(All_Data), " | Genes: ", nrow(All_Data),
        " | Clusters at res ", DEFAULT_RESOLUTION, ": ",
        length(unique(Idents(All_Data))))

# ---- Diagnostic UMAPs ----
message("\n  Diagnostic UMAP figures ...")

p_umap_clusters <- DimPlot(All_Data, reduction = "umap", label = TRUE,
                            label.size = 3, repel = TRUE) +
  pub_theme() +
  labs(title = paste0("Harmony UMAP — Seurat clusters (res ",
                       DEFAULT_RESOLUTION, ")"))
save_fig(p_umap_clusters, "umap_clusters", FIG_QC, width = 8, height = 7)

p_umap_sample <- DimPlot(All_Data, reduction = "umap", group.by = "sample") +
  pub_theme() + labs(title = "Harmony UMAP — by sample (batch check)")
save_fig(p_umap_sample, "umap_by_sample", FIG_QC, width = 8, height = 7)

p_umap_time <- DimPlot(All_Data, reduction = "umap", group.by = "time_label") +
  pub_theme() + labs(title = "Harmony UMAP — by timepoint")
save_fig(p_umap_time, "umap_by_timepoint", FIG_QC, width = 8, height = 7)

# ---- Marker dot plot (pre-annotation; informs stage 02 mapping) ----
DefaultAssay(All_Data) <- "SCT"
markers_present <- intersect(unlist(CELLTYPE_PANELS), rownames(All_Data))
p_dot <- DotPlot(All_Data, features = unique(markers_present),
                  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  pub_theme() +
  labs(title = paste0("Cell-type marker expression by cluster (res ",
                       DEFAULT_RESOLUTION, ")"))
save_fig(p_dot, "dotplot_markers_by_cluster", FIG_QC,
         width = 14, height = 8)

message("\n=== 01_qc_harmony_integration.R complete ===")
