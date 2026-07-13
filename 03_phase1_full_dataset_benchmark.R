# ============================================================
# 03_phase1_full_dataset_benchmark.R   (pipeline_v2)
#
# PHASE 1 — Positive control on the full retinal dataset.
#
# Run the same 7-method × multi-k benchmarking pipeline on the
# WHOLE retinal dataset (~56k cells, 8 known cell types) instead
# of the rod subset. The auto-annotated cell types from stage 02
# serve as ground truth, so we can additionally compute Adjusted
# Rand Index (ARI) and Normalized Mutual Information (NMI) for
# every (method, k) — which Phase 2 cannot, since rods have no
# ground-truth subtype annotation.
#
# Expected result if the pipeline is sound:
#   - At k around 8 (the number of known cell types), some (method, k)
#     pair should achieve ARI close to 1 vs the cell-type annotation.
#   - The four-metric panel (ROGUE / balance / Jaccard / marker_spec)
#     should jointly favor that same configuration without ever seeing
#     the ground truth.
#
# Methods:  all 7 (METHODS in config.R)
# k-range:  POC_K_RANGE (config.R, default 2..12)
# Latents:  Harmony PCA from stage 01 (already trained); plus a
#           freshly-trained scDHA on the full dataset (non-deterministic
#           — trained once and frozen via a sidecar lock file).
#
# Output:
#   output_full_dataset/tables/poc_scores_matrix.csv          (sil + ROGUE + ARI + NMI)
#   output_full_dataset/tables/poc_scores_matrix_4metric.csv  (4-metric panel + ARI + NMI)
#   output_full_dataset/tables/poc_all_cluster_assignments.csv
#   output_full_dataset/tables/poc_harmony_pca_latent.csv
#   output_full_dataset/tables/poc_harmony_scdha_latent.csv
#   output_full_dataset/tables/poc_scdha_latent.csv
#   output_full_dataset/rds/poc_benchmarking.rds
#
# Wall-time on a workstation: ~12 hours (77 partitions × 5 metrics
# + ARI/NMI per partition).
# ============================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix)
  library(cluster); library(mclust); library(igraph)
  library(uwot); library(dplyr); library(tidyr)
  library(ROGUE); library(aricode); library(RhpcBLASctl)
})

set.seed(SEED)
blas_set_num_threads(max(1L, floor(N_CORES / 2)))
omp_set_num_threads(max(1L, floor(N_CORES / 2)))

OUT <- "output_full_dataset"
dir.create(file.path(OUT, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "rds"),     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

T0 <- Sys.time()
log_step <- function(m) message(sprintf("[%s] %s",
                                          format(Sys.time() - T0, digits = 3), m))

POC_K_RANGE <- 2:12   # the proof-of-concept sweep range

# ============================================================
# PHASE 1 — Load full annotated object + ground truth
# ============================================================
log_step("=== PHASE 1: Loading annotated object ===")

obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
log_step(sprintf("  Cells: %d  |  Genes: %d", ncol(obj), nrow(obj)))

# Ground-truth cell-type labels (from stage 02 module-score auto-annotation)
gt_labels <- as.character(obj$cell_type)
names(gt_labels) <- colnames(obj)
n_celltypes <- length(unique(gt_labels))
log_step(sprintf("  Ground-truth cell types: %d", n_celltypes))
print(table(gt_labels))

# ============================================================
# PHASE 2 — Pull existing Harmony PCA, train fresh scDHA, Harmony scDHA
# ============================================================
log_step("=== PHASE 2: Building latents ===")

# Harmony PCA — already computed in stage 01
har_pca <- Embeddings(obj, "harmony")
log_step(sprintf("  Harmony PCA: %d cells x %d dims",
                 nrow(har_pca), ncol(har_pca)))

# Cap to N_PCS_CLUSTERING dims for parity with stage 03's rod benchmark
har_pca <- har_pca[, 1:min(N_PCS_CLUSTERING, ncol(har_pca)), drop = FALSE]
write.csv(data.frame(cell_id = rownames(har_pca), har_pca),
          file.path(OUT, "tables", "poc_harmony_pca_latent.csv"),
          row.names = FALSE)

# scDHA — train fresh on full data, then freeze
scdha_path <- file.path(OUT, "tables", "poc_scdha_latent.csv")
scdha_lock <- file.path(OUT, "tables", "_scdha_latent.lock")
if (file.exists(scdha_path) && file.exists(scdha_lock)) {
  log_step("  scDHA: loading frozen ...")
  df <- read.csv(scdha_path, check.names = FALSE,
                  stringsAsFactors = FALSE)
  scdha_lat <- as.matrix(df[, -1, drop = FALSE])
  rownames(scdha_lat) <- as.character(df[[1]])
  stopifnot(setequal(rownames(scdha_lat), colnames(obj)))
  scdha_lat <- scdha_lat[colnames(obj), , drop = FALSE]
} else {
  log_step("  scDHA: training fresh ...")
  suppressPackageStartupMessages({
    library(torch); library(scDHA)
  })
  data_scdha <- t(as.matrix(GetAssayData(obj, assay = "RNA",
                                          layer = "counts")))
  data_scdha <- log2(data_scdha + 1)
  set.seed(SEED)
  sc_res <- scDHA(data    = data_scdha,
                   k       = NULL,
                   method  = "scDHA",
                   sparse  = FALSE,
                   n       = SCDHA_N_GENES,
                   ncores  = N_CORES,
                   gen_fil = TRUE,
                   do.clus = TRUE,
                   seed    = SEED)
  scdha_lat <- as.matrix(sc_res$latent)
  rownames(scdha_lat) <- rownames(data_scdha)
  rm(data_scdha); gc()

  write.csv(data.frame(cell_id = rownames(scdha_lat), scdha_lat),
            scdha_path, row.names = FALSE)
  writeLines(c(paste0("trained_at=", Sys.time()),
                paste0("seed=", SEED),
                paste0("n_cells=", nrow(scdha_lat)),
                paste0("dims=", ncol(scdha_lat))),
              scdha_lock)
}
log_step(sprintf("  scDHA latent: %d cells x %d dims",
                 nrow(scdha_lat), ncol(scdha_lat)))

# Apply Harmony to scDHA latent
har_meta <- data.frame(sample = obj$sample, row.names = colnames(obj))
set.seed(SEED)
har_scdha <- harmony::RunHarmony(scdha_lat, har_meta, vars_use = "sample",
                                  theta = HARMONY_THETA,
                                  lambda = HARMONY_LAMBDA,
                                  max_iter = HARMONY_MAX_ITER,
                                  verbose = FALSE)
write.csv(data.frame(cell_id = rownames(har_scdha), har_scdha),
          file.path(OUT, "tables", "poc_harmony_scdha_latent.csv"),
          row.names = FALSE)
log_step("  Harmony applied to scDHA latent")

cell_ids   <- colnames(obj)
counts_raw <- GetAssayData(obj, assay = "RNA", layer = "counts")
sample_vec <- as.character(obj$sample)

# ============================================================
# PHASE 3 — ROGUE prefilter
# ============================================================
log_step("=== PHASE 3: ROGUE prefilter ===")
counts_dense <- as.matrix(counts_raw)
rogue_expr   <- ROGUE::matr.filter(counts_dense, min.cells = 10,
                                    min.genes = 10)
ROGUE::SE_fun(rogue_expr)
log_step(sprintf("  ROGUE-filtered: %d genes x %d cells",
                 nrow(rogue_expr), ncol(rogue_expr)))
rm(counts_dense); gc()

# ROGUE baseline on the unsplit population
rogue_baseline <- tryCatch({
  one_lab <- rep("1", ncol(rogue_expr))
  r <- ROGUE::rogue(rogue_expr, labels = one_lab,
                    samples = rep(1, length(one_lab)),
                    platform = "UMI")
  mean(unlist(r), na.rm = TRUE)
}, error = function(e) NA_real_)
log_step(sprintf("  ROGUE baseline (unsplit): %.4f", rogue_baseline))

# ============================================================
# PHASE 4 — Build SNN graphs (for Leiden methods)
# ============================================================
log_step("=== PHASE 4: SNN graphs ===")
g_har_pca   <- build_seurat_snn_igraph(har_pca,   cell_ids,
                                        counts_for_dummy = counts_raw,
                                        reduction_name = "harpca",
                                        key_prefix = "HARPC_",
                                        dims_use = 1:N_DIMS_NEIGHBORS)
log_step(sprintf("  PCA-Harmony SNN: %d edges", igraph::ecount(g_har_pca)))
g_har_scdha <- build_seurat_snn_igraph(har_scdha, cell_ids,
                                        counts_for_dummy = counts_raw,
                                        reduction_name = "harscdha",
                                        key_prefix = "HARSC_",
                                        dims_use = 1:ncol(har_scdha))
log_step(sprintf("  scDHA-Harmony SNN: %d edges",
                 igraph::ecount(g_har_scdha)))

# ============================================================
# Helpers
# ============================================================
gini <- function(x) {
  x <- sort(x); n <- length(x)
  if (n <= 1L || sum(x) == 0) return(0)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}
balance_of <- function(l) 1 - gini(as.numeric(table(l)))

marker_specificity <- function(seurat_obj, labels) {
  Idents(seurat_obj) <- factor(labels)
  m <- tryCatch(
    FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25,
                   logfc.threshold = 0.5, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0) return(0)
  m <- m[m$p_val_adj < 0.05, , drop = FALSE]
  if (nrow(m) == 0) return(0)
  m$score <- (m$pct.1 - m$pct.2) * abs(m$avg_log2FC)
  per_cl <- aggregate(score ~ cluster, data = m, FUN = mean)
  mean(per_cl$score)
}

recluster_subsample <- function(method, latent_sub, target_k) {
  if (target_k == 1) return(rep(1L, nrow(latent_sub)))
  if (method %in% c("PCA_KMeans", "scDHA_KMeans")) {
    set.seed(SEED)
    return(kmeans(latent_sub, centers = target_k,
                  nstart = 10, iter.max = 100)$cluster)
  }
  if (method == "PCA_GMM") {
    g <- tryCatch(
      Mclust(latent_sub[, 1:min(20, ncol(latent_sub))],
              G = target_k, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(g)) return(rep(1L, nrow(latent_sub)))
    return(as.integer(g$classification))
  }
  if (method %in% c("PCA_LL", "scDHA_LL")) {
    obj <- CreateSeuratObject(matrix(0, 1, nrow(latent_sub),
                              dimnames = list("d", rownames(latent_sub))))
    obj[["lat"]] <- CreateDimReducObject(embeddings = latent_sub,
                                          key = "L_", assay = "RNA")
    obj <- FindNeighbors(obj, reduction = "lat",
                          dims = 1:ncol(latent_sub), verbose = FALSE)
    best_lab <- NULL; best_diff <- Inf
    for (r in exp(seq(log(1e-3), log(2), length.out = 30))) {
      set.seed(SEED)
      obj <- FindClusters(obj, resolution = r, algorithm = 1,
                           verbose = FALSE)
      lab <- as.integer(Idents(obj)); k_now <- length(unique(lab))
      if (abs(k_now - target_k) < best_diff) {
        best_diff <- abs(k_now - target_k); best_lab <- lab
      }
      if (best_diff == 0) break
    }
    return(best_lab)
  }
  if (method %in% c("PCA_Leiden", "scDHA_Leiden")) {
    obj <- CreateSeuratObject(matrix(0, 1, nrow(latent_sub),
                              dimnames = list("d", rownames(latent_sub))))
    obj[["lat"]] <- CreateDimReducObject(embeddings = latent_sub,
                                          key = "L_", assay = "RNA")
    obj <- FindNeighbors(obj, reduction = "lat",
                          dims = 1:ncol(latent_sub), verbose = FALSE)
    snn_mat <- methods::as(obj@graphs[["RNA_snn"]], "TsparseMatrix")
    keep    <- snn_mat@i < snn_mat@j
    edges   <- cbind(snn_mat@i[keep] + 1L, snn_mat@j[keep] + 1L)
    g_sub   <- igraph::make_graph(t(edges), n = ncol(snn_mat),
                                    directed = FALSE)
    igraph::E(g_sub)$weight <- snn_mat@x[keep]
    igraph::V(g_sub)$name   <- rownames(snn_mat)
    res <- find_resolution_for_k(g_sub, target_k = target_k,
                                  algorithm = "leiden",
                                  res_min = 1e-6, res_max = 5,
                                  max_iter = 20, seed = SEED)
    return(as.integer(res$labels))
  }
  rep(1L, nrow(latent_sub))
}

per_cluster_jaccard <- function(method, latent, ref_labels, target_k,
                                 n_iter = JACCARD_BOOTSTRAPS, frac = 0.8) {
  n <- nrow(latent)
  ref <- as.character(ref_labels)
  per_iter_max <- vector("list", length(unique(ref)))
  names(per_iter_max) <- unique(ref)
  for (cl in names(per_iter_max)) per_iter_max[[cl]] <- numeric(0)
  set.seed(SEED)
  for (it in seq_len(n_iter)) {
    samp_idx <- sample(seq_len(n), floor(frac * n))
    sub_lat  <- latent[samp_idx, , drop = FALSE]
    new_lab  <- recluster_subsample(method, sub_lat, target_k)
    new_lab  <- as.character(new_lab)
    ref_sub  <- ref[samp_idx]
    for (cl in unique(ref_sub)) {
      ref_set <- which(ref_sub == cl)
      jac <- vapply(unique(new_lab), function(nc) {
        new_set <- which(new_lab == nc)
        length(intersect(ref_set, new_set)) /
          length(union(ref_set, new_set))
      }, numeric(1))
      per_iter_max[[cl]] <- c(per_iter_max[[cl]], max(jac))
    }
  }
  per_cluster_mean <- sapply(per_iter_max, mean)
  if (length(per_cluster_mean) == 0) return(NA_real_)
  min(per_cluster_mean)
}

# ============================================================
# PHASE 5 — Sweep
# ============================================================
log_step("=== PHASE 5: Clustering sweep + metrics + ARI/NMI ===")

method_to_latent <- function(m) {
  if (grepl("scDHA", m)) har_scdha else har_pca
}
method_to_graph <- function(m) {
  if (grepl("scDHA", m)) g_har_scdha else g_har_pca
}

silrogue_rows  <- list()
fourmet_rows   <- list()
all_assignments <- data.frame(cell_id = cell_ids)

checkpoint <- function() {
  if (length(silrogue_rows) > 0) {
    write.csv(do.call(rbind, silrogue_rows),
              file.path(OUT, "tables", "poc_scores_matrix.csv"),
              row.names = FALSE)
  }
  if (length(fourmet_rows) > 0) {
    write.csv(do.call(rbind, fourmet_rows),
              file.path(OUT, "tables", "poc_scores_matrix_4metric.csv"),
              row.names = FALSE)
  }
  write.csv(all_assignments,
            file.path(OUT, "tables", "poc_all_cluster_assignments.csv"),
            row.names = FALSE)
}

for (mname in METHODS) {
  lat <- method_to_latent(mname)
  for (k in POC_K_RANGE) {
    log_step(sprintf("--- %s @ k=%d ---", mname, k))
    labels <- NULL

    tryCatch({
      if (mname == "PCA_KMeans") {
        set.seed(SEED)
        labels <- kmeans(lat, centers = k, nstart = 25,
                         iter.max = 100)$cluster
      } else if (mname == "PCA_GMM") {
        set.seed(SEED)
        gmm <- tryCatch(
          Mclust(lat[, 1:min(20, ncol(lat))], G = k, verbose = FALSE),
          error = function(e) NULL)
        labels <- if (is.null(gmm)) rep(1L, nrow(lat)) else
                  as.integer(gmm$classification)
      } else if (mname == "scDHA_KMeans") {
        set.seed(SEED)
        labels <- kmeans(lat, centers = k, nstart = 25,
                         iter.max = 100)$cluster
      } else if (mname %in% c("PCA_LL", "scDHA_LL")) {
        obj_tmp <- CreateSeuratObject(counts = counts_raw[, cell_ids])
        obj_tmp[["lat"]] <- CreateDimReducObject(embeddings = lat,
                                                   key = "L_",
                                                   assay = "RNA")
        obj_tmp <- FindNeighbors(obj_tmp, reduction = "lat",
                                  dims = 1:ncol(lat), verbose = FALSE)
        best_lab <- NULL; best_diff <- Inf; best_res <- NA_real_
        for (r in exp(seq(log(1e-4), log(10), length.out = 60))) {
          set.seed(SEED)
          obj_tmp <- FindClusters(obj_tmp, resolution = r, algorithm = 1,
                                   verbose = FALSE)
          lab <- as.integer(Idents(obj_tmp))
          nk  <- length(unique(lab))
          if (abs(nk - k) < best_diff) {
            best_diff <- abs(nk - k); best_lab <- lab; best_res <- r
          }
          if (best_diff == 0) break
        }
        labels <- best_lab
        rm(obj_tmp); gc(verbose = FALSE)
      } else if (mname %in% c("PCA_Leiden", "scDHA_Leiden")) {
        g <- method_to_graph(mname)
        res <- find_resolution_for_k(g, target_k = k, algorithm = "leiden",
                                      res_min = 1e-6, res_max = 5,
                                      max_iter = 30, seed = SEED)
        labels <- as.integer(res$labels[cell_ids])
      }
      names(labels) <- cell_ids
      actual_k <- length(unique(labels))
      col <- paste0(mname, "_k", k)
      all_assignments[[col]] <- labels[cell_ids]

      # Metrics
      sil <- compute_silhouette(lat, labels,
                                 max_cells = MAX_SIL_CELLS, seed = SEED)
      rg <- tryCatch({
        rogue_lab <- labels[colnames(rogue_expr)]
        r <- ROGUE::rogue(rogue_expr,
                          labels   = as.character(rogue_lab),
                          samples  = rep(1, length(rogue_lab)),
                          platform = "UMI")
        mean(unlist(r), na.rm = TRUE)
      }, error = function(e) NA_real_)
      bal <- balance_of(labels)
      ms  <- marker_specificity(obj, labels)
      jac <- tryCatch(
        per_cluster_jaccard(mname, lat, labels, target_k = k,
                             n_iter = BENCH_JACCARD_BOOTSTRAPS),
        error = function(e) NA_real_)

      # GROUND-TRUTH agreement metrics
      ari <- tryCatch(
        aricode::ARI(as.character(labels), as.character(gt_labels)),
        error = function(e) NA_real_)
      nmi <- tryCatch(
        aricode::NMI(as.character(labels), as.character(gt_labels)),
        error = function(e) NA_real_)

      log_step(sprintf("  sil=%.3f rogue=%.3f bal=%.3f jacc=%.3f ms=%.3f ARI=%.3f NMI=%.3f",
                        sil, rg, bal, jac, ms, ari, nmi))

      silrogue_rows[[length(silrogue_rows) + 1L]] <- data.frame(
        method = mname, target_k = k, actual_k = actual_k,
        silhouette = round(sil, 4), rogue = round(rg, 4),
        ari = round(ari, 4), nmi = round(nmi, 4),
        stringsAsFactors = FALSE)
      fourmet_rows[[length(fourmet_rows) + 1L]] <- data.frame(
        method = mname, target_k = k, actual_k = actual_k,
        balance = round(bal, 4), jacc_min = round(jac, 4),
        marker_spec = round(ms, 4), rogue = round(rg, 4),
        ari = round(ari, 4), nmi = round(nmi, 4),
        sizes_str = paste(table(labels), collapse = "/"),
        stringsAsFactors = FALSE)
      checkpoint()

    }, error = function(e) {
      log_step(sprintf("  ERROR for %s k=%d: %s", mname, k, e$message))
      silrogue_rows[[length(silrogue_rows) + 1L]] <<- data.frame(
        method = mname, target_k = k, actual_k = NA,
        silhouette = NA, rogue = NA, ari = NA, nmi = NA,
        stringsAsFactors = FALSE)
      fourmet_rows[[length(fourmet_rows) + 1L]] <<- data.frame(
        method = mname, target_k = k, actual_k = NA,
        balance = NA, jacc_min = NA, marker_spec = NA, rogue = NA,
        ari = NA, nmi = NA, sizes_str = NA_character_,
        stringsAsFactors = FALSE)
      checkpoint()
    })
  }
}

# Final saves
write.csv(do.call(rbind, silrogue_rows),
          file.path(OUT, "tables", "poc_scores_matrix.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, fourmet_rows),
          file.path(OUT, "tables", "poc_scores_matrix_4metric.csv"),
          row.names = FALSE)
write.csv(all_assignments,
          file.path(OUT, "tables", "poc_all_cluster_assignments.csv"),
          row.names = FALSE)

bundle <- list(
  scores_silrogue     = do.call(rbind, silrogue_rows),
  scores_4metric      = do.call(rbind, fourmet_rows),
  cluster_assignments = all_assignments,
  har_pca             = har_pca,
  har_scdha           = har_scdha,
  cell_ids            = cell_ids,
  ground_truth        = gt_labels,
  rogue_baseline      = rogue_baseline,
  n_celltypes         = n_celltypes
)
saveRDS(bundle, file.path(OUT, "rds", "poc_benchmarking.rds"))

log_step(sprintf("=== DONE — total %s ===",
                  format(Sys.time() - T0, digits = 3)))
