# ============================================================
# 04b_chosen_partition_stability.R   (pipeline_v2)
#
# BOOTSTRAP STABILITY of the chosen clustering method across k.
# For CHOSEN_METHOD (scDHA + Louvain) at each k in STABILITY_K_RANGE,
# resamples JACCARD_BOOTSTRAPS times (80% without replacement),
# re-clusters, and reports the FULL per-cluster Jaccard (not just the
# minimum) plus mean/sd of the bootstrap-vs-reference ARI. This is the
# within-embedding reproducibility test that feeds the manuscript's
# stability figure; the across-embedding (seed) test is 04c.
#
# Input:  output/tables/03_harmony_scdha_latent.csv  (from stage 04)
# Output: output/tables/04b_chosen_partition_stability.csv
# ============================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat); library(igraph); library(aricode)
})

T0 <- Sys.time()
log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time() - T0, digits = 3), m))

latent_path <- file.path(TABLE_DIR, "03_harmony_scdha_latent.csv")
stopifnot(file.exists(latent_path))
df  <- read.csv(latent_path, check.names = FALSE)
lat <- as.matrix(df[, -1]); rownames(lat) <- as.character(df[[1]])
log_step(sprintf("Loaded harmony-scDHA latent: %d x %d", nrow(lat), ncol(lat)))

# Louvain to a target k on a latent (mirrors stage 04 main-sweep logic)
louvain_to_k <- function(latent, target_k, grid_n = 60, lo = 1e-4, hi = 10) {
  obj <- CreateSeuratObject(matrix(0, 1, nrow(latent), dimnames = list("d", rownames(latent))))
  obj[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA")
  obj <- FindNeighbors(obj, reduction = "lat", dims = 1:ncol(latent), verbose = FALSE)
  bl <- NULL; bd <- Inf
  for (r in exp(seq(log(lo), log(hi), length.out = grid_n))) {
    set.seed(SEED); obj <- FindClusters(obj, resolution = r, algorithm = 1, verbose = FALSE)
    lab <- as.integer(Idents(obj)); nk <- length(unique(lab))
    if (abs(nk - target_k) < bd) { bd <- abs(nk - target_k); bl <- lab }
    if (bd == 0) break
  }
  names(bl) <- rownames(latent); bl
}

stability_at_k <- function(latent, target_k, n_iter = JACCARD_BOOTSTRAPS, frac = 0.8) {
  ref  <- as.character(louvain_to_k(latent, target_k))
  names(ref) <- rownames(latent)
  refsz <- sort(table(ref), decreasing = TRUE)
  n <- nrow(latent)
  per_cl <- setNames(vector("list", length(unique(ref))), unique(ref))
  aris <- numeric(n_iter)
  set.seed(SEED)
  for (it in seq_len(n_iter)) {
    idx <- sample(seq_len(n), floor(frac * n))
    sub <- latent[idx, , drop = FALSE]
    nl  <- as.character(louvain_to_k(sub, target_k, grid_n = 30, lo = 1e-3, hi = 2))
    rs  <- ref[idx]
    aris[it] <- aricode::ARI(rs, nl)
    for (cl in unique(rs)) {
      rset <- which(rs == cl)
      jmax <- max(vapply(unique(nl), function(nc) {
        ns <- which(nl == nc)
        length(intersect(rset, ns)) / length(union(rset, ns))
      }, numeric(1)))
      per_cl[[cl]] <- c(per_cl[[cl]], jmax)
    }
  }
  pcm <- sapply(per_cl, mean)
  # order per-cluster jaccard by descending cluster size
  cl_by_size <- names(sort(table(ref), decreasing = TRUE))
  data.frame(
    method = CHOSEN_METHOD, target_k = target_k,
    sizes = paste(as.integer(refsz), collapse = " / "),
    jacc_per_cluster = paste(round(pcm[cl_by_size], 4), collapse = " / "),
    jacc_min = round(min(pcm), 4),
    mean_ari = round(mean(aris), 4), sd_ari = round(sd(aris), 4),
    n_iter = n_iter, stringsAsFactors = FALSE)
}

rows <- list()
for (k in STABILITY_K_RANGE) {
  log_step(sprintf("=== bootstrap stability @ k=%d (n_iter=%d) ===", k, JACCARD_BOOTSTRAPS))
  rows[[length(rows) + 1]] <- stability_at_k(lat, k)
}
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "04b_chosen_partition_stability.csv"), row.names = FALSE)
log_step("=== STABILITY TABLE ==="); print(out)
log_step("DONE")
