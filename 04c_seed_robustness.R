# ============================================================
# 04c_seed_robustness.R   (pipeline_v2)
#
# scDHA EMBEDDING-STABILITY analysis (methods/robustness section).
# Trains scDHA with multiple random seeds on the SAME cone-cleaned rod
# input, re-applies Harmony, clusters each embedding with Louvain at the
# chosen k (CHOSEN_K), and quantifies how reproducible the partition is
# across seeds (pairwise ARI/NMI + dominant/minor-cluster concordance).
#
# This is the definitive test of whether the rod sub-structure is a
# property of the data or of a single scDHA training run. Only the
# scDHA random seed varies between runs; input, Harmony params, and
# clustering are held fixed.
#
# Input:  output/rds/02_rods.rds
# Output: output/tables/04c_seed_sizes.csv          (k=CHOSEN_K sizes per seed)
#         output/tables/04c_seed_ari_matrix.csv     (pairwise ARI)
#         output/tables/04c_seed_robustness_summary.csv
# ============================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix)
  library(igraph); library(aricode)
})

SEED_LIST <- c(1, 2, 3, 7, 42)   # scDHA training seeds to compare

T0 <- Sys.time()
log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time() - T0, digits = 3), m))

# ---- Load rods + cone-contamination filter (identical to stage 04) ----
log_step("Loading rods + cone filter")
rods <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
DefaultAssay(rods) <- "RNA"
rods[["RNA"]] <- JoinLayers(rods[["RNA"]])
rods <- NormalizeData(rods, verbose = FALSE)

cone_g <- intersect(CONE_PANEL, rownames(rods))
rod_g  <- intersect(ROD_PANEL,  rownames(rods))
rods <- AddModuleScore(rods, features = list(cone_g), name = "cone_score", seed = SEED)
rods <- AddModuleScore(rods, features = list(rod_g),  name = "rod_score",  seed = SEED)
rods$cone_score <- rods$cone_score1; rods$rod_score <- rods$rod_score1
cone_q <- quantile(rods$cone_score, CONE_QUANTILE_CUT, na.rm = TRUE)
contam <- (rods$cone_score >= cone_q & rods$cone_score > rods$rod_score) |
          ((rods$cone_score - rods$rod_score) > CONE_DELTA_CUT)
clean <- subset(rods, cells = colnames(rods)[!contam])
log_step(sprintf("Clean rod cells: %d", ncol(clean)))
rm(rods); gc(verbose = FALSE)

har_meta <- data.frame(sample = clean$sample, row.names = colnames(clean))
counts_log2 <- log2(t(as.matrix(GetAssayData(clean, assay = "RNA", layer = "counts"))) + 1)

louvain_to_k <- function(latent, target_k) {
  obj <- CreateSeuratObject(matrix(0, 1, nrow(latent), dimnames = list("d", rownames(latent))))
  obj[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA")
  obj <- FindNeighbors(obj, reduction = "lat", dims = 1:ncol(latent), verbose = FALSE)
  bl <- NULL; bd <- Inf
  for (r in exp(seq(log(1e-4), log(10), length.out = 60))) {
    set.seed(SEED); obj <- FindClusters(obj, resolution = r, algorithm = 1, verbose = FALSE)
    lab <- as.integer(Idents(obj)); nk <- length(unique(lab))
    if (abs(nk - target_k) < bd) { bd <- abs(nk - target_k); bl <- lab }
    if (bd == 0) break
  }
  names(bl) <- rownames(latent); bl
}

suppressPackageStartupMessages({ library(torch); library(scDHA) })
if (!torch_is_installed()) install_torch()

parts <- list(); size_rows <- list()
for (s in SEED_LIST) {
  log_step(sprintf("=== scDHA seed %d ===", s))
  set.seed(s)
  sc <- scDHA(data = counts_log2, k = NULL, method = "scDHA", sparse = FALSE,
              n = SCDHA_N_GENES, ncores = N_CORES, gen_fil = TRUE,
              do.clus = TRUE, seed = s)
  lat <- as.matrix(sc$latent); rownames(lat) <- rownames(counts_log2)
  set.seed(s)
  har <- harmony::RunHarmony(lat, har_meta, vars_use = "sample",
                             theta = HARMONY_THETA, lambda = HARMONY_LAMBDA,
                             max_iter = HARMONY_MAX_ITER, verbose = FALSE)
  lab <- louvain_to_k(har, CHOSEN_K)
  parts[[as.character(s)]] <- lab
  sz <- sort(as.integer(table(lab)), decreasing = TRUE)
  size_rows[[length(size_rows)+1]] <- data.frame(
    seed = s, k = length(unique(lab)),
    sizes = paste(sz, collapse = " / "),
    smallest = min(sz), smallest_pct = round(100*min(sz)/length(lab), 2))
  log_step(sprintf("  seed %d k=3 sizes: %s", s, paste(sz, collapse = " / ")))
}

size_df <- do.call(rbind, size_rows)
write.csv(size_df, file.path(TABLE_DIR, "04c_seed_sizes.csv"), row.names = FALSE)

# Pairwise ARI / NMI on the common cells (identical for all seeds)
ids <- Reduce(intersect, lapply(parts, names))
sl <- names(parts)
ariM <- matrix(NA, length(sl), length(sl), dimnames = list(sl, sl))
nmiM <- ariM
for (i in seq_along(sl)) for (j in seq_along(sl)) {
  ariM[i,j] <- aricode::ARI(parts[[i]][ids], parts[[j]][ids])
  nmiM[i,j] <- aricode::NMI(parts[[i]][ids], parts[[j]][ids])
}
write.csv(data.frame(seed = sl, round(ariM, 4)),
          file.path(TABLE_DIR, "04c_seed_ari_matrix.csv"), row.names = FALSE)

offdiag <- ariM[upper.tri(ariM)]
summary_df <- data.frame(
  metric = c("n_seeds","mean_pairwise_ARI","min_pairwise_ARI","max_pairwise_ARI","mean_pairwise_NMI"),
  value  = c(length(sl), round(mean(offdiag),4), round(min(offdiag),4),
             round(max(offdiag),4), round(mean(nmiM[upper.tri(nmiM)]),4)))
write.csv(summary_df, file.path(TABLE_DIR, "04c_seed_robustness_summary.csv"), row.names = FALSE)

log_step("=== SEED ROBUSTNESS SUMMARY ===")
print(size_df); cat("\nPairwise ARI:\n"); print(round(ariM,3)); cat("\n"); print(summary_df)
saveRDS(list(parts=parts, sizes=size_df, ari=ariM, nmi=nmiM, summary=summary_df),
        file.path(RDS_DIR, "04c_seed_robustness.rds"))
log_step("DONE")
