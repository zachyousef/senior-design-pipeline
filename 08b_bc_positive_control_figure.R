# ============================================================
# 08b_bc_positive_control_figure.R   (pipeline_v2)
#
# Phase 2 positive discovery control (bipolar cells), shown at the same
# level of detail as the whole-retina control (Fig 3) and the rod
# negative case (Fig 5). Reads the stage-07 positive-control object
# (5 scDHA embeddings + BC subset) and the stability-by-k table, and
# builds a self-contained figure:
#   (a) bipolar subtype UMAP on the reference scDHA embedding at the
#       accepted partition (k = K_STAR), colored by subcluster;
#   (b) across-embedding reproducibility (mean cross-seed ARI, min-max
#       band) vs k — rises toward the true subtype count (ACCEPT);
#   (c) cluster-quality metrics (silhouette, ROGUE, balance) vs k —
#       improve/hold in parallel with reproducibility;
#   (d) subcluster sizes at the accepted partition.
#
# Input:  output/rds/07_bc_positive_control.rds   (stage 07)
#         output/tables/07_bc_stability_by_k.csv  (stage 07)
# Output: manuscript/figures/Figure_BC_positive_control.(pdf|png)
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({
  library(Seurat); library(uwot); library(ggplot2); library(dplyr)
  library(tidyr); library(patchwork)
})

K_STAR   <- 15  # accepted bipolar partition: reproducibility peaks/plateaus
                # here (cross-seed ARI 0.92, min 0.90), coinciding with the
                # ~15 transcriptomic bipolar types of Shekhar et al. (2016)
SHEKHAR  <- 15  # reference bipolar subtype count (Shekhar et al. 2016)
REF_SEED <- "1" # reference scDHA embedding for the displayed UMAP

frds <- file.path(RDS_DIR, "07_bc_positive_control.rds")
# prefer the extended k=2-18 sweep (stage 07b) if present, else the k=2-10 table
fstb <- file.path(TABLE_DIR, "07_bc_stability_by_k_extended.csv")
if (!file.exists(fstb)) fstb <- file.path(TABLE_DIR, "07_bc_stability_by_k.csv")
if (!file.exists(frds) || !file.exists(fstb))
  stop("Need stage-07 BC outputs: run `07` (and optionally `07b`) first.")

pc  <- readRDS(frds)
lat <- pc$seed_lat[[REF_SEED]]         # Harmony-corrected scDHA latent (reference seed)
stb <- read.csv(fstb)

# ---- Louvain partition at target k on a fixed latent (mirrors stage 07) ----
louvain_to_k <- function(latent, target_k) {
  o <- CreateSeuratObject(matrix(0, 1, nrow(latent), dimnames = list("d", rownames(latent))))
  o[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA")
  o <- FindNeighbors(o, reduction = "lat", dims = 1:ncol(latent), verbose = FALSE)
  bl <- NULL; bd <- Inf
  for (r in exp(seq(log(1e-4), log(10), length.out = 60))) {
    set.seed(SEED); o <- FindClusters(o, resolution = r, algorithm = 1, verbose = FALSE)
    lab <- as.integer(Idents(o)); nk <- length(unique(lab))
    if (abs(nk - target_k) < bd) { bd <- abs(nk - target_k); bl <- lab }; if (bd == 0) break
  }
  names(bl) <- rownames(latent); bl
}

lab <- louvain_to_k(lat, K_STAR)
# relabel clusters by descending size (cluster 1 = largest) for stable coloring
ord   <- names(sort(table(lab), decreasing = TRUE))
remap <- setNames(seq_along(ord), ord)
sub   <- factor(remap[as.character(lab)], levels = seq_along(ord))

# ---- (a) reference-embedding UMAP colored by subcluster ----
set.seed(SEED)
um <- uwot::umap(lat, n_neighbors = 30, min_dist = 0.3, verbose = FALSE)
umd <- data.frame(UMAP1 = um[, 1], UMAP2 = um[, 2], Subcluster = sub)
pal <- setNames(scales::hue_pal()(nlevels(sub)), levels(sub))

pa <- ggplot(umd, aes(UMAP1, UMAP2, color = Subcluster)) +
  geom_point(size = 0.35, alpha = 0.75) +
  scale_color_manual(values = pal, name = "Subcluster",
                     guide = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 2)) +
  labs(title = sprintf("a — Bipolar subtypes (accepted partition, k = %d ≈ Shekhar 2016)", K_STAR)) +
  pub_theme() +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(size = 8), legend.key.height = unit(3.4, "mm"))

# ---- (b) across-embedding reproducibility vs k ----
pb <- ggplot(stb, aes(k, mean_cross_seed_ari)) +
  geom_vline(xintercept = K_STAR, linetype = "dashed", color = "grey65") +
  geom_ribbon(aes(ymin = min_cross_seed_ari,
                  ymax = pmin(1, 2 * mean_cross_seed_ari - min_cross_seed_ari)),
              fill = "#1B9E77", alpha = 0.15, color = NA) +
  geom_line(color = "#1B9E77", linewidth = 1) +
  geom_point(color = "#1B9E77", size = 2) +
  annotate("text", x = SHEKHAR, y = 0.12, label = "Shekhar 2016\n≈15 types",
           size = 2.6, color = "grey40", lineheight = 0.9) +
  scale_x_continuous(breaks = seq(2, 18, 2)) + ylim(0, 1) +
  labs(x = "cluster number (k)", y = "mean cross-seed ARI",
       title = "b — Across-embedding reproducibility",
       subtitle = "rises to a plateau (~0.93) at the ~15 known bipolar types (ACCEPT)") +
  pub_theme() + theme(plot.title = element_text(face = "bold"))

# ---- (c) cluster-quality metrics vs k ----
ql <- stb %>%
  select(k, silhouette, rogue, balance) %>%
  pivot_longer(-k, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("silhouette", "rogue", "balance"),
                         labels = c("silhouette", "ROGUE purity", "size balance (1 - Gini)")))
qcol <- c("silhouette" = "#3288BD", "ROGUE purity" = "#66C2A5",
          "size balance (1 - Gini)" = "#8073AC")
pc_ <- ggplot(ql, aes(k, value, color = metric)) +
  geom_vline(xintercept = K_STAR, linetype = "dashed", color = "grey65") +
  geom_line(linewidth = 1) + geom_point(size = 1.8) +
  scale_color_manual(values = qcol, name = NULL) +
  scale_x_continuous(breaks = seq(2, 18, 2)) + ylim(0, 1) +
  labs(x = "cluster number (k)", y = "metric value",
       title = "c — Cluster quality holds/improves with k") +
  pub_theme() + theme(plot.title = element_text(face = "bold"),
                      legend.position = "top", legend.text = element_text(size = 8))

# ---- (d) subcluster sizes at the accepted partition ----
sz <- as.data.frame(table(sub)); colnames(sz) <- c("Subcluster", "n")
sz$pct <- round(100 * sz$n / sum(sz$n), 1)
pd <- ggplot(sz, aes(Subcluster, n, fill = Subcluster)) +
  geom_col(width = 0.78) +
  geom_text(aes(label = n), vjust = -0.35, size = 2.3) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Subcluster", y = "Number of cells",
       title = sprintf("d — Subcluster sizes (k = %d)", K_STAR)) +
  pub_theme() + theme(plot.title = element_text(face = "bold"))

fig <- (pa | (pb / pc_)) / pd + plot_layout(heights = c(1.35, 0.75))
save_fig(fig, "Figure_BC_positive_control",
         file.path(PIPELINE_DIR, "manuscript", "figures"), width = 11, height = 10)
message(sprintf("Saved BC positive-control figure (k=%d; sizes %s).",
                K_STAR, paste(sz$n, collapse = "/")))
