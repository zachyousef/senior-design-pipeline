# ============================================================
# 04d_stability_figures.R   (pipeline_v2)
#
# Supplementary stability figure for the rod analysis, built from the
# permanent stability tables produced by stages 04b and 04c. No one-off
# data here — reads only committed CSV outputs.
#
#   (a) scDHA seed sensitivity: k=CHOSEN_K cluster sizes per training seed
#       (shows the smallest cluster is not reproducible across seeds).
#   (b) Pairwise ARI between seed partitions (heatmap).
#   (c) Within-embedding bootstrap stability: per-k mean ARI + jacc_min.
#
# Input:  output/tables/04c_seed_sizes.csv
#         output/tables/04c_seed_ari_matrix.csv
#         output/tables/04b_chosen_partition_stability.csv
# Output: manuscript/supplemental_figures/Supplemental_Fig_Stability.(pdf|png)
# ============================================================

source("config.R")
source("utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })

tdir <- TABLE_DIR
fig_out <- file.path(PIPELINE_DIR, "manuscript", "supplemental_figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

req <- file.path(tdir, c("04c_seed_sizes.csv", "04c_seed_ari_matrix.csv",
                         "04b_chosen_partition_stability.csv"))
if (!all(file.exists(req))) {
  stop("Missing stability tables; run stages 04b and 04c first:\n  ",
       paste(req[!file.exists(req)], collapse = "\n  "))
}

# ---- (a) seed sensitivity: k=3 sizes per seed ----
sizes <- read.csv(file.path(tdir, "04c_seed_sizes.csv"))
sz_long <- sizes %>%
  tidyr::separate_rows(sizes, sep = " / ") %>%
  group_by(seed) %>% mutate(rank = factor(row_number())) %>% ungroup() %>%
  mutate(n = as.numeric(sizes), seed = factor(seed))
pa <- ggplot(sz_long, aes(x = seed, y = n, fill = rank)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = "cluster\n(size rank)") +
  labs(x = "scDHA training seed", y = "cells",
       title = sprintf("a  Seed sensitivity (k = %d)", CHOSEN_K),
       subtitle = "smallest cluster is not reproducible across seeds") +
  pub_theme()

# ---- (b) pairwise ARI heatmap ----
ariM <- read.csv(file.path(tdir, "04c_seed_ari_matrix.csv"), check.names = FALSE)
seeds <- ariM$seed
mat <- as.matrix(ariM[, -1]); rownames(mat) <- seeds
ari_long <- as.data.frame(as.table(mat)); colnames(ari_long) <- c("seed_a", "seed_b", "ARI")
pb <- ggplot(ari_long, aes(seed_a, seed_b, fill = ARI)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI)), size = 3) +
  scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
                       midpoint = 0.5, limits = c(0, 1), name = "ARI") +
  labs(x = "seed", y = "seed", title = "b  Cross-seed partition agreement",
       subtitle = "pairwise adjusted Rand index") +
  pub_theme() + theme(legend.position = "right")

# ---- (c) bootstrap stability per k ----
boot <- read.csv(file.path(tdir, "04b_chosen_partition_stability.csv"))
boot_long <- boot %>%
  select(target_k, mean_ari, jacc_min) %>%
  tidyr::pivot_longer(c(mean_ari, jacc_min), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, mean_ari = "mean bootstrap ARI",
                         jacc_min = "min per-cluster Jaccard"))
pc <- ggplot(boot_long, aes(factor(target_k), value, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("mean bootstrap ARI" = "#7570B3",
                               "min per-cluster Jaccard" = "#1B9E77"), name = NULL) +
  ylim(0, 1) +
  labs(x = "target k", y = "stability",
       title = "c  Within-embedding bootstrap stability",
       subtitle = sprintf("n_iter = %d, 80%% resample", JACCARD_BOOTSTRAPS)) +
  pub_theme() + theme(legend.position = "top")

fig <- (pa | pb) / pc + plot_layout(heights = c(1, 0.9))
save_fig(fig, "Supplemental_Fig_Stability", fig_out, width = 11, height = 9)
message("Saved stability supplementary figure to ", fig_out)
