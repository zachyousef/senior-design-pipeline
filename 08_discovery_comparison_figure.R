# ============================================================
# 08_discovery_comparison_figure.R   (pipeline_v2)
#
# Figure 5: across-embedding stability discriminates real from
# artifactual sub-structure. Reads the discovery-stability-by-k tables
# produced by stage 07 for bipolar cells (positive control; real
# subtypes) and rods (negative case), and plots (a) mean cross-seed ARI
# vs k and (b) silhouette vs k for both cell types on shared axes.
#
# Input:  output/tables/07_bc_stability_by_k.csv
#         output/tables/07_rods_stability_by_k.csv
# Output: manuscript/figures/Figure_5_discovery_comparison.(pdf|png)
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })

tdir <- TABLE_DIR
# prefer the extended BC sweep (k=2-18, stage 07b) if present
fbc <- file.path(tdir, "07_bc_stability_by_k_extended.csv")
if (!file.exists(fbc)) fbc <- file.path(tdir, "07_bc_stability_by_k.csv")
frd <- file.path(tdir, "07_rods_stability_by_k.csv")
if (!file.exists(fbc) || !file.exists(frd))
  stop("Need stage-07 outputs for both BC and Rods: run `07` and `CELLTYPE=Rods ... 07` first.")

bc <- read.csv(fbc); bc$cell_type <- "Bipolar cells (real subtypes)"
rd <- read.csv(frd); rd$cell_type <- "Rods (no established subtypes)"
d  <- bind_rows(bc, rd)
cols <- c("Bipolar cells (real subtypes)"="#1B9E77", "Rods (no established subtypes)"="#D95F02")

pa <- ggplot(d, aes(k, mean_cross_seed_ari, color=cell_type)) +
  geom_ribbon(data=subset(d, !is.na(min_cross_seed_ari)),
              aes(ymin=min_cross_seed_ari, ymax=pmin(1, 2*mean_cross_seed_ari-min_cross_seed_ari), fill=cell_type),
              alpha=0.12, color=NA, inherit.aes=TRUE) +
  geom_line(linewidth=1) + geom_point(size=2) +
  scale_color_manual(values=cols, name=NULL) + scale_fill_manual(values=cols, name=NULL, guide="none") +
  scale_x_continuous(breaks=seq(2,18,2)) + ylim(0,1) +
  labs(x="cluster number (k)", y="mean cross-seed ARI",
       title="a  Across-embedding reproducibility vs granularity",
       subtitle="rises toward the true subtype count for real sub-structure only") +
  pub_theme() + theme(legend.position="top")

pb <- ggplot(d, aes(k, silhouette, color=cell_type)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  scale_color_manual(values=cols, name=NULL) + scale_x_continuous(breaks=seq(2,18,2)) +
  labs(x="cluster number (k)", y="silhouette width",
       title="b  Cluster quality vs granularity",
       subtitle="improves with k only when real sub-structure is present") +
  pub_theme() + theme(legend.position="none")

fig <- pa / pb + plot_layout(heights=c(1,0.9))
save_fig(fig, "Figure_5_discovery_comparison", file.path(PIPELINE_DIR,"manuscript","figures"), width=8, height=8)
message("Saved Figure 5 (discovery comparison).")
