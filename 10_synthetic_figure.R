# ============================================================
# 10_synthetic_figure.R   (pipeline_v2)
#
# Figure 6: ground-truth validation on synthetic (Splatter) data.
# (a) UMAPs of representative simulated datasets colored by the TRUE
#     group (null = one population; structured = three groups; rare =
#     a 5% subpopulation) — what the pipeline is asked to detect.
# (b) Recovery of the true partition (ARI vs ground truth) vs k.
# (c) Across-embedding reproducibility (cross-seed ARI) vs k.
# Structured data peak at the true k (=3); null data recover nothing.
#
# Input:  output/tables/11_synthetic_umap_coords.csv (stage 11)
#         output/tables/09_synthetic_validation.csv  (stage 09)
# Output: manuscript/figures/Figure_6_synthetic_validation.(pdf|png)
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

fu <- file.path(TABLE_DIR, "11_synthetic_umap_coords.csv")
fv <- file.path(TABLE_DIR, "09_synthetic_validation.csv")
if (!file.exists(fv)) stop("Run stage 09 first.")

# ---- (a) UMAPs of the simulated data ----
gcols <- c("group 1"="#B0B7BE","group 2"="#3288BD","group 3"="#E7298A")
if (file.exists(fu)) {
  um <- read.csv(fu)
  um$scenario <- factor(um$scenario, levels=unique(um$scenario[order(um$order)]))
  pa <- ggplot(um, aes(UMAP1, UMAP2, color=group)) +
    geom_point(size=0.25, alpha=0.6) +
    facet_wrap(~scenario, nrow=1, scales="free") +
    scale_color_manual(values=gcols, name="true group", guide=guide_legend(override.aes=list(size=2.5,alpha=1))) +
    labs(title="a  Synthetic datasets (colored by ground-truth group)") +
    pub_theme() + theme(axis.text=element_blank(), axis.ticks=element_blank(),
      axis.title=element_blank(), legend.position="right",
      panel.border=element_rect(color="grey85", fill=NA), strip.text=element_text(size=9))
} else { pa <- NULL }

# ---- (b,c) recovery + reproducibility curves ----
d <- read.csv(fv)
lab <- c(null_A="Null (no structure)", null_B="Null (no structure)",
         weak_3="Weak (de=0.03)", medium_3="Medium (de=0.10)",
         strong_3="Strong (de=0.30)", rare_3="Rare 5% subpop (de=0.15)")
d$label <- factor(lab[d$scenario],
  levels=c("Null (no structure)","Weak (de=0.03)","Rare 5% subpop (de=0.15)","Medium (de=0.10)","Strong (de=0.30)"))
cols <- c("Null (no structure)"="grey60","Weak (de=0.03)"="#66C2A5","Rare 5% subpop (de=0.15)"="#E7298A",
          "Medium (de=0.10)"="#3288BD","Strong (de=0.30)"="#1B9E77")

pb <- ggplot(d, aes(k, ari_vs_truth, color=label, group=scenario)) +
  geom_vline(xintercept=3, linetype="dashed", color="grey70") +
  geom_line(linewidth=0.9, alpha=0.9) + geom_point(size=1.8) +
  scale_color_manual(values=cols, name=NULL, guide=guide_legend(nrow=1)) + scale_x_continuous(breaks=2:6) + ylim(0,1) +
  labs(x="cluster number (k)", y="ARI vs. ground truth",
       title="b  Recovery of true structure",
       subtitle="structured data peak at the true k = 3; null data recover nothing") +
  pub_theme() + theme(legend.text=element_text(size=8))

pc <- ggplot(d, aes(k, mean_cross_seed_ari, color=label, group=scenario)) +
  geom_vline(xintercept=3, linetype="dashed", color="grey70") +
  geom_line(linewidth=0.9, alpha=0.9) + geom_point(size=1.8) +
  scale_color_manual(values=cols, name=NULL, guide=guide_legend(nrow=1)) + scale_x_continuous(breaks=2:6) + ylim(0,1) +
  labs(x="cluster number (k)", y="across-embedding cross-seed ARI",
       title="c  Across-embedding reproducibility", subtitle="peaks at the true k for real structure") +
  pub_theme()

# Collect b/c into a single shared legend spanning the full width of the row
# (a per-panel top legend overflowed the half-width panel and clipped "Null").
bc <- (pb | pc) + plot_layout(guides="collect") &
      theme(legend.position="bottom", legend.text=element_text(size=8))
fig <- if (!is.null(pa)) pa / bc + plot_layout(heights=c(0.8,1)) else bc
save_fig(fig, "Figure_6_synthetic_validation", file.path(PIPELINE_DIR,"manuscript","figures"), width=11, height=9)
message("Saved Figure 6 (synthetic validation: data UMAPs + recovery).")
