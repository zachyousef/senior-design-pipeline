# ============================================================
# 06_publication_figures.R   (Final Pipeline)
#
# Generate ALL publication figures for the manuscript from cached
# pipeline outputs. Five main figures:
#
#   Figure 1 — QC + integrated clusters + canonical markers + cell-type ID
#              (combines previous Figure_2_QC_Cluster_Markers and
#              Figure_3_CellTypes_Composition)
#   Figure 2 — Phase 1 positive control: full-dataset benchmark heatmap
#              + chosen-partition evidence (combines previous Figure_4 + 5)
#   Figure 3 — Phase 2 rod application: rod-subset benchmark heatmap
#              + rod subclusters (combines previous Figure_6 + 7)
#   Figure 4 — IPA functional analysis (carry-over image)
#   Figure 5 — Pipeline overview schematic (carry-over image)
#
# ALL figures are saved to manuscript/figures/ with the simple naming
# scheme Figure_1.{png,pdf} ... Figure_5.{png,pdf}.
#
# Stage 06 is self-sufficient: it recomputes Phase 1 ARI/NMI fresh from
# the cached cluster assignments + current obj$cell_type. This means a
# fresh stage 02 with an updated cell-type override will flow through
# to figures via stage 06 alone, without needing to re-run the long
# Phase 1 benchmark in stage 03.
# ============================================================
source("config.R"); source("utils.R")
suppressPackageStartupMessages({
  library(Seurat); library(uwot)
  library(ggplot2); library(patchwork)
  library(dplyr); library(tidyr)
  library(RColorBrewer)
  library(aricode)
  library(cowplot)
})
set.seed(SEED)

OUT <- "manuscript/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SUPP <- "manuscript/supplemental_figures"
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

T0 <- Sys.time()
log_step <- function(m) message(sprintf("[%s] %s",
                                          format(Sys.time() - T0, digits = 3), m))

save_both <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".png")), p,
          width = w, height = h, units = "in", dpi = 600, bg = "white",
          limitsize = FALSE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
          width = w, height = h, units = "in", bg = "white",
          limitsize = FALSE)
  log_step(sprintf("  Saved %s/%s.{png,pdf}", dir, name))
}

read_lat <- function(p) {
  d <- read.csv(p, check.names = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- d[[1]]; m
}

# ---- Shared aesthetics ----
clean_names <- c("PCA_KMeans"   = "PCA + K-Means",
                 "PCA_GMM"      = "PCA + GMM",
                 "scDHA_KMeans" = "scDHA + K-Means",
                 "scDHA_LL"     = "scDHA + Louvain",
                 "PCA_LL"       = "PCA + Louvain",
                 "scDHA_Leiden" = "scDHA + Leiden",
                 "PCA_Leiden"   = "PCA + Leiden")
method_order <- c("PCA + K-Means", "PCA + GMM",
                   "PCA + Louvain", "PCA + Leiden",
                   "scDHA + K-Means",
                   "scDHA + Louvain", "scDHA + Leiden")
unified_palette <- c("#440154", "#3B528B", "#21908C",
                      "#5DC863", "#FDE725")
ct_pal <- c("Rods" = "#E41A1C", "Cones" = "#377EB8",
             "BC" = "#4DAF4A", "Muller Glia" = "#984EA3",
             "RGC/Neuron" = "#FF7F00", "AC" = "#A65628",
             "RPE" = "#F781BF", "Microglia" = "#999999",
             "Unassigned" = "#CCCCCC")
cl_pal <- setNames(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                       "#A65628","#F781BF","#999999","#66C2A5"),
                    as.character(1:9))

size_order_remap <- function(lab) {
  tbl <- sort(table(lab), decreasing = TRUE)
  remap <- setNames(seq_along(tbl), names(tbl))
  factor(remap[as.character(lab)], levels = seq_along(tbl))
}

# Phase 1 chosen partition (positive control): PCA + Louvain @ k = 8.
POC_CHOSEN_METHOD <- "PCA_LL"
POC_CHOSEN_K      <- 8

# ============================================================
# Load shared inputs
# ============================================================
log_step("\n========== Loading shared inputs ==========")
obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
DefaultAssay(obj) <- "SCT"

ump <- Embeddings(obj, "umap")
d_meta <- data.frame(UMAP1 = ump[, 1], UMAP2 = ump[, 2],
                      sample    = obj$sample,
                      time_label = obj$time_label,
                      cluster   = factor(obj$SCT_snn_res.0.6),
                      cell_type = as.character(obj$cell_type))
d_meta <- d_meta[sample(nrow(d_meta)), ]

sample_pal <- setNames(
  RColorBrewer::brewer.pal(min(8, length(unique(d_meta$sample))), "Set2")[
    seq_along(unique(d_meta$sample))],
  sort(unique(d_meta$sample)))

n_cl <- length(levels(d_meta$cluster))
cluster_pal <- colorRampPalette(
  RColorBrewer::brewer.pal(11, "Spectral"))(n_cl)
names(cluster_pal) <- levels(d_meta$cluster)

# ============================================================
# FIGURE 1 — QC violins + Cluster UMAP + Marker dotplot
#            + Cell-type UMAP + Composition bars
# ============================================================
log_step("\n========== FIGURE 1: QC, clusters, markers, cell types ==========")

# Panel a: QC violins on the analysis dataset (UMI / Genes / % Ribosomal).
# NOTE: cell-level filtering used UMI, gene, % mitochondrial and % hemoglobin
# thresholds upstream (stage 01); the mitochondrial and hemoglobin genes are
# then removed, so those percentages are ~0 on the saved object and cannot be
# shown here. % ribosomal remains informative and is displayed instead.
qc_long <- obj@meta.data %>%
  select(sample, nCount_RNA, nFeature_RNA, percent.rb) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("nCount_RNA", "nFeature_RNA", "percent.rb"),
                          labels = c("UMI Counts", "Genes Detected", "% Ribosomal")))

f1a <- ggplot(qc_long, aes(sample, value, fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.size = 0.15, alpha = 0.55, fill = "white") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = sample_pal, guide = "none") +
  pub_theme(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text  = element_text(face = "bold", size = 11),
        plot.title  = element_text(face = "bold", size = 13),
        plot.margin = margin(2, 4, 2, 4, "mm")) +
  labs(title = "a", x = NULL, y = NULL)

# Panel b: Cluster UMAP
f1b <- ggplot(d_meta, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 0.18, alpha = 0.7) +
  scale_color_manual(values = cluster_pal, name = "Cluster", drop = FALSE) +
  pub_theme(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 9),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.30, "cm"),
        plot.margin = margin(2, 4, 2, 4, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1),
                               ncol = 2)) +
  labs(title = "b", x = "UMAP1", y = "UMAP2")

# Panel c: Marker dotplot
ct_levels_dp <- c("Rods", "Cones", "BC", "Muller Glia",
                   "RGC/Neuron", "AC", "RPE", "Microglia")
markers_use <- c("Rho", "Nrl", "Pde6h", "Opn1sw",
                  "Vsx2", "Isl1", "Snhg11", "Slc32a1",
                  "Pax6", "Vim", "Apoe", "Rpe65", "Tmem119")
markers_use <- intersect(markers_use, rownames(obj))

obj_dp <- obj
obj_dp@meta.data$cell_type_dp <-
  factor(as.character(obj_dp$cell_type), levels = rev(ct_levels_dp))
Idents(obj_dp) <- "cell_type_dp"
obj_dp <- subset(obj_dp, idents = ct_levels_dp)
DefaultAssay(obj_dp) <- "SCT"

f1c <- DotPlot(obj_dp, features = markers_use, assay = "SCT",
                cols = c("lightgrey", "#3B287C"), dot.scale = 5) +
  RotatedAxis() +
  pub_theme(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y   = element_text(size = 9),
        axis.title.x  = element_text(face = "bold", size = 10),
        axis.title.y  = element_text(face = "bold", size = 10),
        legend.position = "right",
        legend.title  = element_text(face = "bold", size = 9),
        legend.text   = element_text(size = 8),
        legend.box    = "vertical",
        plot.title    = element_text(face = "bold", size = 13),
        plot.margin   = margin(2, 4, 2, 4, "mm")) +
  labs(title = "c", x = "Marker Genes", y = "Cell Type")

# Panel d: Cell-type UMAP
ct_levels <- c("Rods", "Cones", "BC", "Muller Glia",
                "RGC/Neuron", "AC", "RPE", "Microglia", "Unassigned")
d_meta$cell_type <- factor(d_meta$cell_type, levels = ct_levels)

f1d <- ggplot(d_meta, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.20, alpha = 0.75) +
  scale_color_manual(values = ct_pal, drop = FALSE, name = "Cell Type") +
  pub_theme(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 10),
        legend.text     = element_text(size = 9),
        plot.margin     = margin(2, 4, 2, 4, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "d", x = "UMAP1", y = "UMAP2")

# Panel e: Composition with Overall column
prop_df <- obj@meta.data %>%
  count(sample, cell_type) %>%
  group_by(sample) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
overall <- obj@meta.data %>%
  count(cell_type) %>%
  mutate(sample = "Overall", pct = 100 * n / sum(n))
prop_full <- bind_rows(prop_df, overall) %>%
  mutate(cell_type = factor(cell_type, levels = ct_levels),
         sample = factor(sample,
                          levels = c(sort(unique(obj$sample)), "Overall")))
n_samples <- length(unique(obj$sample))
dash_x <- n_samples + 0.5

f1e <- ggplot(prop_full, aes(sample, pct, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.78) +
  geom_vline(xintercept = dash_x, linetype = "dashed",
              color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = ct_pal, drop = FALSE, name = "Cell Type") +
  scale_y_continuous(expand = c(0, 0)) +
  pub_theme(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        plot.title   = element_text(face = "bold", size = 13),
        legend.title = element_text(face = "bold", size = 10),
        legend.text  = element_text(size = 9),
        plot.margin  = margin(2, 4, 2, 4, "mm")) +
  labs(title = "e", x = NULL, y = "Proportion (%)")

# Assemble Figure 1
fig1 <- f1a / (f1b | f1c) / (f1d | f1e) +
  plot_layout(heights = c(0.85, 1.35, 1.10))
save_both(fig1, OUT, "Figure_1", 14, 16)

# ============================================================
# FIGURE 2 — Phase 1 positive control: heatmap + evidence
# ============================================================
log_step("\n========== FIGURE 2: Phase 1 (heatmap + evidence) ==========")

poc4 <- read.csv("output_full_dataset/tables/poc_scores_matrix_4metric.csv",
                  stringsAsFactors = FALSE)
pocR <- read.csv("output_full_dataset/tables/poc_scores_matrix.csv",
                  stringsAsFactors = FALSE)
poc <- merge(poc4, pocR[, c("method","target_k","silhouette")],
              by = c("method","target_k"))

# Mask cells where the resolution sweep failed (target_k != actual_k)
poc_miss <- !is.na(poc$actual_k) & poc$actual_k != poc$target_k
metric_cols_poc <- c("rogue", "balance", "jacc_min", "marker_spec",
                      "silhouette", "ari", "nmi")
for (mc in metric_cols_poc) {
  if (mc %in% colnames(poc)) poc[poc_miss, mc] <- NA
}
log_step(sprintf("  Masked %d / %d Phase-1 cells (target_k != actual_k)",
                  sum(poc_miss), nrow(poc)))

# ---- Load Phase 1 partitions (for the UMAP / confusion-matrix panels) ----
parts_poc <- read.csv("output_full_dataset/tables/poc_all_cluster_assignments.csv",
                       check.names = FALSE)
gt <- as.character(obj$cell_type); names(gt) <- colnames(obj)
ids <- parts_poc$cell_id
gt_aligned <- gt[ids]

# ---- Phase 1 ARI/NMI come straight from the benchmark table ----
# (poc_scores_matrix_4metric.csv). These are the canonical scored values
# reported in the manuscript; do NOT recompute here, because re-scoring
# against a later copy of obj$cell_type introduces small annotation drift
# and would desynchronise the figure from the benchmark table / text.
log_step(sprintf("  Using benchmark-table ARI/NMI (PCA_LL k=%d ARI = %.4f)",
                  POC_CHOSEN_K,
                  suppressWarnings(poc$ari[poc$method == "PCA_LL" &
                                           poc$target_k == POC_CHOSEN_K][1])))

# ---- Heatmap panels ----
mk_panel_poc <- function(df, metric_col, title_str) {
  d <- df %>%
    mutate(method_label = factor(clean_names[method],
                                  levels = rev(method_order)),
           k_label = factor(paste0("k=", target_k),
                             levels = paste0("k=",
                                              sort(unique(target_k)))),
           is_chosen = method == POC_CHOSEN_METHOD &
                        target_k == POC_CHOSEN_K,
           value = .data[[metric_col]],
           label = ifelse(is.na(value), "—",
                           sprintf("%.2f", value)))
  vmin <- min(d$value, na.rm = TRUE)
  vmax <- max(d$value, na.rm = TRUE)
  d$value_norm <- (d$value - vmin) / max((vmax - vmin), 1e-9)
  d$txt_col <- ifelse(!is.na(d$value_norm) & d$value_norm > 0.55,
                       "grey10", "white")
  ggplot(d, aes(k_label, method_label)) +
    geom_tile(aes(fill = value), color = "grey20", linewidth = 0.3) +
    geom_tile(data = subset(d, is_chosen), fill = NA,
              color = "black", linewidth = 1.4) +
    geom_text(aes(label = label, color = txt_col),
              size = 3.4, fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradientn(colors = unified_palette, na.value = "grey85",
                          name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    pub_theme(base_size = 11) +
    theme(axis.title  = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "bottom",
          legend.key.height = unit(0.30, "cm"),
          legend.key.width  = unit(1.20, "cm"),
          legend.text = element_text(size = 9),
          panel.grid = element_blank(),
          plot.title = element_text(size = 13, face = "bold",
                                       hjust = 0.5, margin = margin(b = 3)),
          plot.margin = margin(1, 1, 1, 1, "mm")) +
    labs(title = title_str)
}

p1_a <- mk_panel_poc(poc, "rogue",       "e — ROGUE purity (baseline 0.47)")
p1_b <- mk_panel_poc(poc, "balance",     "f — Cluster-size balance (1 - Gini)")
p1_c <- mk_panel_poc(poc, "jacc_min",    "g — Per-cluster Jaccard min stability")
p1_d <- mk_panel_poc(poc, "marker_spec", "h — Marker specificity")
p1_e <- mk_panel_poc(poc, "silhouette",  "i — Silhouette width")
p1_f <- mk_panel_poc(poc, "ari",         "j — ARI vs cell-type ground truth")
p1_g <- mk_panel_poc(poc, "nmi",         "k — NMI vs cell-type ground truth")

# ---- Evidence panels ----
ump_pca <- read_lat("output_full_dataset/tables/poc_harmony_pca_umap.csv")
ump_pca <- ump_pca[ids, , drop = FALSE]
raw_pll8 <- parts_poc$PCA_LL_k8; names(raw_pll8) <- parts_poc$cell_id
lab_pll8 <- size_order_remap(raw_pll8[ids])
gt_v <- factor(gt_aligned, levels = ct_levels)

mk_umap <- function(ump, grp, palette, title, subtitle) {
  d <- data.frame(UMAP1 = ump[, 1], UMAP2 = ump[, 2], group = grp)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = group)) +
    geom_point(size = 0.30, alpha = 0.75) +
    scale_color_manual(values = palette, drop = FALSE) +
    pub_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey25"),
          legend.title = element_blank(),
          legend.text  = element_text(size = 8),
          legend.position = "bottom",
          legend.box = "horizontal",
          plot.margin = margin(2, 2, 2, 2, "mm")) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                 nrow = 1)) +
    labs(title = title, subtitle = subtitle)
}

f2a <- mk_umap(ump_pca, lab_pll8, cl_pal, "a",
                "Pipeline partition: PCA + Louvain @ k = 8")
f2b <- mk_umap(ump_pca, gt_v, ct_pal, "b",
                "Manual cell-type annotation (ground truth)")

# Confusion matrix
ct_v <- gt_aligned
tab <- table(GroundTruth = ct_v, Cluster = lab_pll8)
pct <- as.data.frame(round(100 * t(t(tab) / colSums(tab)), 1))
pct$Cluster <- factor(pct$Cluster, levels = sort(unique(pct$Cluster)))
pct$GroundTruth <- factor(pct$GroundTruth, levels = ct_levels)

f2c <- ggplot(pct, aes(Cluster, GroundTruth, fill = Freq)) +
  geom_tile(color = "grey50", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", Freq),
                  color = ifelse(Freq > 50, "white", "grey20")),
            size = 2.8, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradient(low = "white", high = "#3B528B",
                       limits = c(0, 100), name = "% of cluster") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_equal() +
  pub_theme(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 9, margin = margin(t = 3)),
        axis.title.y = element_text(face = "bold", size = 9, margin = margin(r = 3)),
        legend.position = "right",
        legend.key.height = unit(1.0, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  labs(title = "c — Confusion matrix",
       x = "Pipeline cluster (1 = largest)", y = "Cell type (ground truth)")

# ARI bar chart at chosen k=8
ari_at_k <- poc %>%
  filter(target_k == POC_CHOSEN_K, !is.na(ari)) %>%
  mutate(method_label = factor(clean_names[method], levels = method_order))

f2d <- ggplot(ari_at_k,
               aes(reorder(method_label, ari), ari, fill = method_label)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("ARI = %.3f", ari)),
            hjust = -0.05, size = 3.0, fontface = "bold") +
  geom_hline(yintercept = 0.95, linetype = "dashed",
              color = "darkgreen", linewidth = 0.5) +
  annotate("text", x = 0.7, y = 0.95, label = "perfect recovery",
            color = "darkgreen", size = 2.6, vjust = -0.3, hjust = 0) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.18), expand = c(0, 0),
                      breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2, "mm")) +
  labs(title = sprintf("d — Cell-type recovery at chosen k = %d", POC_CHOSEN_K),
       y = sprintf("ARI vs cell-type annotation (k = %d)", POC_CHOSEN_K))

# Re-label panel titles (heatmap a-g, evidence h-k)
p1_a <- p1_a + labs(title = "a — ROGUE purity (baseline 0.47)")
p1_b <- p1_b + labs(title = "b — Cluster-size balance (1 - Gini)")
p1_c <- p1_c + labs(title = "c — Per-cluster Jaccard min stability")
p1_d <- p1_d + labs(title = "d — Marker specificity")
p1_f <- p1_f + labs(title = "e — ARI vs cell-type ground truth")
p1_g <- p1_g + labs(title = "f — NMI vs cell-type ground truth")
p1_e <- p1_e + labs(title = "g — Silhouette width")
f2a  <- f2a  + labs(title = "h", subtitle = "Pipeline partition: PCA + Louvain @ k = 8")
f2b  <- f2b  + labs(title = "i",
                      subtitle = "Manual cell-type annotation (ground truth)")
f2c  <- f2c  + labs(title = "j — Confusion matrix",
                      x = "Pipeline cluster (1 = largest)",
                      y = "Cell type (ground truth)")
f2d  <- f2d  + labs(title = sprintf("k — Cell-type recovery at chosen k = %d",
                                        POC_CHOSEN_K),
                      y = sprintf("ARI vs cell-type annotation (k = %d)",
                                    POC_CHOSEN_K))

# Build Figure 2 with cowplot::plot_grid throughout for rigid
# equal-cell sizing. Patchwork's design-string + widths argument
# wasn't enough to keep right-column heatmaps the same width as the
# left column (axis-label/legend chrome leaks into col-width
# allocation). cowplot allocates each grid cell strictly by
# rel_widths regardless of the panel's chrome.
hm_row1 <- cowplot::plot_grid(p1_a, p1_b, ncol = 2, align = "hv", axis = "tblr")
hm_row2 <- cowplot::plot_grid(p1_c, p1_d, ncol = 2, align = "hv", axis = "tblr")
hm_row3 <- cowplot::plot_grid(p1_f, p1_g, ncol = 2, align = "hv", axis = "tblr")
hm_row4 <- cowplot::plot_grid(NULL, p1_e, NULL,
                                ncol = 3, rel_widths = c(1, 2, 1))
heatmap_block <- cowplot::plot_grid(hm_row1, hm_row2, hm_row3, hm_row4,
                                      ncol = 1, align = "v", axis = "lr")

ev_row1 <- cowplot::plot_grid(f2a, f2b, ncol = 2, align = "hv", axis = "tblr")
ev_row2 <- cowplot::plot_grid(f2c, f2d, ncol = 2, align = "hv", axis = "tblr")
evidence_block <- cowplot::plot_grid(ev_row1, ev_row2, ncol = 1)

fig2 <- cowplot::plot_grid(
  heatmap_block, evidence_block,
  ncol = 1, rel_heights = c(26, 16)
)
save_both(fig2, OUT, "Figure_2", 16, 42)

# ============================================================
# FIGURE 3 — Phase 2 rod application: subclusters + heatmap
# ============================================================
log_step("\n========== FIGURE 3: Phase 2 (subclusters + heatmap) ==========")

rod4 <- read.csv("output/tables/03_scores_matrix_4metric.csv",
                  stringsAsFactors = FALSE)
rodR <- read.csv("output/tables/03_scores_matrix.csv",
                  stringsAsFactors = FALSE)
rod <- merge(rod4, rodR[, c("method","target_k","silhouette")],
              by = c("method","target_k"))
rod_miss <- !is.na(rod$actual_k) & rod$actual_k != rod$target_k
metric_cols_rod <- c("rogue", "balance", "jacc_min", "marker_spec", "silhouette")
for (mc in metric_cols_rod) {
  if (mc %in% colnames(rod)) rod[rod_miss, mc] <- NA
}
log_step(sprintf("  Masked %d / %d Phase-2 cells (target_k != actual_k)",
                  sum(rod_miss), nrow(rod)))

mk_panel_rod <- function(df, metric_col, title_str) {
  d <- df %>%
    mutate(method_label = factor(clean_names[method],
                                  levels = rev(method_order)),
           k_label = factor(paste0("k=", target_k),
                             levels = paste0("k=",
                                              sort(unique(target_k)))),
           is_chosen = method == CHOSEN_METHOD & target_k == CHOSEN_K,
           value = .data[[metric_col]],
           label = ifelse(is.na(value), "—",
                           sprintf("%.2f", value)))
  vmin <- min(d$value, na.rm = TRUE)
  vmax <- max(d$value, na.rm = TRUE)
  d$value_norm <- (d$value - vmin) / max((vmax - vmin), 1e-9)
  d$txt_col <- ifelse(!is.na(d$value_norm) & d$value_norm > 0.55,
                       "grey10", "white")
  ggplot(d, aes(k_label, method_label)) +
    geom_tile(aes(fill = value), color = "grey20", linewidth = 0.3) +
    geom_tile(data = subset(d, is_chosen), fill = NA,
              color = "black", linewidth = 1.6) +
    geom_text(aes(label = label, color = txt_col),
              size = 3.8, fontface = "bold") +
    scale_color_identity() +
    scale_fill_gradientn(colors = unified_palette, na.value = "grey75",
                          name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_equal() +
    pub_theme(base_size = 12) +
    theme(axis.title = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 11),
          legend.position = "bottom",
          legend.key.height = unit(0.30, "cm"),
          legend.key.width  = unit(1.20, "cm"),
          legend.text = element_text(size = 10),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14, face = "bold",
                                       hjust = 0.5, margin = margin(b = 3)),
          plot.margin = margin(1, 1, 1, 1, "mm")) +
    labs(title = title_str)
}

p2_a <- mk_panel_rod(rod, "rogue",       "h — ROGUE purity")
p2_b <- mk_panel_rod(rod, "balance",     "i — Cluster-size balance (1 - Gini)")
p2_c <- mk_panel_rod(rod, "jacc_min",    "j — Per-cluster Jaccard min stability")
p2_d <- mk_panel_rod(rod, "marker_spec", "k — Marker specificity")
p2_e <- mk_panel_rod(rod, "silhouette",  "l — Silhouette width")

# rod_heatmap_block built inline below in the Fig 3 assembly

# ---- Rod subcluster panels ----
har_rod <- read_lat("output/tables/03_harmony_scdha_latent.csv")
parts_rod <- read.csv("output/tables/03_all_cluster_assignments.csv",
                       check.names = FALSE)
rods <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
DefaultAssay(rods) <- "RNA"
rods[["RNA"]] <- JoinLayers(rods[["RNA"]])
rod_ids <- parts_rod$cell_id
rods_clean <- subset(rods, cells = rod_ids); rods_clean <- rods_clean[, rod_ids]
rm(rods); gc()

rod_ump_path <- "output/tables/03_harmony_scdha_umap.csv"
if (file.exists(rod_ump_path)) {
  ump_rod <- read_lat(rod_ump_path)[rod_ids, , drop = FALSE]
} else {
  set.seed(SEED)
  ump_rod <- uwot::umap(har_rod[rod_ids, ],
                          n_neighbors = 30, min_dist = 0.3,
                          n_threads = N_CORES, verbose = FALSE)
  rownames(ump_rod) <- rod_ids; colnames(ump_rod) <- c("UMAP1","UMAP2")
  write.csv(data.frame(cell_id = rownames(ump_rod), ump_rod),
            rod_ump_path, row.names = FALSE)
}

partition_col <- paste0(CHOSEN_METHOD, "_k", CHOSEN_K)
raw_lab <- parts_rod[[partition_col]]; names(raw_lab) <- parts_rod$cell_id
lab_v <- size_order_remap(raw_lab[rod_ids])
time_v <- rods_clean$time_label[rod_ids]
xlim_all <- range(ump_rod[, 1]) + c(-0.5, 0.5)
ylim_all <- range(ump_rod[, 2]) + c(-0.5, 0.5)

mk_rod_umap_big <- function(idx, lab, title) {
  d <- data.frame(UMAP1 = ump_rod[idx, 1], UMAP2 = ump_rod[idx, 2],
                   cluster = lab)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(size = 0.45, alpha = 0.8) +
    scale_color_manual(values = cl_pal, drop = FALSE, name = "Subcluster") +
    coord_fixed(xlim = xlim_all, ylim = ylim_all) +
    pub_theme(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          legend.position = "right",
          legend.title  = element_text(face = "bold", size = 10),
          legend.text   = element_text(size = 9),
          plot.margin   = margin(3, 3, 3, 3, "mm")) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    labs(title = title, x = "UMAP1", y = "UMAP2")
}

mk_rod_umap_small <- function(idx, lab, title) {
  d <- data.frame(UMAP1 = ump_rod[idx, 1], UMAP2 = ump_rod[idx, 2],
                   cluster = lab)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(size = 0.30, alpha = 0.75) +
    scale_color_manual(values = cl_pal, drop = FALSE, guide = "none") +
    coord_fixed(xlim = xlim_all, ylim = ylim_all) +
    pub_theme(base_size = 8) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          axis.title    = element_text(size = 7),
          axis.text     = element_text(size = 6),
          plot.margin   = margin(2, 2, 2, 2, "mm")) +
    labs(title = title, x = "UMAP1", y = "UMAP2")
}

cluster_sizes <- as.numeric(table(lab_v))
n_clusters <- length(cluster_sizes)

f3a <- mk_rod_umap_big(seq_along(rod_ids), lab_v,
   sprintf("a — Rod Photoreceptor Subclusters (k = %d)", CHOSEN_K))

tps <- intersect(c("P15", "P49", "P75", "P91"), unique(time_v))
tp_panels <- lapply(seq_along(tps), function(i) {
  tp <- tps[i]
  idx <- which(time_v == tp)
  mk_rod_umap_small(idx, lab_v[idx], sprintf("%s — %s", letters[i + 1], tp))
})

# Per-subcluster DEG counts (compute or load cached)
deg_path <- file.path(TABLE_DIR, "05_deg_counts.csv")
if (!file.exists(deg_path)) deg_path <- file.path(TABLE_DIR, "04_deg_counts.csv")
if (file.exists(deg_path)) {
  deg_counts <- read.csv(deg_path, stringsAsFactors = FALSE)
  deg_counts$cluster <- as.character(deg_counts$cluster)
} else {
  log_step("  Computing FindAllMarkers inline (no cached DEG counts)")
  rcd <- rods_clean
  DefaultAssay(rcd) <- "RNA"
  rcd[["RNA"]] <- JoinLayers(rcd[["RNA"]])
  rcd <- NormalizeData(rcd, verbose = FALSE)
  cl_ids <- as.character(lab_v)
  names(cl_ids) <- rod_ids
  rcd@meta.data$cluster_v <-
    factor(cl_ids[colnames(rcd)],
            levels = as.character(seq_along(cluster_sizes)))
  Idents(rcd) <- "cluster_v"
  markers_inline <- FindAllMarkers(rcd, only.pos = TRUE,
                                    min.pct = 0.10, logfc.threshold = 0.25,
                                    verbose = FALSE)
  markers_sig <- markers_inline %>% filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))
  deg_counts <- markers_sig %>% count(cluster, name = "n_sig_markers") %>%
    mutate(cluster = as.character(cluster)) %>% arrange(cluster)
  write.csv(deg_counts, file.path(TABLE_DIR, "05_deg_counts.csv"),
            row.names = FALSE)
  write.csv(markers_sig, file.path(TABLE_DIR, "05_all_markers.csv"),
            row.names = FALSE)
  rm(rcd); gc()
}
deg_full <- data.frame(cluster = as.character(seq_len(n_clusters)))
deg_full <- merge(deg_full, deg_counts, by = "cluster", all.x = TRUE)
deg_full$n_sig_markers[is.na(deg_full$n_sig_markers)] <- 0
deg_full$cluster <- factor(deg_full$cluster,
                            levels = as.character(seq_len(n_clusters)))

f3_size <- ggplot(
  data.frame(cluster = factor(seq_along(cluster_sizes)),
              n = cluster_sizes,
              pct = round(100 * cluster_sizes / length(rod_ids), 1)),
  aes(cluster, n, fill = cluster)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%s\n(%s%%)",
                                 format(n, big.mark = ","), pct)),
            vjust = -0.25, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cl_pal, guide = "none") +
  pub_theme(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.margin = margin(2, 3, 2, 3, "mm")) +
  expand_limits(y = max(cluster_sizes) * 1.22) +
  labs(title = "f — Subcluster sizes",
       x = "Subcluster", y = "Number of cells")

f3_deg <- ggplot(deg_full, aes(cluster, n_sig_markers, fill = cluster)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = format(n_sig_markers, big.mark = ",")),
            vjust = -0.25, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cl_pal, guide = "none") +
  pub_theme(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.margin = margin(2, 3, 2, 3, "mm")) +
  expand_limits(y = max(deg_full$n_sig_markers) * 1.22) +
  labs(title = "g — Per-subcluster marker counts",
       x = "Subcluster", y = "Significant cluster-defining DEGs")

# Re-label panels (heatmap a-e, subclusters f-l)
p2_a <- p2_a + labs(title = "a — ROGUE purity")
p2_b <- p2_b + labs(title = "b — Cluster-size balance (1 - Gini)")
p2_c <- p2_c + labs(title = "c — Per-cluster Jaccard min stability")
p2_d <- p2_d + labs(title = "d — Marker specificity")
p2_e <- p2_e + labs(title = "e — Silhouette width")
f3a  <- f3a  + labs(title = sprintf("f — Rod Photoreceptor Subclusters (k = %d)",
                                        CHOSEN_K))
tp_panels[[1]] <- tp_panels[[1]] + labs(title = sprintf("g — %s", tps[1]))
tp_panels[[2]] <- tp_panels[[2]] + labs(title = sprintf("h — %s", tps[2]))
tp_panels[[3]] <- tp_panels[[3]] + labs(title = sprintf("i — %s", tps[3]))
tp_panels[[4]] <- tp_panels[[4]] + labs(title = sprintf("j — %s", tps[4]))
f3_size <- f3_size + labs(title = "k — Subcluster sizes")
f3_deg  <- f3_deg  + labs(title = "l — Per-subcluster marker counts")

# Build Figure 3 with cowplot::plot_grid throughout (same reason as
# Figure 2 — keep heatmap columns rigidly equal width).
hm2_row1 <- cowplot::plot_grid(p2_a, p2_b, ncol = 2, align = "hv", axis = "tblr")
hm2_row2 <- cowplot::plot_grid(p2_c, p2_d, ncol = 2, align = "hv", axis = "tblr")
hm2_row3 <- cowplot::plot_grid(NULL, p2_e, NULL,
                                 ncol = 3, rel_widths = c(1, 2, 1))
heatmap_block_p2 <- cowplot::plot_grid(hm2_row1, hm2_row2, hm2_row3,
                                         ncol = 1, align = "v", axis = "lr")

mid_grid <- cowplot::plot_grid(
  tp_panels[[1]], tp_panels[[2]],
  tp_panels[[3]], tp_panels[[4]],
  ncol = 2)
right_col <- cowplot::plot_grid(f3_size, f3_deg, ncol = 1)
subcluster_block <- cowplot::plot_grid(
  f3a, mid_grid, right_col,
  ncol = 3, rel_widths = c(2.4, 2.0, 1.1)
)

fig3 <- cowplot::plot_grid(
  heatmap_block_p2, subcluster_block,
  ncol = 1, rel_heights = c(20, 8)
)
save_both(fig3, OUT, "Figure_3", 16, 28)

# ============================================================
# FIGURE 4 — IPA functional analysis (carry-over)
# ============================================================
log_step("\n========== FIGURE 4: IPA functional ==========")
ipa_src <- "assets/IPA_functional_main.png"
if (file.exists(ipa_src)) {
  file.copy(ipa_src, file.path(OUT, "Figure_4.png"), overwrite = TRUE)
  log_step(sprintf("  Copied IPA placeholder from %s", ipa_src))
} else {
  log_step("  WARNING: assets/IPA_functional_main.png not found; Figure 4 not generated")
}

# ============================================================
# FIGURE 5 — Pipeline schematic
# ============================================================
log_step("\n========== FIGURE 5: Pipeline overview ==========")
pipeline_src <- "assets/pipeline_schematic.png"
if (file.exists(pipeline_src)) {
  file.copy(pipeline_src, file.path(OUT, "Figure_5.png"), overwrite = TRUE)
  log_step(sprintf("  Copied pipeline schematic from %s", pipeline_src))
} else {
  log_step("  WARNING: assets/pipeline_schematic.png not found; Figure 5 not generated")
}

# ============================================================
# SUPPLEMENTARY FIGURES & TABLES
# ============================================================
log_step("\n========== SUPPLEMENTAL: figures + tables ==========")

# ------------------------------------------------------------
# Supplemental Figure S1 — Extended QC + Harmony batch validation
# ------------------------------------------------------------
log_step("Supplemental_Fig_S1: extended QC + by-sample UMAP")

# (a) cells per sample after filtering
sa_cells <- obj@meta.data %>% count(sample, name = "n_cells")
sa_cells$sample <- factor(sa_cells$sample, levels = sort(sa_cells$sample))
s1a <- ggplot(sa_cells, aes(sample, n_cells, fill = sample)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n_cells, big.mark = ",")),
             vjust = -0.3, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = sample_pal, guide = "none") +
  expand_limits(y = max(sa_cells$n_cells) * 1.12) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "a — Cells per sample (post-filter)",
       x = NULL, y = "Cells retained after QC")

# (b) doublet flag counts per sample (singlets retained vs doublets removed)
if ("scDblFinder.class" %in% colnames(obj@meta.data)) {
  dbl_long <- obj@meta.data %>%
    count(sample, scDblFinder.class) %>%
    mutate(scDblFinder.class = factor(scDblFinder.class,
                                        levels = c("singlet", "doublet")))
  s1b <- ggplot(dbl_long, aes(sample, n, fill = scDblFinder.class)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(singlet = "#3B528B", doublet = "#FDE725"),
                       name = "scDblFinder") +
    pub_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "b — Doublet flags per sample",
         x = NULL, y = "Cell count")
} else {
  s1b <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, size = 4,
              label = "scDblFinder.class not in metadata") +
    labs(title = "b — Doublet flags per sample")
}

# (c) Per-sample QC violins on the analysis dataset.
# Mitochondrial and hemoglobin genes are removed during stage 01 (after they
# are used for filtering), so percent.mt/percent.hb are ~0 on the saved object
# and are not shown; UMI, genes, and % ribosomal remain informative.
qc_metrics_long <- obj@meta.data %>%
  select(sample, nCount_RNA, nFeature_RNA, percent.rb) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                          levels = c("nCount_RNA", "nFeature_RNA", "percent.rb"),
                          labels = c("UMI counts", "Genes detected",
                                      "% ribosomal")))

s1c <- ggplot(qc_metrics_long, aes(sample, value, fill = sample)) +
  geom_violin(scale = "width", alpha = 0.85, trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.size = 0.15, alpha = 0.55, fill = "white") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = sample_pal, guide = "none") +
  pub_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "c — Per-sample QC distributions (UMI, genes, % ribosomal)",
       x = NULL, y = NULL)

# (d) Harmony UMAP colored by sample (batch-mixing check)
s1d <- ggplot(d_meta, aes(UMAP1, UMAP2, color = sample)) +
  geom_point(size = 0.18, alpha = 0.6) +
  scale_color_manual(values = sample_pal, name = "Sample") +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "d — Harmony-corrected UMAP, colored by sample",
       x = "UMAP1", y = "UMAP2")

s_fig_s1 <- (s1a | s1b) / s1c / s1d +
  plot_layout(heights = c(0.9, 0.9, 1.4))
save_both(s_fig_s1, SUPP, "Supplemental_Fig_S1", 14, 14)

# ------------------------------------------------------------
# Supplemental Figure S2 — Robustness suite for the chosen rod partition
# ------------------------------------------------------------
log_step("Supplemental_Fig_S2: robustness suite")

knn_sum <- read.csv("output/tables/_knn_robustness_summary.csv",
                     stringsAsFactors = FALSE)
seed_sum_path <- "output/tables/_seed_robustness_summary.csv"
seed_sum <- if (file.exists(seed_sum_path)) {
  read.csv(seed_sum_path, stringsAsFactors = FALSE)
} else NULL

# (a) kNN sensitivity at k=3: cluster sizes across k.param values
knn_k3 <- knn_sum %>% filter(target_k == 3) %>%
  mutate(parts = strsplit(sizes, " / ")) %>%
  rowwise() %>%
  mutate(c1 = as.integer(gsub("[^0-9]", "", parts[1])),
         c2 = as.integer(gsub("[^0-9]", "", parts[2])),
         c3 = as.integer(gsub("[^0-9]", "", parts[3]))) %>%
  ungroup() %>%
  select(k_param, c1, c2, c3) %>%
  pivot_longer(-k_param, names_to = "cluster", values_to = "size") %>%
  mutate(cluster = factor(gsub("c", "", cluster),
                            levels = as.character(seq_len(3))))

s2a <- ggplot(knn_k3, aes(factor(k_param), size, fill = cluster)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = format(size, big.mark = ",")),
             position = position_dodge(width = 0.78),
             vjust = -0.3, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = cl_pal, name = "Subcluster") +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "a — kNN sensitivity (scDHA + Louvain @ k = 3)",
       x = "SNN graph k.param", y = "Cluster size (cells)")

# (b) overlap with default kNN=20 partition
overlap_col <- grep("overlap", colnames(knn_sum), value = TRUE)[1]
knn_overlap <- knn_sum %>% filter(target_k == 3) %>%
  select(k_param, all_of(overlap_col)) %>%
  rename(overlap = !!overlap_col)

s2b <- ggplot(knn_overlap, aes(factor(k_param), overlap)) +
  geom_col(width = 0.7, fill = "#3B528B") +
  geom_text(aes(label = sprintf("%.1f%%", overlap)),
             vjust = -0.3, size = 3.0, fontface = "bold") +
  geom_hline(yintercept = 95, linetype = "dashed", color = "darkgreen") +
  expand_limits(y = 105) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "b — Cell-overlap to default kNN = 20",
       x = "SNN graph k.param",
       y = "% cells with same cluster label")

# (c) per-sample distribution of subclusters
psamp_path <- "output/tables/_per_sample_clusters_scdha_LL.csv"
if (file.exists(psamp_path)) {
  psamp <- read.csv(psamp_path, stringsAsFactors = FALSE)
  # Expect cols: sample, cluster, n, pct (or similar) — adapt if needed
  if (!all(c("sample", "cluster") %in% colnames(psamp))) {
    psamp <- NULL
  }
}
if (is.null(psamp_path) || !file.exists(psamp_path)) psamp <- NULL
if (!is.null(psamp) && nrow(psamp) > 0) {
  psamp$cluster <- factor(as.character(psamp$cluster),
                            levels = as.character(seq_len(9)))
  s2c <- ggplot(psamp, aes(sample, fill = cluster)) +
    {if ("pct" %in% colnames(psamp)) {
        geom_col(aes(y = pct), position = "stack", width = 0.7)
     } else {
        geom_col(aes(y = n), position = "fill", width = 0.7)
     }} +
    scale_fill_manual(values = cl_pal, drop = FALSE, name = "Subcluster") +
    pub_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "c — Per-sample subcluster composition",
         x = NULL, y = "Proportion (%)")
} else {
  s2c <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, size = 4,
              label = "Per-sample subcluster table not found") +
    labs(title = "c — Per-sample subcluster composition")
}

# (d) per-timepoint distribution
ptime_path <- "output/tables/_per_timepoint_clusters_scdha_LL.csv"
if (file.exists(ptime_path)) {
  ptime <- read.csv(ptime_path, stringsAsFactors = FALSE)
  if (all(c("time_label", "cluster") %in% colnames(ptime))) {
    ptime$cluster <- factor(as.character(ptime$cluster),
                              levels = as.character(seq_len(9)))
    s2d <- ggplot(ptime, aes(time_label, fill = cluster)) +
      {if ("pct" %in% colnames(ptime)) {
          geom_col(aes(y = pct), position = "stack", width = 0.7)
       } else {
          geom_col(aes(y = n), position = "fill", width = 0.7)
       }} +
      scale_fill_manual(values = cl_pal, drop = FALSE, name = "Subcluster") +
      pub_theme(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 12)) +
      labs(title = "d — Per-timepoint subcluster composition",
           x = NULL, y = "Proportion (%)")
  } else {
    s2d <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, size = 4,
                label = "Per-timepoint subcluster table not found") +
      labs(title = "d — Per-timepoint subcluster composition")
  }
} else {
  s2d <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, size = 4,
              label = "Per-timepoint subcluster table not found") +
    labs(title = "d — Per-timepoint subcluster composition")
}

s_fig_s2 <- (s2a | s2b) / (s2c | s2d) +
  plot_layout(heights = c(1, 1))
save_both(s_fig_s2, SUPP, "Supplemental_Fig_S2", 14, 10)

# ------------------------------------------------------------
# Supplemental Figure S3 — Phase 1 supporting (resolution plateaus + ARI sweep)
# ------------------------------------------------------------
log_step("Supplemental_Fig_S3: Phase 1 supporting")

# (a) ARI by (method, k) — line plot (alternative view of heatmap)
ari_long <- poc %>%
  filter(!is.na(ari)) %>%
  mutate(method_label = factor(clean_names[method], levels = method_order))

s3a <- ggplot(ari_long, aes(target_k, ari, color = method_label,
                              group = method_label)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  geom_hline(yintercept = 0.95, linetype = "dashed",
              color = "darkgreen", linewidth = 0.4) +
  scale_color_brewer(palette = "Set2", name = "Method") +
  scale_x_continuous(breaks = sort(unique(ari_long$target_k))) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right") +
  labs(title = "a — Phase 1 ARI vs cell-type ground truth, by (method, k)",
       x = "Target k", y = "ARI")

# (b) Resolution plateau diagnostic: target_k vs actual_k for graph methods
plateau <- pocR %>%
  filter(method %in% c("PCA_LL", "PCA_Leiden", "scDHA_LL", "scDHA_Leiden")) %>%
  mutate(method_label = clean_names[method],
         hit = target_k == actual_k)

s3b <- ggplot(plateau, aes(target_k, actual_k, color = method_label,
                             shape = hit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Set2", name = "Method") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                      name = "Hit target_k") +
  scale_x_continuous(breaks = sort(unique(plateau$target_k))) +
  scale_y_continuous(breaks = sort(unique(plateau$actual_k))) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right") +
  labs(title = "b — Resolution-plateau diagnostic",
       x = "Target k", y = "Actual k achieved")

# (c) cluster sizes at chosen partition
sizes_chosen <- as.numeric(table(lab_pll8))
size_df <- data.frame(cluster = factor(seq_along(sizes_chosen)),
                       n = sizes_chosen,
                       pct = round(100 * sizes_chosen / sum(sizes_chosen), 1))

s3c <- ggplot(size_df, aes(cluster, n, fill = cluster)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%s%%)",
                                  format(n, big.mark = ","), pct)),
             vjust = -0.25, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cl_pal, guide = "none") +
  expand_limits(y = max(sizes_chosen) * 1.18) +
  pub_theme(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12)) +
  labs(title = "c — Cluster sizes at chosen Phase 1 partition (PCA + Louvain @ k = 8)",
       x = "Pipeline cluster (1 = largest)", y = "Cells")

s_fig_s3 <- s3a / (s3b | s3c) +
  plot_layout(heights = c(0.9, 1.0))
save_both(s_fig_s3, SUPP, "Supplemental_Fig_S3", 14, 10)

# ------------------------------------------------------------
# Supplemental Figure S4 — IPA extended (regulators + divergent pathways)
# ------------------------------------------------------------
log_step("Supplemental_Fig_S4: IPA extended (carryover placeholder)")
ipa_extended_src <- "assets/IPA_extended.png"
if (file.exists(ipa_extended_src)) {
  file.copy(ipa_extended_src,
             file.path(SUPP, "Supplemental_Fig_S4.png"), overwrite = TRUE)
  log_step(sprintf("  Copied IPA-extended placeholder from %s",
                    ipa_extended_src))
} else {
  log_step("  WARNING: assets/IPA_extended.png not found; placeholder skipped")
}

# ============================================================
# SUPPLEMENTAL TABLES
# ============================================================
log_step("\nSupplemental tables")

# S1: Phase 1 (method, k) full score table
poc_export <- poc %>%
  mutate(method_label = clean_names[method]) %>%
  select(method, method_label, target_k, actual_k,
          rogue, balance, jacc_min, marker_spec, silhouette, ari, nmi,
          any_of("sizes_str"))
write.csv(poc_export, file.path(SUPP, "Supplemental_Table_S1.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (Phase 1 scores, %d rows)",
                  SUPP, "Supplemental_Table_S1.csv", nrow(poc_export)))

# S2: Phase 2 (method, k) full score table
rod_export <- rod %>%
  mutate(method_label = clean_names[method]) %>%
  select(method, method_label, target_k, actual_k,
          rogue, balance, jacc_min, marker_spec, silhouette,
          any_of("sizes_str"))
write.csv(rod_export, file.path(SUPP, "Supplemental_Table_S2.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (Phase 2 scores, %d rows)",
                  SUPP, "Supplemental_Table_S2.csv", nrow(rod_export)))

# S3: Per-cluster cluster-defining DEGs (chosen rod partition)
all_markers_path <- file.path(TABLE_DIR, "05_all_markers.csv")
if (file.exists(all_markers_path)) {
  file.copy(all_markers_path, file.path(SUPP, "Supplemental_Table_S3.csv"),
             overwrite = TRUE)
  ms <- read.csv(all_markers_path)
  log_step(sprintf("  Wrote %s/%s (%d cluster-defining DEGs)",
                    SUPP, "Supplemental_Table_S3.csv", nrow(ms)))
}

# S4: Cell-type composition per sample
ct_per_sample <- obj@meta.data %>%
  count(sample, cell_type) %>%
  group_by(sample) %>%
  mutate(pct = round(100 * n / sum(n), 2)) %>%
  ungroup() %>%
  arrange(sample, desc(n))
write.csv(ct_per_sample, file.path(SUPP, "Supplemental_Table_S4.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (cell-type composition per sample)",
                  SUPP, "Supplemental_Table_S4.csv"))

# S5: Robustness suite (kNN + seed combined)
robust_combined <- knn_sum %>% mutate(test = "kNN")
if (!is.null(seed_sum) && nrow(seed_sum) > 0) {
  seed_sum$test <- "seed"
  shared <- intersect(colnames(seed_sum), colnames(robust_combined))
  robust_combined <- bind_rows(robust_combined[, shared, drop = FALSE],
                                seed_sum[, shared, drop = FALSE])
}
write.csv(robust_combined, file.path(SUPP, "Supplemental_Table_S5.csv"),
          row.names = FALSE)
log_step(sprintf("  Wrote %s/%s (robustness suite: kNN + seed)",
                  SUPP, "Supplemental_Table_S5.csv"))

log_step(sprintf("\n=== ALL FIGURES + SUPPLEMENTAL DONE — %s ===",
                  format(Sys.time() - T0, digits = 3)))
