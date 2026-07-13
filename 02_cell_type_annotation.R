# ============================================================
# 02_cell_type_annotation.R   (pipeline_v2)
#
# Auto-annotate Seurat clusters using AddModuleScore on the
# CELLTYPE_PANELS marker lists in config.R, then subset to rods.
#
# For each cell:
#   - score = mean module score across the cell's cluster, computed
#     for each cell-type panel
#   - cell-type label = panel with the highest cluster-mean score,
#     provided that score exceeds CELLTYPE_MIN_SCORE; else "Unassigned"
#
# This replaces v1's hard-coded CELLTYPE_MAP which assumed specific
# cluster-IDs from CCA-integrated clustering. Harmony-integrated cluster
# IDs are not stable across runs, so a marker-based assignment is
# more portable.
#
# Input:  output/rds/01_integrated.rds
# Output: output/rds/02_annotated.rds
#         output/rds/02_rods.rds
#         output/tables/02_celltype_assignments.csv
#         output/tables/02_celltype_summary.csv
#         output/figures/fig2_celltypes/
# ============================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})
set.seed(SEED)

# ============================================================
# PHASE 1 — Load + module-score every cell-type panel
# ============================================================
message("=== PHASE 1: Loading object + scoring marker panels ===")

obj <- readRDS(file.path(RDS_DIR, "01_integrated.rds"))
DefaultAssay(obj) <- "SCT"
message("  Cells: ", ncol(obj), " | clusters: ",
        length(levels(Idents(obj))))

panel_score_cols <- character(0)
for (ct in names(CELLTYPE_PANELS)) {
  genes_present <- intersect(CELLTYPE_PANELS[[ct]], rownames(obj))
  if (length(genes_present) < 2) {
    message("  Skipping ", ct, " — only ", length(genes_present),
            " markers present")
    next
  }
  score_col <- paste0("score_", gsub("[^A-Za-z0-9]", "_", ct))
  obj <- AddModuleScore(obj, features = list(genes_present),
                         name = paste0(score_col, "_"),
                         seed = SEED)
  obj@meta.data[[score_col]] <- obj@meta.data[[paste0(score_col, "_1")]]
  obj@meta.data[[paste0(score_col, "_1")]] <- NULL
  panel_score_cols <- c(panel_score_cols, score_col)
  message("  Scored ", ct, " (", length(genes_present), " markers)")
}

# ============================================================
# PHASE 2 — Cluster-mean assignment
# ============================================================
message("\n=== PHASE 2: Cluster-mean assignment ===")

cluster_means <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(panel_score_cols), mean, .names = "{.col}"),
             n = n(), .groups = "drop")

# For each cluster, find the panel with highest score
panel_to_ct <- setNames(names(CELLTYPE_PANELS),
                         paste0("score_",
                                 gsub("[^A-Za-z0-9]", "_",
                                       names(CELLTYPE_PANELS))))

cluster_assignment <- cluster_means %>%
  rowwise() %>%
  mutate(top_score = max(c_across(all_of(panel_score_cols)), na.rm = TRUE),
         top_panel = panel_score_cols[which.max(c_across(
                       all_of(panel_score_cols)))],
         cell_type = ifelse(top_score >= CELLTYPE_MIN_SCORE,
                             panel_to_ct[top_panel], "Unassigned")) %>%
  ungroup() %>%
  select(seurat_clusters, n, top_panel, top_score, cell_type,
          everything())

write.csv(cluster_assignment,
          file.path(TABLE_DIR, "02_celltype_assignments.csv"),
          row.names = FALSE)
message("  Cluster assignments:")
print(cluster_assignment %>% select(seurat_clusters, n, cell_type, top_score),
       n = Inf)

# Apply to cells. Note: bypass Seurat's `obj$col <- vec` setter, which
# requires named vectors and can fail with "No cell overlap" when the
# RHS has no names. Write directly to @meta.data and then refactor.
ct_map <- setNames(cluster_assignment$cell_type,
                    as.character(cluster_assignment$seurat_clusters))
ct_vec <- ct_map[as.character(obj$seurat_clusters)]
names(ct_vec) <- colnames(obj)
obj <- AddMetaData(obj, metadata = ct_vec, col.name = "cell_type")
obj$cell_type <- factor(obj$cell_type,
                         levels = c(names(CELLTYPE_PANELS), "Unassigned"))

# ============================================================
# PHASE 3 — Save annotated object + rod subset
# ============================================================
message("\n=== PHASE 3: Saving annotated object + rod subset ===")

# Cell-type composition table
# ============================================================
# Manual cell-type override
# ------------------------------------------------------------
# If a manual_celltype_override.csv is present in the project root, the
# automated module-score assignments above are *replaced* by the manual
# annotation produced by the v1 pipeline (Louvain res 0.6 -> 34 fine
# clusters -> manual cluster-collapse against canonical markers + DEGs).
# This is the canonical ground truth used for the Phase 1 positive control.
# Cells in the v2 dataset that do not appear in the override CSV (typically
# small numbers due to QC differences) keep their automated label.
# ============================================================
override_path <- "manual_celltype_override.csv"
if (file.exists(override_path)) {
  message("\n=== Applying manual cell-type override from ", override_path, " ===")
  ov <- read.csv(override_path, stringsAsFactors = FALSE)
  ov_lookup <- setNames(ov$cell_type_manual, ov$cell_id)
  cur <- as.character(obj$cell_type)
  names(cur) <- colnames(obj)
  matched <- intersect(names(cur), names(ov_lookup))
  cur[matched] <- ov_lookup[matched]
  # Direct @meta.data assignment to avoid Seurat v5 "no cell overlap"
  obj@meta.data$cell_type_auto   <- as.character(obj$cell_type)
  obj@meta.data$cell_type        <- cur[colnames(obj)]
  message(sprintf("  Overrode %d / %d cells with manual annotation",
                   length(matched), ncol(obj)))
  message(sprintf("  %d cells kept automated label (no manual annotation found)",
                   ncol(obj) - length(matched)))
} else {
  message("\n  No manual override CSV found; keeping automated cell_type")
}

ct_counts <- table(obj$cell_type)
ct_summary <- data.frame(cell_type = names(ct_counts),
                          n         = as.numeric(ct_counts),
                          pct       = round(100 * as.numeric(ct_counts) /
                                              sum(ct_counts), 2))
write.csv(ct_summary, file.path(TABLE_DIR, "02_celltype_summary.csv"),
          row.names = FALSE)
message("\n  Cell-type composition (after any override):")
print(ct_summary)

saveRDS(obj, file.path(RDS_DIR, "02_annotated.rds"))
message("\n  Saved: output/rds/02_annotated.rds")

# Rod subset (uses possibly-overridden cell_type)
rods <- subset(obj, cells = colnames(obj)[obj$cell_type == "Rods"])
saveRDS(rods, file.path(RDS_DIR, "02_rods.rds"))
message("  Saved: output/rds/02_rods.rds  (", ncol(rods), " rod cells)")

# ============================================================
# PHASE 4 — Cell-type figures
# ============================================================
message("\n=== PHASE 4: Figures ===")

# UMAP colored by cell type
p_ct <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
                 cols = CELLTYPE_COLORS, label = TRUE, label.size = 3,
                 repel = TRUE) +
  pub_theme() +
  labs(title = paste0("Harmony UMAP — annotated cell types (n = ",
                       format(ncol(obj), big.mark = ","), ")"))
save_fig(p_ct, "umap_celltypes", FIG_CT, width = 8.5, height = 7)

# Stacked bar of cell-type proportions per sample
prop_df <- obj@meta.data %>%
  count(sample, cell_type) %>%
  group_by(sample) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
all_row <- obj@meta.data %>%
  count(cell_type) %>%
  mutate(sample = "Overall", pct = 100 * n / sum(n))
prop_full <- bind_rows(prop_df, all_row)

p_bar <- ggplot(prop_full, aes(x = sample, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = CELLTYPE_COLORS, drop = FALSE) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Cell-type proportions by sample",
       x = "Sample", y = "Percent of cells")
save_fig(p_bar, "celltype_proportions", FIG_CT, width = 9, height = 6)

# Marker dot plot organized by cell type
DefaultAssay(obj) <- "SCT"
ordered_markers <- unique(unlist(CELLTYPE_PANELS))
ordered_markers <- intersect(ordered_markers, rownames(obj))
p_dot <- DotPlot(obj, features = ordered_markers, group.by = "cell_type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  pub_theme() +
  labs(title = "Canonical markers by annotated cell type",
       x = "Gene", y = "Cell type")
save_fig(p_dot, "dotplot_markers_by_celltype", FIG_CT, width = 16, height = 6)

message("\n=== 02_cell_type_annotation.R complete ===")
