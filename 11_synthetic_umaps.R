# ============================================================
# 11_synthetic_umaps.R   (pipeline_v2)
#
# Visual of the synthetic datasets used for the ground-truth validation
# (stage 09): re-simulates the representative scenarios with the same
# seeds, computes a PCA-based UMAP, and saves the 2D coordinates colored
# by the TRUE simulated group. Feeds the top row of Figure 6.
#
# Output: output/tables/11_synthetic_umap_coords.csv
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({
  library(splatter); library(SingleCellExperiment); library(Seurat); library(uwot)
})
NG <- 5000; NGENES <- 2000
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0,digits=3), m))

scenarios <- list(
  list(id="null",     type="null",       group.prob=1,              de.prob=0,    seed=101, lab="Null (no structure)"),
  list(id="weak",     type="structured", group.prob=c(.34,.33,.33), de.prob=0.03, seed=505, lab="Weak (de=0.03)"),
  list(id="strong",   type="structured", group.prob=c(.34,.33,.33), de.prob=0.30, seed=303, lab="Strong (de=0.30)"),
  list(id="rare",     type="structured", group.prob=c(.80,.15,.05), de.prob=0.15, seed=606, lab="Rare 5% subpopulation")
)

out <- list()
for (sc in scenarios) {
  log_step(sprintf("simulate %s", sc$id))
  params <- newSplatParams(nGenes=NGENES, batchCells=NG, seed=sc$seed)
  if (sc$type=="null") { sim <- splatSimulate(params, verbose=FALSE); truth <- rep("G1", ncol(sim)) }
  else { sim <- splatSimulateGroups(params, group.prob=sc$group.prob, de.prob=sc$de.prob, verbose=FALSE)
         truth <- as.character(colData(sim)$Group) }
  cts <- as.matrix(counts(sim)); colnames(cts) <- paste0("c", seq_len(ncol(cts)))
  so <- CreateSeuratObject(cts) |> NormalizeData(verbose=FALSE) |>
        FindVariableFeatures(nfeatures=2000, verbose=FALSE) |> ScaleData(verbose=FALSE) |>
        RunPCA(npcs=30, verbose=FALSE)
  set.seed(SEED)
  um <- uwot::umap(Embeddings(so,"pca")[,1:30], n_neighbors=30, min_dist=0.3, verbose=FALSE)
  # relabel groups by size (Group1 = largest) for consistent coloring
  tb <- sort(table(truth), decreasing=TRUE); remap <- setNames(seq_along(tb), names(tb))
  grp <- paste0("group ", remap[truth])
  out[[sc$id]] <- data.frame(scenario=sc$lab, order=which(sapply(scenarios,function(z)z$id)==sc$id),
                             UMAP1=um[,1], UMAP2=um[,2], group=grp)
}
res <- do.call(rbind, out)
write.csv(res, file.path(TABLE_DIR, "11_synthetic_umap_coords.csv"), row.names=FALSE)
log_step(sprintf("Saved UMAP coords: %d cells across %d scenarios", nrow(res), length(scenarios)))
