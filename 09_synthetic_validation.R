# ============================================================
# 09_synthetic_validation.R   (pipeline_v2)
#
# Ground-truth validation of the across-embedding stability test on
# synthetic scRNA-seq data (Splatter). We simulate:
#   - NULL datasets: one homogeneous group (no sub-structure).
#   - STRUCTURED datasets: real groups at varying separation (de.prob),
#     including a rare-subpopulation case (unbalanced group sizes).
# For each dataset we train five independent scDHA embeddings and, at
# each k, compute (i) across-embedding reproducibility (mean pairwise
# ARI between the five partitions) and (ii) ARI against the TRUE labels.
# The test should ACCEPT structured data (reproducibility high at the
# true k and recovering the truth) and REJECT null data (reproducibility
# low at all k) -> sensitivity and specificity of the criterion.
#
# Output: output/tables/09_synthetic_validation.csv
#         output/rds/09_synthetic_validation.rds
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({
  library(splatter); library(SingleCellExperiment); library(Seurat)
  library(aricode); library(cluster)
})
NG <- 5000; NGENES <- 2000; SEEDS <- SCDHA_SEED_LIST; KS <- 2:6
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0,digits=3), m))

# ---- Simulation scenarios ----
scenarios <- list(
  list(id="null_A",     type="null",       true_k=1, group.prob=1,               de.prob=0,    seed=101),
  list(id="null_B",     type="null",       true_k=1, group.prob=1,               de.prob=0,    seed=202),
  list(id="strong_3",   type="structured", true_k=3, group.prob=c(.34,.33,.33),  de.prob=0.30, seed=303),
  list(id="medium_3",   type="structured", true_k=3, group.prob=c(.34,.33,.33),  de.prob=0.10, seed=404),
  list(id="weak_3",     type="structured", true_k=3, group.prob=c(.34,.33,.33),  de.prob=0.03, seed=505),
  list(id="rare_3",     type="structured", true_k=3, group.prob=c(.80,.15,.05),  de.prob=0.15, seed=606)
)

louvain_to_k <- function(latent, target_k) {
  o <- CreateSeuratObject(matrix(0,1,nrow(latent), dimnames=list("d", rownames(latent))))
  o[["lat"]] <- CreateDimReducObject(embeddings=latent, key="L_", assay="RNA")
  o <- FindNeighbors(o, reduction="lat", dims=1:ncol(latent), verbose=FALSE)
  bl<-NULL; bd<-Inf
  for (r in exp(seq(log(1e-4), log(10), length.out=50))) {
    set.seed(SEED); o<-FindClusters(o, resolution=r, algorithm=1, verbose=FALSE)
    lab<-as.integer(Idents(o)); nk<-length(unique(lab))
    if (abs(nk-target_k)<bd){bd<-abs(nk-target_k); bl<-lab}; if(bd==0)break
  }
  names(bl)<-rownames(latent); bl
}

suppressPackageStartupMessages({ library(torch); library(scDHA) })
if (!torch_is_installed()) install_torch()

rows <- list()
for (sc in scenarios) {
  log_step(sprintf("=== simulate %s (%s, true_k=%d, de=%.2f) ===", sc$id, sc$type, sc$true_k, sc$de.prob))
  params <- newSplatParams(nGenes=NGENES, batchCells=NG, seed=sc$seed)
  if (sc$type=="null") {
    sim <- splatSimulate(params, verbose=FALSE); truth <- rep("G1", ncol(sim))
  } else {
    sim <- splatSimulateGroups(params, group.prob=sc$group.prob, de.prob=sc$de.prob, verbose=FALSE)
    truth <- as.character(colData(sim)$Group)
  }
  cts <- as.matrix(counts(sim)); colnames(cts) <- paste0("c", seq_len(ncol(cts)))
  names(truth) <- colnames(cts)
  data_scdha <- log2(t(cts) + 1)

  # 5 scDHA embeddings
  lat <- list()
  for (s in SEEDS) { set.seed(s)
    r <- scDHA(data=data_scdha, k=NULL, method="scDHA", sparse=FALSE, n=min(NGENES,3000),
               ncores=N_CORES, gen_fil=TRUE, do.clus=TRUE, seed=s)
    m <- as.matrix(r$latent); rownames(m) <- rownames(data_scdha); lat[[as.character(s)]] <- m
  }
  ids <- Reduce(intersect, lapply(lat, rownames))
  for (k in KS) {
    labs <- lapply(lat, function(L) louvain_to_k(L, k))
    aris <- c(); sn <- names(labs)
    for (i in 1:(length(sn)-1)) for (j in (i+1):length(sn)) aris <- c(aris, aricode::ARI(labs[[i]][ids], labs[[j]][ids]))
    ari_truth <- aricode::ARI(labs[["1"]][ids], truth[ids])
    sil <- tryCatch(compute_silhouette(lat[["1"]][ids,,drop=FALSE], labs[["1"]][ids],
                    max_cells=MAX_SIL_CELLS, seed=SEED), error=function(e) NA_real_)
    rows[[length(rows)+1]] <- data.frame(scenario=sc$id, type=sc$type, true_k=sc$true_k, de_prob=sc$de.prob,
        k=k, mean_cross_seed_ari=round(mean(aris),4), ari_vs_truth=round(ari_truth,4), silhouette=round(sil,4))
    log_step(sprintf("  k=%d: cross-seed ARI=%.3f | ARI vs truth=%.3f | sil=%.3f", k, mean(aris), ari_truth, sil))
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "09_synthetic_validation.csv"), row.names=FALSE)
saveRDS(out, file.path(RDS_DIR, "09_synthetic_validation.rds"))

# ---- Sensitivity / specificity of the accept rule ----
# Accept if reproducibility at the true k >= 0.7 (structured) OR max across k >= 0.7 (null -> should stay low)
peak <- aggregate(mean_cross_seed_ari ~ scenario + type + true_k, data=out, FUN=max)
at_truek <- do.call(rbind, lapply(split(out, out$scenario), function(df) df[df$k==df$true_k[1] | (df$true_k[1]==1 & df$k==2),][1,]))
log_step("=== SUMMARY: peak cross-seed ARI per scenario ==="); print(peak)
log_step("DONE")
