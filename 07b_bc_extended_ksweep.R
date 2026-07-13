# ============================================================
# 07b_bc_extended_ksweep.R   (pipeline_v2)
#
# Extends the bipolar positive-control k-sweep from k=2-10 to k=2-18 to
# characterize how across-embedding reproducibility behaves relative to the
# ~15 transcriptomic bipolar types reported by Shekhar et al. (2016): does it
# keep rising, plateau, or decline (over-split)?
#
# Reuses the FIVE cached scDHA+Harmony embeddings already trained in
# stage 07 (output/rds/07_bc_positive_control.rds) — NO scDHA retraining.
# Metric computations mirror stage 07 exactly for consistency.
#
# Input:  output/rds/07_bc_positive_control.rds  (seed_lat + bc)
# Output: output/tables/07_bc_stability_by_k_extended.csv
#         output/tables/07_bc_seed_sizes_extended.csv
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({ library(Seurat); library(aricode); library(ROGUE); library(tibble); library(dplyr) })

K_EXT <- as.integer(strsplit(Sys.getenv("K_EXT", "2:18"), ":")[[1]])
K_EXT <- seq(K_EXT[1], K_EXT[2])
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0, digits=3), m))

pc <- readRDS(file.path(RDS_DIR, "07_bc_positive_control.rds"))
seed_lat <- pc$seed_lat
bc <- pc$bc
log_step(sprintf("Loaded %d cached embeddings (seeds %s); %d cells",
                 length(seed_lat), paste(names(seed_lat), collapse=","), ncol(bc)))

# ---- Louvain-to-target-k on a fixed latent (identical to stage 07) ----
louvain_to_k <- function(latent, target_k) {
  o <- CreateSeuratObject(matrix(0, 1, nrow(latent), dimnames=list("d", rownames(latent))))
  o[["lat"]] <- CreateDimReducObject(embeddings=latent, key="L_", assay="RNA")
  o <- FindNeighbors(o, reduction="lat", dims=1:ncol(latent), verbose=FALSE)
  bl<-NULL; bd<-Inf
  for (r in exp(seq(log(1e-4), log(10), length.out=60))) {
    set.seed(SEED); o<-FindClusters(o, resolution=r, algorithm=1, verbose=FALSE)
    lab<-as.integer(Idents(o)); nk<-length(unique(lab))
    if (abs(nk-target_k)<bd){bd<-abs(nk-target_k); bl<-lab}; if(bd==0)break
  }
  names(bl)<-rownames(latent); bl
}

# ---- ROGUE prefilter (identical to stage 07) ----
counts_dense <- as.matrix(GetAssayData(bc, assay="RNA", layer="counts"))
rogue_expr <- ROGUE::matr.filter(counts_dense, min.cells=10, min.genes=10); ROGUE::SE_fun(rogue_expr)
gini <- function(x){x<-sort(x);n<-length(x); if(n<=1||sum(x)==0)return(0); sum((2*seq_len(n)-n-1)*x)/(n*sum(x))}

# ---- Sweep k: cross-seed ARI + seed-1 metrics ----
rows<-list(); size_rows<-list()
ref_ids <- Reduce(intersect, lapply(seed_lat, rownames))
for (k in K_EXT) {
  labs <- lapply(seed_lat, function(L) louvain_to_k(L, k))
  aris<-c(); sn<-names(labs)
  for (i in 1:(length(sn)-1)) for (j in (i+1):length(sn))
    aris<-c(aris, aricode::ARI(labs[[i]][ref_ids], labs[[j]][ref_ids]))
  l1 <- labs[["1"]]
  sil <- tryCatch(compute_silhouette(seed_lat[["1"]], l1, max_cells=MAX_SIL_CELLS, seed=SEED), error=function(e)NA)
  rg  <- tryCatch({ r<-ROGUE::rogue(rogue_expr, labels=as.character(l1[colnames(rogue_expr)]),
                    samples=rep(1,ncol(rogue_expr)), platform="UMI"); mean(unlist(r),na.rm=TRUE)}, error=function(e)NA)
  bal <- 1 - gini(as.numeric(table(l1)))
  rows[[length(rows)+1]] <- data.frame(k=k, actual_k=length(unique(l1)),
      mean_cross_seed_ari=round(mean(aris),4), min_cross_seed_ari=round(min(aris),4),
      silhouette=round(sil,4), rogue=round(rg,4), balance=round(bal,4))
  size_rows[[length(size_rows)+1]] <- data.frame(k=k, seed1_sizes=paste(sort(as.integer(table(l1)),decreasing=TRUE),collapse="/"))
  log_step(sprintf("  k=%d: mean cross-seed ARI=%.3f (min=%.3f sil=%.3f rogue=%.3f bal=%.3f)",
                   k, mean(aris), min(aris), sil, rg, bal))
}
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, "07_bc_stability_by_k_extended.csv"), row.names=FALSE)
write.csv(do.call(rbind,size_rows), file.path(TABLE_DIR, "07_bc_seed_sizes_extended.csv"), row.names=FALSE)
log_step(sprintf("Done. Peak mean cross-seed ARI = %.3f at k=%d (range k=%d-%d).",
                 max(out$mean_cross_seed_ari), out$k[which.max(out$mean_cross_seed_ari)], min(K_EXT), max(K_EXT)))
