# ============================================================
# 07_positive_control_bc.R   (pipeline_v2)
#
# POSITIVE DISCOVERY CONTROL. The rod application (stages 04/04c) is a
# NEGATIVE case: the metric panel selects a k=3 split that is not
# reproducible across scDHA embeddings. To show that the pipeline's
# across-embedding stability test ACCEPTS true sub-structure when it
# exists, we apply the same machinery to bipolar cells (BC), which have
# well-established reproducible subtypes (Shekhar et al. 2016).
#
# For BC we: build PCA + scDHA embeddings (5 independent scDHA seeds),
# re-apply Harmony, sweep k, and at each k compute (i) the mean pairwise
# ARI across the 5 seed embeddings (across-embedding reproducibility) and
# (ii) the seed-1 quality metrics (silhouette, ROGUE, balance). The key
# comparison is BC cross-seed ARI (expected HIGH) vs rod cross-seed ARI
# (0.69 at k=3).
#
# Input:  output/rds/02_annotated.rds
# Output: output/tables/07_bc_stability_by_k.csv
#         output/tables/07_bc_seed_sizes.csv
#         output/rds/07_bc_positive_control.rds
# ============================================================

source("config.R"); source("utils.R")
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix)
  library(igraph); library(aricode); library(cluster); library(ROGUE)
  library(tibble); library(dplyr)
})
SEED_LIST <- SCDHA_SEED_LIST          # c(1,2,3,7,42)
K_SWEEP   <- 2:10
T0 <- Sys.time(); log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time()-T0,digits=3), m))

# ---- Load target cell population (CELLTYPE env var; default BC) ----
CELLTYPE <- Sys.getenv("CELLTYPE", "BC")
TAG <- tolower(gsub("[^A-Za-z]", "", CELLTYPE))
log_step(sprintf("Target cell population: %s", CELLTYPE))
if (CELLTYPE == "Rods") {
  # match the stage-04 rod input: rods + cone-contamination filter
  bc <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
  DefaultAssay(bc) <- "RNA"; bc[["RNA"]] <- JoinLayers(bc[["RNA"]]); bc <- NormalizeData(bc, verbose=FALSE)
  cone_g <- intersect(CONE_PANEL, rownames(bc)); rod_g <- intersect(ROD_PANEL, rownames(bc))
  bc <- AddModuleScore(bc, features=list(cone_g), name="cone_score", seed=SEED)
  bc <- AddModuleScore(bc, features=list(rod_g),  name="rod_score",  seed=SEED)
  bc$cone_score <- bc$cone_score1; bc$rod_score <- bc$rod_score1
  cq <- quantile(bc$cone_score, CONE_QUANTILE_CUT, na.rm=TRUE)
  contam <- (bc$cone_score>=cq & bc$cone_score>bc$rod_score) | ((bc$cone_score-bc$rod_score)>CONE_DELTA_CUT)
  bc <- subset(bc, cells=colnames(bc)[!contam])
} else {
  obj <- readRDS(file.path(RDS_DIR, "02_annotated.rds"))
  DefaultAssay(obj) <- "RNA"; obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  ct_col <- if ("cell_type" %in% colnames(obj@meta.data)) "cell_type" else "celltype"
  bc <- subset(obj, cells = colnames(obj)[obj@meta.data[[ct_col]] == CELLTYPE])
  rm(obj)
}
gc(verbose=FALSE)
log_step(sprintf("%s cells: %d", CELLTYPE, ncol(bc)))

bc <- NormalizeData(bc, verbose=FALSE) |> FindVariableFeatures(nfeatures=3000, verbose=FALSE) |>
      ScaleData(verbose=FALSE) |> RunPCA(npcs=N_PCS_CLUSTERING, verbose=FALSE)
pca <- Embeddings(bc, "pca")
har_meta <- data.frame(sample=bc$sample, row.names=colnames(bc))
counts_log2 <- log2(t(as.matrix(GetAssayData(bc, assay="RNA", layer="counts")))+1)

louvain_to_k <- function(latent, target_k) {
  o <- CreateSeuratObject(matrix(0,1,nrow(latent), dimnames=list("d", rownames(latent))))
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

# ---- Train 5 scDHA embeddings + Harmony (PCA-Harmony as seed-independent reference too) ----
suppressPackageStartupMessages({ library(torch); library(scDHA) })
if (!torch_is_installed()) install_torch()
set.seed(SEED)
har_pca <- harmony::RunHarmony(pca, har_meta, vars_use="sample", theta=HARMONY_THETA,
                               lambda=HARMONY_LAMBDA, max_iter=HARMONY_MAX_ITER, verbose=FALSE)
seed_lat <- list()
for (s in SEED_LIST) {
  log_step(sprintf("scDHA seed %d", s)); set.seed(s)
  sc <- scDHA(data=counts_log2, k=NULL, method="scDHA", sparse=FALSE, n=SCDHA_N_GENES,
              ncores=N_CORES, gen_fil=TRUE, do.clus=TRUE, seed=s)
  lat <- as.matrix(sc$latent); rownames(lat)<-rownames(counts_log2)
  set.seed(s)
  seed_lat[[as.character(s)]] <- harmony::RunHarmony(lat, har_meta, vars_use="sample",
      theta=HARMONY_THETA, lambda=HARMONY_LAMBDA, max_iter=HARMONY_MAX_ITER, verbose=FALSE)
}

# ---- ROGUE prefilter (for purity) ----
counts_dense <- as.matrix(GetAssayData(bc, assay="RNA", layer="counts"))
rogue_expr <- ROGUE::matr.filter(counts_dense, min.cells=10, min.genes=10); ROGUE::SE_fun(rogue_expr)
gini <- function(x){x<-sort(x);n<-length(x); if(n<=1||sum(x)==0)return(0); sum((2*seq_len(n)-n-1)*x)/(n*sum(x))}

# ---- Sweep k: cross-seed ARI + seed-1 metrics ----
rows<-list(); size_rows<-list()
ref_ids <- Reduce(intersect, lapply(seed_lat, rownames))
for (k in K_SWEEP) {
  labs <- lapply(seed_lat, function(L) louvain_to_k(L, k))
  # pairwise ARI across seeds
  aris<-c()
  sn<-names(labs)
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
  log_step(sprintf("  k=%d: mean cross-seed ARI=%.3f (sil=%.3f rogue=%.3f bal=%.3f)", k, mean(aris), sil, rg, bal))
}
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, sprintf("07_%s_stability_by_k.csv", TAG)), row.names=FALSE)
write.csv(do.call(rbind,size_rows), file.path(TABLE_DIR, sprintf("07_%s_seed_sizes.csv", TAG)), row.names=FALSE)
saveRDS(list(stability=out, seed_lat=seed_lat, bc=bc), file.path(RDS_DIR, sprintf("07_%s_positive_control.rds", TAG)))
log_step(sprintf("=== %s DISCOVERY-STABILITY ===", CELLTYPE)); print(out)
log_step("DONE")
