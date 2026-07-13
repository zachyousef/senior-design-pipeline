# ============================================================
# utils.R — Shared utility functions for the pipeline
# ============================================================

# ---- Publication-quality ggplot2 theme ----
pub_theme <- function(base_size = PUB_BASE_SIZE) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text          = element_text(color = "black"),
      plot.title    = element_text(size = base_size + 2, face = "bold", hjust = 0,
                                   margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, hjust = 0, color = "grey30"),
      axis.text     = element_text(size = base_size - 1, color = "black"),
      axis.title    = element_text(size = base_size, color = "black"),
      axis.line     = element_line(color = "black", linewidth = 0.4),
      axis.ticks    = element_line(color = "black", linewidth = 0.3),
      legend.text   = element_text(size = base_size - 1),
      legend.title  = element_text(size = base_size, face = "bold"),
      legend.key.size = unit(0.4, "cm"),
      strip.text    = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      panel.grid    = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin   = margin(8, 8, 8, 8, "mm")
    )
}

# ---- Save figure as BOTH PDF + PNG (600 DPI) ----
save_fig <- function(plot, filename, dir, width = 7, height = 5) {
  pdf_path <- file.path(dir, paste0(filename, ".pdf"))
  png_path <- file.path(dir, paste0(filename, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height,
         device = "pdf")
  ggsave(png_path, plot = plot, width = width, height = height,
         dpi = FIG_DPI, device = "png", bg = "white")
  message("  Saved: ", filename, " (.pdf + .png)")
}

# ---- Save figure as TIFF (high-res raster) ----
save_tiff <- function(plot, filename, dir, width = 7, height = 5) {
  filepath <- file.path(dir, paste0(filename, ".tiff"))
  ggsave(filepath, plot = plot, width = width, height = height,
         dpi = FIG_DPI, device = "tiff", compression = "lzw")
  message("  Saved: ", filename, ".tiff")
}

# ---- Read 10X data robustly ----
read10x_robust <- function(path) {
  tryCatch(
    Seurat::Read10X(data.dir = path, gene.column = 2),
    error = function(e1) {
      tryCatch(
        Seurat::Read10X(data.dir = path, gene.column = 1),
        error = function(e2) stop("Read10X failed for ", path, "\n", e1$message)
      )
    }
  )
}

# ---- Compute mean silhouette width ----
compute_silhouette <- function(embedding, labels, max_cells = MAX_SIL_CELLS,
                                seed = SEED) {
  k <- length(unique(labels))
  if (k < 2) return(NA_real_)

  n <- nrow(embedding)
  if (n > max_cells) {
    set.seed(seed)
    idx <- sample(n, max_cells)
    embedding <- embedding[idx, , drop = FALSE]
    labels    <- labels[idx]
  }

  d   <- dist(embedding)
  sil <- cluster::silhouette(as.integer(factor(labels)), d)
  mean(sil[, "sil_width"])
}

# ---- Compute mean ROGUE score ----
compute_rogue_score <- function(counts_mat, labels) {
  requireNamespace("ROGUE", quietly = TRUE)
  k <- length(unique(labels))
  if (k < 1) return(NA_real_)

  expr_filt <- ROGUE::matr.filter(counts_mat, min.cells = 10, min.genes = 10)
  ROGUE::SE_fun(expr_filt)

  rogue_res <- tryCatch(
    ROGUE::rogue(expr_filt, labels = labels,
                 samples = rep(1, length(labels)), platform = "UMI"),
    error = function(e) { message("  ROGUE error: ", e$message); return(NA) }
  )
  if (all(is.na(rogue_res))) return(NA_real_)
  mean(as.numeric(rogue_res), na.rm = TRUE)
}

# ---- Build SNN graph via Seurat::FindNeighbors (denser, Leiden-friendly) ----
# Uses Seurat's SNN construction (which the Louvain methods already use)
# and converts the resulting sparse SNN matrix into an undirected
# weighted igraph. Far more robust for Leiden than the raw kNN+Jaccard
# build_snn_graph below.
build_seurat_snn_igraph <- function(embedding, cell_ids, counts_for_dummy,
                                     reduction_name, key_prefix, dims_use) {
  requireNamespace("Seurat", quietly = TRUE)
  requireNamespace("igraph", quietly = TRUE)
  requireNamespace("Matrix", quietly = TRUE)

  tmp_obj <- Seurat::CreateSeuratObject(counts = counts_for_dummy[, cell_ids])
  tmp_obj[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = embedding, key = key_prefix, assay = "RNA")
  tmp_obj <- Seurat::FindNeighbors(tmp_obj, reduction = reduction_name,
                                    dims = dims_use, verbose = FALSE)
  snn_mat <- tmp_obj@graphs[["RNA_snn"]]
  rm(tmp_obj); gc(verbose = FALSE)

  snn_mat <- methods::as(snn_mat, "TsparseMatrix")
  keep    <- snn_mat@i < snn_mat@j
  edges   <- cbind(snn_mat@i[keep] + 1L, snn_mat@j[keep] + 1L)
  weights <- snn_mat@x[keep]

  g <- igraph::make_graph(t(edges), n = ncol(snn_mat), directed = FALSE)
  igraph::E(g)$weight <- weights
  igraph::V(g)$name   <- rownames(snn_mat)
  g
}

# ---- Build SNN graph on an embedding (raw kNN+Jaccard) ----
build_snn_graph <- function(embedding, cell_ids, k = KNN_K, prune = SNN_PRUNE) {
  requireNamespace("FNN",    quietly = TRUE)
  requireNamespace("igraph", quietly = TRUE)

  n <- nrow(embedding)
  message("  Building SNN graph: ", n, " cells, kNN k=", k)

  knn    <- FNN::get.knn(embedding, k = k)
  nn_idx <- knn$nn.index
  storage.mode(nn_idx) <- "integer"

  from_idx <- rep(seq_len(n), each = k)
  to_idx   <- as.vector(nn_idx)
  ok       <- !is.na(to_idx) & to_idx > 0L & from_idx != to_idx

  a <- pmin(from_idx[ok], to_idx[ok])
  b <- pmax(from_idx[ok], to_idx[ok])
  cand <- unique(data.frame(a = a, b = b))

  nbrs <- vector("list", n)
  for (i in seq_len(n)) {
    v <- nn_idx[i, ]
    nbrs[[i]] <- sort(unique(v[!is.na(v) & v > 0L & v != i]))
  }

  jacc <- vapply(seq_len(nrow(cand)), function(idx) {
    ni <- nbrs[[cand$a[idx]]]; nj <- nbrs[[cand$b[idx]]]
    inter <- length(intersect(ni, nj))
    denom <- length(ni) + length(nj) - inter
    if (denom > 0) inter / denom else 0
  }, numeric(1))

  keep_e <- jacc > prune
  if (sum(keep_e) == 0) stop("All SNN edges pruned. Lower SNN_PRUNE in config.R.")

  edge_df <- data.frame(from = cell_ids[cand$a[keep_e]],
                         to = cell_ids[cand$b[keep_e]],
                         weight = jacc[keep_e], stringsAsFactors = FALSE)
  vtx_df  <- data.frame(name = cell_ids, stringsAsFactors = FALSE)

  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = vtx_df)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                        edge.attr.comb = list(weight = "max"))
  message("  SNN graph: ", igraph::vcount(g), " vertices, ",
          igraph::ecount(g), " edges")
  g
}

# ---- Tune resolution for target k clusters ----
find_resolution_for_k <- function(g, target_k, algorithm = "louvain",
                                  res_min = 1e-6, res_max = 50,
                                  max_iter = 20, seed = SEED) {
  requireNamespace("igraph", quietly = TRUE)

  if (target_k == 1) {
    labs <- rep(1L, igraph::vcount(g))
    names(labs) <- igraph::V(g)$name
    return(list(labels = labs, resolution = 0, actual_k = 1, exact = TRUE))
  }

  cluster_at_res <- function(res) {
    set.seed(seed)
    w <- igraph::E(g)$weight
    if (algorithm == "louvain") {
      comm <- igraph::cluster_louvain(g, weights = w, resolution = res)
    } else {
      # CPM (Constant Potts Model) is Leiden's default and avoids the
      # modularity resolution limit (Fortunato & Barthélemy 2007). On
      # graphs with one dominant community, modularity-Leiden peels off
      # tiny outliers instead of splitting the bulk; CPM gives balanced
      # partitions across the resolution range.
      comm <- igraph::cluster_leiden(g, weights = w,
                                     resolution_parameter = res,
                                     objective_function = "CPM",
                                     n_iterations = 15L)
    }
    memb <- as.integer(factor(igraph::membership(comm)))
    names(memb) <- igraph::V(g)$name
    list(labels = memb, n_clusters = length(unique(memb)))
  }

  best      <- NULL
  best_diff <- Inf
  lo <- res_min; hi <- res_max

  for (iter in seq_len(max_iter)) {
    mid <- sqrt(lo * hi)
    res <- cluster_at_res(mid)
    diff <- abs(res$n_clusters - target_k)
    if (diff < best_diff) {
      best <- list(labels = res$labels, resolution = mid,
                   actual_k = res$n_clusters, exact = (diff == 0))
      best_diff <- diff
    }
    if (diff == 0) break
    if (res$n_clusters < target_k) lo <- mid else hi <- mid
  }

  if (best_diff > 0) {
    grid <- exp(seq(log(res_min), log(res_max), length.out = 60))
    for (r in grid) {
      res <- cluster_at_res(r)
      diff <- abs(res$n_clusters - target_k)
      if (diff < best_diff) {
        best <- list(labels = res$labels, resolution = r,
                     actual_k = res$n_clusters, exact = (diff == 0))
        best_diff <- diff
      }
      if (diff == 0) break
    }
  }
  best
}
