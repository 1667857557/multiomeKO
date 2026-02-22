# ---- extracting T/A/X from Seurat+Signac ----

#' Extract Regulator (T), ATAC (A), RNA (X) matrices on shared metacells
#'
#' Key design choices (objective):
#' 1) Same metacell ids are used for RNA and ATAC to preserve paired structure.
#' 2) Uses counts for Poisson/NB-style modeling and provides offsets (library sizes).
#' 3) Provides log1p(A_counts) as a stable ATAC feature for downstream RNA modeling.
#'
#' @param obj Seurat object.
#' @param metacell_ids Factor/character vector length ncol(obj), or NULL to create with make_metacell_ids().
#' @param group.by Used only if metacell_ids is NULL.
#' @param rna_assay RNA assay name.
#' @param atac_assay ATAC assay name.
#' @param genes_use Character vector of genes to keep in X (optional).
#' @param peaks_use Character vector of peaks to keep in A (optional).
#' @param tf_mode "RNA", "chromvar", or "both".
#' @param tf_genes TF genes (required if tf_mode includes RNA).
#' @param covariates Optional metadata columns to aggregate as unpenalized covariates.
#' @return list(T, A_counts, A_log1p, X_counts, offsets, covariates, metacell_ids)
extract_TAX_metacell <- function(
  obj,
  metacell_ids = NULL,
  group.by = NULL,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  genes_use = NULL,
  peaks_use = NULL,
  tf_mode = c("RNA", "chromvar", "both"),
  tf_genes = NULL,
  covariates = NULL
) {
  tf_mode <- match.arg(tf_mode)

  if (is.null(metacell_ids)) {
    if (is.null(group.by)) stop("Provide metacell_ids or group.by")
    metacell_ids <- make_metacell_ids(obj, group.by = group.by)
  } else {
    if (length(metacell_ids) != ncol(obj)) stop("metacell_ids length must equal ncol(obj)")
    metacell_ids <- as.factor(metacell_ids)
  }

  # RNA counts
  Xc <- Seurat::GetAssayData(obj, assay = rna_assay, slot = "counts")
  if (!is.null(genes_use)) {
    genes_use <- intersect(genes_use, rownames(Xc))
    Xc <- Xc[genes_use, , drop=FALSE]
  }

  # ATAC counts
  Ac <- Seurat::GetAssayData(obj, assay = atac_assay, slot = "counts")
  if (!is.null(peaks_use)) {
    peaks_use <- intersect(peaks_use, rownames(Ac))
    Ac <- Ac[peaks_use, , drop=FALSE]
  }

  # Aggregate to metacells (SUM for count models)
  X_mc <- .aggregate_by_group(Xc, metacell_ids, fun = "sum")
  A_mc <- .aggregate_by_group(Ac, metacell_ids, fun = "sum")

  # Offsets
  off_rna <- log(Matrix::colSums(X_mc) + 1)
  off_atac <- log(Matrix::colSums(A_mc) + 1)

  # ATAC feature for stage2 (stable across prediction): log1p(count)
  A_log1p <- log1p(A_mc)

  # Regulators T
  T_list <- list()
  if (tf_mode %in% c("RNA", "both")) {
    if (is.null(tf_genes)) stop("tf_genes is required when tf_mode includes RNA")
    tf_genes <- intersect(tf_genes, rownames(Xc))
    if (length(tf_genes) == 0) stop("No tf_genes found in RNA counts")
    # Use log1p CPM-like scale for TF RNA as covariates (not as response)
    X_tf_log <- .log1p_cpm(Xc[tf_genes, , drop=FALSE])
    T_rna <- .aggregate_by_group(X_tf_log, metacell_ids, fun = "mean")
    rownames(T_rna) <- paste0("RNA:", rownames(T_rna))
    T_list[["RNA"]] <- T_rna
  }

  if (tf_mode %in% c("chromvar", "both")) {
    if (!"chromvar" %in% names(obj@assays)) {
      stop("chromvar assay not found. RunChromVAR() or use tf_mode='RNA'.")
    }
    T_chr <- Seurat::GetAssayData(obj, assay = "chromvar", slot = "data")
    T_chr <- .aggregate_by_group(T_chr, metacell_ids, fun = "mean")
    rownames(T_chr) <- paste0("MOTIF:", rownames(T_chr))
    T_list[["chromvar"]] <- T_chr
  }

  T <- do.call(rbind, T_list)

  # Optional covariates aggregated to metacell (unpenalized)
  C <- NULL
  if (!is.null(covariates)) {
    md <- obj@meta.data
    covariates <- intersect(covariates, colnames(md))
    if (length(covariates) > 0) {
      # numeric only
      cov_df <- md[, covariates, drop=FALSE]
      keep_num <- vapply(cov_df, is.numeric, logical(1))
      cov_df <- cov_df[, keep_num, drop=FALSE]
      if (ncol(cov_df) > 0) {
        # aggregate mean per metacell
        cov_mat <- t(as.matrix(cov_df))
        C <- .aggregate_by_group(cov_mat, metacell_ids, fun = "mean")
        rownames(C) <- paste0("COV:", rownames(C))
      }
    }
  }

  list(
    T = T,
    A_counts = A_mc,
    A_log1p = A_log1p,
    X_counts = X_mc,
    offsets = list(rna = off_rna, atac = off_atac),
    covariates = C,
    metacell_ids = metacell_ids
  )
}
