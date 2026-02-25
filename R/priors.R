# ---- motif and peak2gene priors ----

#' Extract peak2gene links from Signac Links()
#' @param obj Seurat object
#' @param atac_assay ATAC assay name
#' @return data.frame with columns peak,gene and optional score
get_peak2gene_links <- function(obj, atac_assay = "ATAC") {
  lk <- Signac::Links(obj[[atac_assay]])
  if (nrow(lk) == 0) stop("No Links found in atac assay. Run LinkPeaks() first.")
  df <- as.data.frame(lk)

  # robust column detection
  peak_col <- NULL
  if ("peak" %in% names(df)) peak_col <- "peak"
  if (is.null(peak_col) && "query_region" %in% names(df)) peak_col <- "query_region"
  if (is.null(peak_col)) stop("Cannot find peak column in Links data.frame; expected one of: peak, query_region")

  gene_col <- NULL
  if ("gene" %in% names(df)) gene_col <- "gene"
  if (is.null(gene_col) && "gene_name" %in% names(df)) gene_col <- "gene_name"
  if (is.null(gene_col)) stop("Cannot find gene column in Links data.frame")

  out <- df[, c(peak_col, gene_col), drop=FALSE]
  colnames(out) <- c("peak", "gene")

  # optional score
  score_col <- NULL
  for (cand in c("score", "zscore", "cor", "correlation")) {
    if (cand %in% names(df)) { score_col <- cand; break }
  }
  if (!is.null(score_col)) out$score <- df[[score_col]]

  out
}

#' Extract peak x motif map and motif-to-TF symbol mapping
#'
#' This uses Signac motif data matrix (usually binary) and parses TF symbol
#' from motif names (e.g., MA0139.1_FOSL1 -> FOSL1).
#'
#' @param obj Seurat object
#' @param atac_assay ATAC assay name
#' @return list(motif_mat, motif_names, tf_symbols)
get_peak_motif_map <- function(obj, atac_assay = "ATAC") {
  mm <- Signac::GetMotifData(obj, assay = atac_assay, slot = "data")
  peaks <- rownames(obj[[atac_assay]])

  # ensure peaks are rows
  if (!all(peaks %in% rownames(mm)) && all(peaks %in% colnames(mm))) {
    mm <- t(mm)
  }
  motif_names <- colnames(mm)
  tf_symbols <- sub(".*_", "", motif_names)
  list(motif_mat = mm, motif_names = motif_names, tf_symbols = tf_symbols)
}


#' Build stage1 priors with pluggable modes
#'
#' @param obj Seurat object
#' @param atac_assay ATAC assay name
#' @param stage1_mode one of "motif" (TF only), "chip" (external peak sets), "hybrid" (both)
#' @param chip_peak_map optional named list for non-TF regulators (regulator -> peak vector)
#' @return list with mode-specific prior maps and diagnostics
build_stage1_priors <- function(
  obj,
  atac_assay = "ATAC",
  stage1_mode = c("motif", "chip", "hybrid"),
  chip_peak_map = NULL
) {
  stage1_mode <- match.arg(stage1_mode)
  out <- list(stage1_mode = stage1_mode)

  if (stage1_mode %in% c("motif", "hybrid")) {
    out$motif_map <- get_peak_motif_map(obj, atac_assay = atac_assay)
  }

  if (stage1_mode %in% c("chip", "hybrid")) {
    if (is.null(chip_peak_map) || !is.list(chip_peak_map) || is.null(names(chip_peak_map))) {
      stop("chip_peak_map must be a named list(regulator -> peak character vector) for chip/hybrid mode")
    }
    peaks <- rownames(obj[[atac_assay]])
    chip_peak_map <- lapply(chip_peak_map, function(pk) intersect(unique(as.character(pk)), peaks))
    out$chip_peak_map <- chip_peak_map
    out$chip_summary <- data.frame(
      regulator = names(chip_peak_map),
      n_peaks = vapply(chip_peak_map, length, integer(1)),
      stringsAsFactors = FALSE
    )
  }

  out
}
