# ---- metacell construction ----

#' Create metacell ids by fixed-size pseudobulk within groups
#'
#' Biological rationale:
#' - scRNA/scATAC are sparse; per-cell counts are noisy.
#' - Pseudobulk/metacell reduces measurement noise while preserving state specificity.
#' - Fixed-size buckets reduce confounding from unequal cell numbers and depth.
#'
#' @param obj Seurat object with paired RNA+ATAC.
#' @param group.by Character vector of metadata columns to stratify (e.g. celltype, SampleID).
#' @param n_cells Target number of cells per metacell within each stratum.
#' @param min_cells Minimum cells required to keep a stratum; smaller strata are kept as one metacell.
#' @param seed Random seed.
#' @return A factor vector of metacell ids aligned to cells (colnames(obj)).
make_metacell_ids <- function(obj, group.by, n_cells = 50, min_cells = 25, seed = 1) {
  md <- obj@meta.data
  if (!all(group.by %in% colnames(md))) {
    stop("group.by columns not found in obj@meta.data: ", paste(setdiff(group.by, colnames(md)), collapse=", "))
  }
  set.seed(seed)
  base_group <- do.call(paste, c(md[, group.by, drop=FALSE], sep = "__"))
  base_group <- as.character(base_group)
  cells <- colnames(obj)

  metacell <- rep(NA_character_, length(cells))
  names(metacell) <- cells

  for (g in unique(base_group)) {
    idx <- which(base_group == g)
    if (length(idx) == 0) next
    # if small, keep one metacell
    if (length(idx) < min_cells) {
      metacell[idx] <- paste0(g, "__mc1")
      next
    }
    # partition into fixed-size buckets
    idx <- sample(idx)
    k <- ceiling(length(idx) / n_cells)
    split_ids <- rep(seq_len(k), each = n_cells)[seq_along(idx)]
    metacell[idx] <- paste0(g, "__mc", split_ids)
  }

  factor(metacell)
}
