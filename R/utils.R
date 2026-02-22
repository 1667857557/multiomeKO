# ---- utilities ----

.check_cols <- function(df, cols) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse=", "))
}

# sparse group design matrix: cells x groups
.group_design <- function(groups) {
  groups <- as.factor(groups)
  Matrix::sparse.model.matrix(~0 + groups)
}

# aggregate features x cells sparse/dense matrix into features x groups
.aggregate_by_group <- function(mat, groups, fun = c("sum", "mean")) {
  fun <- match.arg(fun)
  G <- .group_design(groups)   # cells x groups
  out <- mat %*% G             # features x groups
  if (fun == "mean") {
    n <- Matrix::colSums(G)
    n[n == 0] <- 1
    out <- out %*% Matrix::Diagonal(x = 1 / n)
  }
  colnames(out) <- colnames(G)
  out
}

# log1p(CPM) for sparse counts (features x cells)
.log1p_cpm <- function(counts, scale = 1e4) {
  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- 1
  norm <- counts %*% Matrix::Diagonal(x = scale / lib)
  log1p(norm)
}

# safe exp for sparse/dense matrices (returns dense for numeric stability)
.exp_mat <- function(x) {
  # Matrix::exp exists but can be slow; base exp works on dense.
  if (inherits(x, "dgCMatrix") || inherits(x, "dgTMatrix")) {
    x <- as.matrix(x)
  }
  exp(x)
}

# safe align vector to matrix rownames
.align_vec <- function(v, rn) {
  if (is.null(names(v))) return(rep(v, length(rn)))
  v[rn]
}

# keep only finite values
.finite0 <- function(x, fill = 0) {
  x[!is.finite(x)] <- fill
  x
}
