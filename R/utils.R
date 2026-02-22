# ---- utilities ----

.check_cols <- function(df, cols) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse=", "))
}

# sparse group design matrix: cells x groups
.group_design <- function(groups) {
  groups <- as.factor(groups)
  G <- Matrix::sparse.model.matrix(~0 + groups)
  # keep group labels consistent with factor levels (avoid "groups" prefix leakage)
  colnames(G) <- levels(groups)
  G
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

# choose safe nfolds for cv.glmnet based on sample size
.safe_nfolds <- function(n, max_folds = 5, min_folds = 3) {
  if (is.na(n) || n < min_folds) return(NA_integer_)
  as.integer(max(min_folds, min(max_folds, n)))
}


# cross-platform parallel lapply (PSOCK on all OS)
.parallel_lapply <- function(X, FUN, n_cores = 1, seed = NULL) {
  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores <= 1 || length(X) <= 1) return(lapply(X, FUN))
  n_cores <- min(n_cores, length(X))
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # improve cross-platform reproducibility when randomness exists in workers
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = as.integer(seed))

  parallel::parLapply(cl, X, FUN)
}
