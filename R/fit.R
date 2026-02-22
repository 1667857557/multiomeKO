# ---- model fitting ----

#' Estimate a global lambda for stage1 (TF->ATAC) by sampling peaks
#' @param TAX output of extract_TAX_metacell
#' @param motif_map output of get_peak_motif_map
#' @param peaks vector of peaks to consider
#' @param n_sample number of peaks to sample
#' @param alpha ridge=0
#' @param n_cores number of cores for parallel loop
#' @return numeric lambda
estimate_lambda_stage1 <- function(TAX, motif_map, peaks, n_sample = 50, alpha = 0, seed = 1, n_cores = 1) {
  set.seed(seed)
  peaks <- intersect(peaks, rownames(TAX$A_counts))
  if (length(peaks) == 0) stop("No peaks available for lambda estimation")
  peaks_s <- sample(peaks, min(n_sample, length(peaks)))

  T <- TAX$T
  C <- TAX$covariates
  Xm <- t(T)
  Cm <- if (!is.null(C)) t(C) else NULL
  off <- TAX$offsets$atac

  mm <- motif_map$motif_mat
  if (!all(peaks_s %in% rownames(mm)) && all(peaks_s %in% colnames(mm))) mm <- t(mm)
  mm <- mm[intersect(peaks_s, rownames(mm)), , drop=FALSE]

  lam_list <- .parallel_lapply(rownames(mm), n_cores = n_cores, FUN = function(pk) {
    hit <- which(mm[pk, ] > 0)
    if (length(hit) < 1) return(NA_real_)
    motifs <- colnames(mm)[hit]
    tf_sym <- sub(".*_", "", motifs)

    regs <- unique(c(
      intersect(paste0("RNA:", tf_sym), rownames(T)),
      intersect(paste0("MOTIF:", motifs), rownames(T))
    ))
    if (length(regs) < 3) return(NA_real_)

    X <- Xm[, regs, drop=FALSE]
    if (!is.null(Cm)) {
      X <- cbind(X, Cm)
      pen <- c(rep(1, length(regs)), rep(0, ncol(Cm)))
    } else {
      pen <- rep(1, length(regs))
    }
    y <- as.numeric(TAX$A_counts[pk, ])

    nfolds <- .safe_nfolds(nrow(X))
    if (is.na(nfolds)) return(NA_real_)

    cv <- glmnet::cv.glmnet(
      x = as.matrix(X), y = y, family = "poisson",
      alpha = alpha, offset = off,
      penalty.factor = pen,
      nfolds = nfolds,
      standardize = TRUE
    )
    as.numeric(cv$lambda.1se)
  })

  lambdas <- unlist(lam_list, use.names = FALSE)
  lambdas <- lambdas[is.finite(lambdas)]
  if (length(lambdas) == 0) {
    warning("Lambda estimation failed; fallback to 1")
    return(1)
  }
  stats::median(lambdas)
}

#' Estimate a global lambda for stage2 (ATAC->RNA) by sampling genes
#' @param TAX output of extract_TAX_metacell
#' @param p2g_df peak2gene links
#' @param genes vector of genes
#' @param n_sample number of genes to sample
#' @param n_cores number of cores for parallel loop
#' @return numeric lambda
estimate_lambda_stage2 <- function(TAX, p2g_df, genes, n_sample = 50, seed = 1, n_cores = 1) {
  set.seed(seed)
  genes <- intersect(genes, rownames(TAX$X_counts))
  if (length(genes) == 0) stop("No genes available for lambda estimation")
  genes_s <- sample(genes, min(n_sample, length(genes)))

  Afeat <- TAX$A_log1p
  off <- TAX$offsets$rna

  lam_list <- .parallel_lapply(genes_s, n_cores = n_cores, FUN = function(g) {
    pks <- unique(p2g_df$peak[p2g_df$gene == g])
    pks <- intersect(pks, rownames(Afeat))
    if (length(pks) < 5) return(NA_real_)

    X <- t(Afeat[pks, , drop=FALSE])
    y <- as.numeric(TAX$X_counts[g, ])

    nfolds <- .safe_nfolds(nrow(X))
    if (is.na(nfolds)) return(NA_real_)

    cv <- glmnet::cv.glmnet(
      x = as.matrix(X), y = y, family = "poisson",
      alpha = 0, offset = off,
      nfolds = nfolds,
      standardize = TRUE
    )
    as.numeric(cv$lambda.1se)
  })

  lambdas <- unlist(lam_list, use.names = FALSE)
  lambdas <- lambdas[is.finite(lambdas)]
  if (length(lambdas) == 0) {
    warning("Lambda estimation failed; fallback to 1")
    return(1)
  }
  stats::median(lambdas)
}

#' Fit stage1: TF/Regulator -> peak accessibility (Poisson ridge)
#' @param TAX output of extract_TAX_metacell
#' @param motif_map output of get_peak_motif_map
#' @param peaks peaks to fit
#' @param lambda global lambda
#' @param n_cores number of cores for parallel loop
#' @return list(W_A, b0_A, diagnostics)
fit_stage1_T_to_A <- function(TAX, motif_map, peaks, lambda = 1, n_cores = 1) {
  T <- TAX$T
  C <- TAX$covariates
  Xm <- t(T)
  Cm <- if (!is.null(C)) t(C) else NULL
  off <- TAX$offsets$atac

  mm <- motif_map$motif_mat
  if (!all(peaks %in% rownames(mm)) && all(peaks %in% colnames(mm))) mm <- t(mm)
  peaks <- intersect(peaks, rownames(mm))
  peaks <- intersect(peaks, rownames(TAX$A_counts))
  if (length(peaks) == 0) stop("No peaks to fit in stage1")

  res_list <- .parallel_lapply(peaks, n_cores = n_cores, FUN = function(pk) {
    hit <- which(mm[pk, ] > 0)
    if (length(hit) < 1) {
      b0_null <- .poisson_intercept_with_offset(TAX$A_counts[pk, ], off)
      return(list(pk = pk, b0 = b0_null, regs = character(), coefs = numeric(), status = "no_motif"))
    }

    motifs <- colnames(mm)[hit]
    tf_sym <- sub(".*_", "", motifs)
    regs <- unique(c(
      intersect(paste0("RNA:", tf_sym), rownames(T)),
      intersect(paste0("MOTIF:", motifs), rownames(T))
    ))
    if (length(regs) < 3) {
      b0_null <- .poisson_intercept_with_offset(TAX$A_counts[pk, ], off)
      return(list(pk = pk, b0 = b0_null, regs = character(), coefs = numeric(), status = "few_regs"))
    }

    X <- Xm[, regs, drop=FALSE]
    pen <- rep(1, length(regs))
    if (!is.null(Cm)) {
      X <- cbind(X, Cm)
      pen <- c(pen, rep(0, ncol(Cm)))
    }

    y <- as.numeric(TAX$A_counts[pk, ])
    fit <- glmnet::glmnet(
      x = as.matrix(X), y = y, family = "poisson",
      alpha = 0, lambda = lambda,
      offset = off,
      penalty.factor = pen,
      standardize = TRUE
    )

    cf <- glmnet::coef.glmnet(fit)
    coef_vals <- as.numeric(cf[match(regs, rownames(cf)), 1])
    list(pk = pk, b0 = as.numeric(cf[1, 1]), regs = regs, coefs = coef_vals, status = "fit")
  })

  W <- Matrix::Matrix(0, nrow = nrow(T), ncol = length(peaks), sparse = TRUE,
                      dimnames = list(rownames(T), peaks))
  b0 <- setNames(rep(0, length(peaks)), peaks)
  fitted_peaks <- 0L
  skipped_no_motif <- 0L
  skipped_few_regs <- 0L

  for (x in res_list) {
    b0[x$pk] <- x$b0
    if (x$status == "fit") {
      if (length(x$regs) > 0) W[x$regs, x$pk] <- x$coefs
      fitted_peaks <- fitted_peaks + 1L
    } else if (x$status == "no_motif") {
      skipped_no_motif <- skipped_no_motif + 1L
    } else if (x$status == "few_regs") {
      skipped_few_regs <- skipped_few_regs + 1L
    }
  }

  list(
    W_A = W,
    b0_A = b0,
    diagnostics = list(
      total_peaks = length(peaks),
      fitted_peaks = fitted_peaks,
      skipped_no_motif = skipped_no_motif,
      skipped_few_regs = skipped_few_regs
    )
  )
}

#' Fit stage2: peaks -> gene expression (Poisson ridge) with adaptive penalties
#' @param TAX output of extract_TAX_metacell
#' @param p2g_df output of get_peak2gene_links
#' @param genes genes to fit
#' @param lambda global lambda
#' @param obj Seurat object (for distance-to-TSS)
#' @param atac_assay ATAC assay name
#' @param peak_weights optional named numeric GWAS peak weights
#' @param include_T_direct whether to include direct regulator terms T in stage2
#' @param motif_map optional motif_map; required if include_T_direct=TRUE
#' @param n_cores number of cores for parallel loop
#' @return list(V, W_X, b0_X, diagnostics)
fit_stage2_A_to_X <- function(
  TAX, p2g_df, genes, lambda = 1,
  obj, atac_assay = "ATAC",
  peak_weights = NULL,
  include_T_direct = FALSE,
  motif_map = NULL,
  n_cores = 1
) {
  Afeat <- TAX$A_log1p
  off <- TAX$offsets$rna
  T <- TAX$T
  Xm_T <- t(T)
  C <- TAX$covariates
  Cm <- if (!is.null(C)) t(C) else NULL

  genes <- intersect(genes, rownames(TAX$X_counts))
  if (length(genes) == 0) stop("No genes to fit in stage2")

  peaks_all <- unique(p2g_df$peak[p2g_df$gene %in% genes])
  peaks_all <- intersect(peaks_all, rownames(Afeat))
  if (length(peaks_all) == 0) stop("No linked peaks found for genes")

  mm <- NULL
  if (include_T_direct) {
    if (is.null(motif_map)) stop("motif_map required when include_T_direct=TRUE")
    mm <- motif_map$motif_mat
    if (!all(peaks_all %in% rownames(mm)) && all(peaks_all %in% colnames(mm))) mm <- t(mm)
  }

  res_list <- .parallel_lapply(genes, n_cores = n_cores, FUN = function(g) {
    pks <- unique(p2g_df$peak[p2g_df$gene == g])
    pks <- intersect(pks, peaks_all)
    if (length(pks) < 5) {
      b0_null <- .poisson_intercept_with_offset(TAX$X_counts[g, ], off)
      return(list(g = g, b0 = b0_null, pks = character(), v = numeric(), reg_keep = character(), w = numeric(), status = "few_peaks"))
    }

    Xp <- t(Afeat[pks, , drop=FALSE])
    pen_p <- build_stage2_peak_penalty(
      obj = obj, atac_assay = atac_assay,
      gene = g, peaks = pks,
      p2g_df = p2g_df,
      peak_weights = peak_weights,
      alpha_total = 2,
      min_penalty = 0.2
    )

    Xdesign <- Xp
    pen <- pen_p
    reg_keep <- character()

    if (include_T_direct) {
      hit <- which(Matrix::colSums(mm[pks, , drop=FALSE] > 0) > 0)
      motifs <- colnames(mm)[hit]
      tf_sym <- sub(".*_", "", motifs)
      reg_keep <- unique(c(
        intersect(paste0("RNA:", tf_sym), rownames(T)),
        intersect(paste0("MOTIF:", motifs), rownames(T))
      ))
      if (length(reg_keep) > 50) reg_keep <- reg_keep[seq_len(50)]

      if (length(reg_keep) > 0) {
        Xt <- Xm_T[, reg_keep, drop=FALSE]
        Xdesign <- cbind(Xdesign, Xt)
        pen <- c(pen, rep(1, length(reg_keep)))
      }
    }

    if (!is.null(Cm)) {
      Xdesign <- cbind(Xdesign, Cm)
      pen <- c(pen, rep(0, ncol(Cm)))
    }

    y <- as.numeric(TAX$X_counts[g, ])
    fit <- glmnet::glmnet(
      x = as.matrix(Xdesign), y = y, family = "poisson",
      alpha = 0, lambda = lambda,
      offset = off,
      penalty.factor = pen,
      standardize = TRUE
    )

    cf <- glmnet::coef.glmnet(fit)
    v <- as.numeric(cf[match(colnames(Xp), rownames(cf)), 1])
    w <- numeric()
    if (include_T_direct && length(reg_keep) > 0) {
      w <- as.numeric(cf[match(reg_keep, rownames(cf)), 1])
    }

    list(g = g, b0 = as.numeric(cf[1, 1]), pks = pks, v = v, reg_keep = reg_keep, w = w, status = "fit")
  })

  V <- Matrix::Matrix(0, nrow = length(peaks_all), ncol = length(genes), sparse = TRUE,
                      dimnames = list(peaks_all, genes))
  b0 <- setNames(rep(0, length(genes)), genes)

  Wxt <- NULL
  if (include_T_direct) {
    Wxt <- Matrix::Matrix(0, nrow = nrow(T), ncol = length(genes), sparse = TRUE,
                          dimnames = list(rownames(T), genes))
  }

  fitted_genes <- 0L
  skipped_few_peaks <- 0L
  for (x in res_list) {
    b0[x$g] <- x$b0
    if (x$status == "few_peaks") {
      skipped_few_peaks <- skipped_few_peaks + 1L
      next
    }
    if (length(x$pks) > 0) V[x$pks, x$g] <- x$v
    if (include_T_direct && length(x$reg_keep) > 0) Wxt[x$reg_keep, x$g] <- x$w
    fitted_genes <- fitted_genes + 1L
  }

  list(
    V = V,
    W_X = Wxt,
    b0_X = b0,
    diagnostics = list(
      total_genes = length(genes),
      fitted_genes = fitted_genes,
      skipped_few_peaks = skipped_few_peaks
    )
  )
}
