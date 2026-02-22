# ---- model fitting ----

#' Estimate a global lambda for stage1 (TF->ATAC) by sampling peaks
#' @param TAX output of extract_TAX_metacell
#' @param motif_map output of get_peak_motif_map
#' @param peaks vector of peaks to consider
#' @param n_sample number of peaks to sample
#' @param alpha ridge=0
#' @return numeric lambda
estimate_lambda_stage1 <- function(TAX, motif_map, peaks, n_sample = 50, alpha = 0, seed = 1) {
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
  # ensure peaks are rows
  if (!all(peaks_s %in% rownames(mm)) && all(peaks_s %in% colnames(mm))) mm <- t(mm)
  mm <- mm[intersect(peaks_s, rownames(mm)), , drop=FALSE]

  lambdas <- c()
  for (pk in rownames(mm)) {
    # candidate regulators based on motif hits
    hit <- which(mm[pk, ] > 0)
    if (length(hit) < 1) next
    motifs <- colnames(mm)[hit]
    tf_sym <- sub(".*_", "", motifs)

    regs <- character()
    # RNA: TF symbols
    rna_regs <- paste0("RNA:", tf_sym)
    regs <- c(regs, intersect(rna_regs, rownames(T)))
    # MOTIF: motif names
    motif_regs <- paste0("MOTIF:", motifs)
    regs <- c(regs, intersect(motif_regs, rownames(T)))
    regs <- unique(regs)

    if (length(regs) < 3) next

    X <- Xm[, regs, drop=FALSE]
    if (!is.null(Cm)) {
      X <- cbind(X, Cm)
      pen <- c(rep(1, length(regs)), rep(0, ncol(Cm)))
    } else {
      pen <- rep(1, length(regs))
    }
    y <- as.numeric(TAX$A_counts[pk, ])

    cv <- glmnet::cv.glmnet(
      x = as.matrix(X), y = y, family = "poisson",
      alpha = alpha, offset = off,
      penalty.factor = pen,
      nfolds = min(5, nrow(X)),
      standardize = TRUE
    )
    lambdas <- c(lambdas, cv$lambda.1se)
  }

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
#' @return numeric lambda
estimate_lambda_stage2 <- function(TAX, p2g_df, genes, n_sample = 50, seed = 1) {
  set.seed(seed)
  genes <- intersect(genes, rownames(TAX$X_counts))
  if (length(genes) == 0) stop("No genes available for lambda estimation")
  genes_s <- sample(genes, min(n_sample, length(genes)))

  Afeat <- TAX$A_log1p
  off <- TAX$offsets$rna

  lambdas <- c()
  for (g in genes_s) {
    pks <- unique(p2g_df$peak[p2g_df$gene == g])
    pks <- intersect(pks, rownames(Afeat))
    if (length(pks) < 5) next

    X <- t(Afeat[pks, , drop=FALSE])
    y <- as.numeric(TAX$X_counts[g, ])

    cv <- glmnet::cv.glmnet(
      x = as.matrix(X), y = y, family = "poisson",
      alpha = 0, offset = off,
      nfolds = min(5, nrow(X)),
      standardize = TRUE
    )
    lambdas <- c(lambdas, cv$lambda.1se)
  }

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
#' @return list(W_A, b0_A)
fit_stage1_T_to_A <- function(TAX, motif_map, peaks, lambda = 1) {
  T <- TAX$T
  C <- TAX$covariates
  Xm <- t(T)
  Cm <- if (!is.null(C)) t(C) else NULL
  off <- TAX$offsets$atac

  mm <- motif_map$motif_mat
  # ensure peaks are rows
  if (!all(peaks %in% rownames(mm)) && all(peaks %in% colnames(mm))) mm <- t(mm)
  peaks <- intersect(peaks, rownames(mm))
  peaks <- intersect(peaks, rownames(TAX$A_counts))
  if (length(peaks) == 0) stop("No peaks to fit in stage1")

  W <- Matrix::Matrix(0, nrow = nrow(T), ncol = length(peaks), sparse = TRUE,
                      dimnames = list(rownames(T), peaks))
  b0 <- setNames(rep(0, length(peaks)), peaks)

  for (pk in peaks) {
    hit <- which(mm[pk, ] > 0)
    if (length(hit) < 1) next
    motifs <- colnames(mm)[hit]
    tf_sym <- sub(".*_", "", motifs)

    regs <- character()
    regs <- c(regs, intersect(paste0("RNA:", tf_sym), rownames(T)))
    regs <- c(regs, intersect(paste0("MOTIF:", motifs), rownames(T)))
    regs <- unique(regs)

    if (length(regs) < 3) next

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
    # cf: sparse matrix with rownames (Intercept) + predictors
    b0[pk] <- as.numeric(cf[1, 1])
    # coefficients for regs (exclude covariates)
    if (length(regs) > 0) {
      W[regs, pk] <- as.numeric(cf[match(regs, rownames(cf)), 1])
    }
  }

  list(W_A = W, b0_A = b0)
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
#' @return list(V, W_X, b0_X)
fit_stage2_A_to_X <- function(
  TAX, p2g_df, genes, lambda = 1,
  obj, atac_assay = "ATAC",
  peak_weights = NULL,
  include_T_direct = FALSE,
  motif_map = NULL
) {
  Afeat <- TAX$A_log1p
  off <- TAX$offsets$rna
  T <- TAX$T
  Xm_T <- t(T)
  C <- TAX$covariates
  Cm <- if (!is.null(C)) t(C) else NULL

  genes <- intersect(genes, rownames(TAX$X_counts))
  if (length(genes) == 0) stop("No genes to fit in stage2")

  # collect all peaks needed
  peaks_all <- unique(p2g_df$peak[p2g_df$gene %in% genes])
  peaks_all <- intersect(peaks_all, rownames(Afeat))
  if (length(peaks_all) == 0) stop("No linked peaks found for genes")

  V <- Matrix::Matrix(0, nrow = length(peaks_all), ncol = length(genes), sparse = TRUE,
                      dimnames = list(peaks_all, genes))
  b0 <- setNames(rep(0, length(genes)), genes)

  Wxt <- NULL
  if (include_T_direct) {
    Wxt <- Matrix::Matrix(0, nrow = nrow(T), ncol = length(genes), sparse = TRUE,
                          dimnames = list(rownames(T), genes))
    if (is.null(motif_map)) stop("motif_map required when include_T_direct=TRUE")
  }

  mm <- NULL
  if (include_T_direct) {
    mm <- motif_map$motif_mat
    if (!all(peaks_all %in% rownames(mm)) && all(peaks_all %in% colnames(mm))) mm <- t(mm)
  }

  for (g in genes) {
    pks <- unique(p2g_df$peak[p2g_df$gene == g])
    pks <- intersect(pks, peaks_all)
    if (length(pks) < 5) next

    Xp <- t(Afeat[pks, , drop=FALSE]) # metacells x peaks
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
      # Candidate regulators: those with motifs present in linked peaks
      hit <- which(Matrix::colSums(mm[pks, , drop=FALSE] > 0) > 0)
      motifs <- colnames(mm)[hit]
      tf_sym <- sub(".*_", "", motifs)
      reg_keep <- unique(c(
        intersect(paste0("RNA:", tf_sym), rownames(T)),
        intersect(paste0("MOTIF:", motifs), rownames(T))
      ))
      # keep at most 50 to avoid overfitting
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
    b0[g] <- as.numeric(cf[1, 1])

    # peaks coefficients
    V[pks, g] <- as.numeric(cf[match(colnames(Xp), rownames(cf)), 1])

    # direct regulator coefficients
    if (include_T_direct && length(reg_keep) > 0) {
      Wxt[reg_keep, g] <- as.numeric(cf[match(reg_keep, rownames(cf)), 1])
    }
  }

  list(V = V, W_X = Wxt, b0_X = b0)
}
