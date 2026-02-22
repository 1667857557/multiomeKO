# ---- prediction / virtual KO ----

#' Predict virtual KO counterfactual for given regulators
#'
#' For Poisson models, predictions are expected counts (mu). Deltas are computed
#' on log1p(mu) scale for numerical stability and interpretability.
#'
#' @param TAX output of extract_TAX_metacell
#' @param fit1 stage1 fit list(W_A,b0_A)
#' @param fit2 stage2 fit list(V,W_X,b0_X)
#' @param ko_regulators character vector of regulator rownames(TAX$T)
#' @param ko_value value after KO (0 by default)
#' @param ko_mode "set" sets regulator to ko_value; "scale" multiplies by ko_value
#' @return list(muA_wt, muA_cf, dA, muX_wt, muX_cf, dX)
predict_virtual_ko <- function(TAX, fit1, fit2, ko_regulators, ko_value = 0, ko_mode = c("set", "scale")) {
  ko_mode <- match.arg(ko_mode)
  T <- TAX$T
  if (!all(ko_regulators %in% rownames(T))) {
    stop("Some ko_regulators not found in TAX$T: ", paste(setdiff(ko_regulators, rownames(T)), collapse=", "))
  }

  T_cf <- T
  if (ko_mode == "set") {
    T_cf[ko_regulators, ] <- ko_value
  } else {
    T_cf[ko_regulators, ] <- T_cf[ko_regulators, , drop=FALSE] * ko_value
  }

  offA <- TAX$offsets$atac
  offX <- TAX$offsets$rna

  # ---- Stage1: muA = exp(b0 + W^T T + offset) ----
  W_A <- fit1$W_A
  b0A <- .align_vec(fit1$b0_A, colnames(W_A))

  etaA_wt <- (Matrix::t(W_A) %*% T)              # peaks x metacells
  etaA_cf <- (Matrix::t(W_A) %*% T_cf)

  # add intercept (recycle across columns)
  etaA_wt <- etaA_wt + b0A
  etaA_cf <- etaA_cf + b0A

  # add offsets per metacell (column-wise)
  etaA_wt <- t(t(as.matrix(etaA_wt)) + offA)
  etaA_cf <- t(t(as.matrix(etaA_cf)) + offA)

  muA_wt <- exp(etaA_wt)
  muA_cf <- exp(etaA_cf)

  dA <- log1p(muA_cf) - log1p(muA_wt)

  # Stage2 features: log1p(muA)
  Afeat_wt <- log1p(muA_wt)
  Afeat_cf <- log1p(muA_cf)

  # ---- Stage2: muX = exp(b0 + V^T Afeat + W_X^T T + offset) ----
  V <- fit2$V
  genes <- colnames(V)

  b0X <- .align_vec(fit2$b0_X, genes)

  etaX_wt <- (Matrix::t(V) %*% Matrix::Matrix(Afeat_wt, sparse = FALSE))
  etaX_cf <- (Matrix::t(V) %*% Matrix::Matrix(Afeat_cf, sparse = FALSE))

  # add direct regulator term if present
  if (!is.null(fit2$W_X)) {
    etaX_wt <- etaX_wt + (Matrix::t(fit2$W_X) %*% T)
    etaX_cf <- etaX_cf + (Matrix::t(fit2$W_X) %*% T_cf)
  }

  etaX_wt <- etaX_wt + b0X
  etaX_cf <- etaX_cf + b0X

  etaX_wt <- t(t(as.matrix(etaX_wt)) + offX)
  etaX_cf <- t(t(as.matrix(etaX_cf)) + offX)

  muX_wt <- exp(etaX_wt)
  muX_cf <- exp(etaX_cf)

  dX <- log1p(muX_cf) - log1p(muX_wt)

  list(
    muA_wt = muA_wt, muA_cf = muA_cf, dA = dA,
    muX_wt = muX_wt, muX_cf = muX_cf, dX = dX
  )
}

#' Rank features by average absolute effect across metacells
#' @param dmat matrix features x metacells (delta on log1p scale)
#' @param top_n number to return
#' @return named numeric vector sorted
rank_by_effect <- function(dmat, top_n = 50) {
  if (inherits(dmat, "dgCMatrix")) {
    s <- Matrix::rowMeans(abs(dmat))
  } else {
    s <- rowMeans(abs(dmat))
  }
  s <- sort(s, decreasing = TRUE)
  if (!is.null(top_n)) s <- head(s, top_n)
  s
}
