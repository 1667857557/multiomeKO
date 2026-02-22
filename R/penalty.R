# ---- adaptive penalties for interpretability + accuracy ----

# helper: get gene TSS from annotation
.get_gene_tss <- function(obj, gene) {
  ann <- tryCatch(Signac::Annotation(obj), error = function(e) NULL)
  if (is.null(ann) || length(ann) == 0) return(NULL)
  df <- as.data.frame(ann)
  # find gene name column
  gene_col <- NULL
  for (cand in c("gene_name", "gene", "symbol", "name")) {
    if (cand %in% names(df)) { gene_col <- cand; break }
  }
  if (is.null(gene_col)) return(NULL)
  idx <- which(as.character(df[[gene_col]]) == gene)
  if (length(idx) == 0) return(NULL)
  # Use strand-aware TSS if strand present
  if ("strand" %in% names(df)) {
    # take first transcript (could be multiple)
    s <- as.character(df$strand[idx[1]])
    if (s == "-" && "end" %in% names(df)) return(df$end[idx[1]])
  }
  if ("start" %in% names(df)) return(df$start[idx[1]])
  NULL
}

#' Build stage2 penalty.factor for peaks linked to a gene
#'
#' Combines three biologically grounded evidences:
#' 1) GWAS fine-mapping: peaks harboring high-PIP variants are more plausible causal CREs.
#' 2) Peak2gene link strength: stronger peak-gene association is more plausible.
#' 3) Distance-to-TSS prior: nearer elements have higher baseline probability (soft, not hard cutoff).
#'
#' Implementation: convert evidences to a [0,1]-like score and map to penalty via exp(-alpha*score).
#'
#' @param obj Seurat object
#' @param gene gene symbol
#' @param peaks character vector predictors
#' @param p2g_df data.frame from get_peak2gene_links (peak,gene, optional score)
#' @param peak_weights named numeric (GWAS peak weights in [0,1] typical) or NULL
#' @param alpha_total overall strength
#' @param min_penalty lower bound
#' @param dist_scale distance scale in bp for exp(-d/scale)
#' @return penalty.factor aligned to peaks
build_stage2_peak_penalty <- function(
  obj,
  atac_assay = "ATAC",
  gene,
  peaks,
  p2g_df,
  peak_weights = NULL,
  alpha_total = 2,
  min_penalty = 0.2,
  dist_scale = 50000
) {
  peaks <- as.character(peaks)
  # ---- 1) GWAS weight ----
  w_gwas <- rep(0, length(peaks)); names(w_gwas) <- peaks
  if (!is.null(peak_weights) && length(peak_weights) > 0) {
    idx <- intersect(names(peak_weights), peaks)
    w_gwas[idx] <- peak_weights[idx]
  }

  # ---- 2) link score ----
  w_link <- rep(0, length(peaks)); names(w_link) <- peaks
  if (!is.null(p2g_df) && nrow(p2g_df) > 0) {
    sub <- p2g_df[p2g_df$gene == gene & p2g_df$peak %in% peaks, , drop=FALSE]
    if (nrow(sub) > 0 && "score" %in% colnames(sub)) {
      s <- as.numeric(sub$score)
      s <- .finite0(s, 0)
      # robust min-max scaling
      if (length(unique(s)) > 1) s <- (s - min(s)) / (max(s) - min(s)) else s <- rep(0.5, length(s))
      w_link[sub$peak] <- s
    } else if (nrow(sub) > 0) {
      # if no score, use binary link
      w_link[sub$peak] <- 0.5
    }
  }

  # ---- 3) distance to TSS ----
  w_dist <- rep(0, length(peaks)); names(w_dist) <- peaks
  tss <- .get_gene_tss(obj, gene)
  if (!is.null(tss)) {
    peak_gr <- Signac::granges(obj[[atac_assay]])
    # subset peaks present
    pk <- intersect(peaks, names(peak_gr))
    if (length(pk) > 0) {
      gr <- peak_gr[pk]
      mid <- floor((IRanges::start(gr) + IRanges::end(gr)) / 2)
      d <- abs(mid - tss)
      w_dist[pk] <- exp(-d / dist_scale)
    }
  }

  # combine evidences (bounded-ish)
  score <- (w_gwas + w_link + w_dist) / 3
  score <- pmin(1, pmax(0, score))

  pen <- exp(-alpha_total * score)
  pen <- pmax(min_penalty, pen)
  unname(pen)
}
