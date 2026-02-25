# ---- end-to-end optimized runner ----

#' Run optimized two-stage virtual KO on a Seurat multiome object
#'
#' Objective logic checks:
#' - Uses shared metacells for RNA and ATAC (paired structure).
#' - Fits Poisson ridge with library-size offsets (depth confounding control).
#' - Uses biologically grounded sparsity (motif for TF->peak; LinkPeaks for peak->gene).
#' - Adds adaptive penalties (GWAS PIP + link strength + distance) without breaking interpretability.
#'
#' @param obj Seurat object (RNA+ATAC; optional chromvar assay)
#' @param ko_regulator single regulator in rownames(TAX$T)
#' @param group.by metadata columns for metacell stratification
#' @param n_cells metacell size
#' @param genes_use genes to model (must include target genes of interest)
#' @param tf_mode RNA/chromvar/both
#' @param tf_genes TF genes if needed
#' @param include_T_direct include direct regulator terms in stage2 (default FALSE)
#' @param covariates optional metadata covariates to aggregate; default NULL (avoid depth double-correction)
#' @param stage1_mode prior mode for stage1: motif/chip/hybrid
#' @param chip_peak_map optional named list for chip/hybrid mode (regulator -> peaks)
#' @param ko_value KO value (used by predict_virtual_ko)
#' @param alpha KO effect strength in [0,1]; or 'auto' using externally estimated knockdown fraction
#' @param alpha_grid optional numeric grid in [0,1] to report uncertainty / nearest-grid alpha
#' @param alpha_external optional externally estimated knockdown fraction in [0,1] for alpha='auto'
#' @param allow_motif_ko whether to allow perturbing MOTIF:* regulators (default FALSE)
#' @param couple_rna_motif whether to synchronize RNA KO to matching MOTIF features (requires allow_motif_ko=TRUE; proxy TF-activity coupling)
#' @param ko_mode KO mode: 'set' or 'scale'
#' @param finemap_snps NULL, data.frame(chr,pos,pip), GRanges, or path
#' @param pip_floor minimum PIP
#' @param lambda_A,lambda_X optional; if NULL estimated on subset
#' @param n_cores number of CPU cores for parallelizable loops
#' @param seed random seed
#' @return list(setup, fits, diagnostics, pred, rank_genes, rank_peaks)
run_virtual_ko_optimized <- function(
  obj,
  ko_regulator,
  group.by = c("celltype", "SampleID"),
  n_cells = 50,
  genes_use,
  tf_mode = c("RNA", "chromvar", "both"),
  tf_genes = NULL,
  include_T_direct = FALSE,
  covariates = NULL,
  stage1_mode = c("motif", "chip", "hybrid"),
  chip_peak_map = NULL,
  ko_value = 0,
  alpha = 1,
  alpha_grid = NULL,
  alpha_external = NULL,
  allow_motif_ko = FALSE,
  couple_rna_motif = FALSE,
  ko_mode = c("set", "scale"),
  finemap_snps = NULL,
  pip_floor = 0,
  lambda_A = NULL,
  lambda_X = NULL,
  n_cores = 1,
  seed = 1,
  rna_assay = "RNA",
  atac_assay = "ATAC"
) {
  tf_mode <- match.arg(tf_mode)
  ko_mode <- match.arg(ko_mode)
  stage1_mode <- match.arg(stage1_mode)

  if (!startsWith(ko_regulator, "RNA:") && !isTRUE(allow_motif_ko)) {
    stop("KO semantics are restricted to RNA:* by default. For motif perturbation, set allow_motif_ko=TRUE.")
  }
  if (isTRUE(couple_rna_motif) && !isTRUE(allow_motif_ko)) {
    stop("couple_rna_motif=TRUE requires allow_motif_ko=TRUE.")
  }

  if (ko_mode == "set" && (!is.null(alpha_grid) || (is.character(alpha) && alpha == "auto") || (is.numeric(alpha) && length(alpha) == 1 && !isTRUE(all.equal(alpha, 1))))) {
    warning("alpha/alpha_grid are only used when ko_mode='scale'; current run uses ko_mode='set'.")
  }

  # ---- priors: peak2gene ----
  p2g <- get_peak2gene_links(obj, atac_assay = atac_assay)
  genes_use <- intersect(genes_use, unique(p2g$gene))
  if (length(genes_use) == 0) stop("genes_use not found in Links gene column")
  p2g <- p2g[p2g$gene %in% genes_use, , drop=FALSE]
  peaks_use <- unique(p2g$peak)

  # ---- GWAS weights (optional) ----
  peak_w <- NULL
  snp_peak <- NULL
  if (!is.null(finemap_snps)) {
    if (is.character(finemap_snps) && length(finemap_snps) == 1 && file.exists(finemap_snps)) {
      snp_df <- read_finemap_snps(finemap_snps)
      snps_gr <- snps_to_granges(snp_df)
    } else if (inherits(finemap_snps, "GRanges")) {
      snps_gr <- finemap_snps
    } else if (is.data.frame(finemap_snps)) {
      if (!all(c("chr","pos","pip") %in% colnames(finemap_snps))) {
        stop("finemap_snps data.frame must have chr,pos,pip")
      }
      snps_gr <- snps_to_granges(finemap_snps)
    } else {
      stop("Unsupported finemap_snps type")
    }
    snp_peak <- map_snps_to_peaks(obj, snps_gr, atac_assay = atac_assay)
    peak_w <- peak_weights_from_snps(snp_peak, method = "max", pip_floor = pip_floor)
  }

  # ---- metacells ----
  mc <- make_metacell_ids(obj, group.by = group.by, n_cells = n_cells, seed = seed)

  # ---- stage1 priors (pluggable) ----
  stage1_prior <- build_stage1_priors(
    obj,
    atac_assay = atac_assay,
    stage1_mode = stage1_mode,
    chip_peak_map = chip_peak_map
  )
  motif_map <- if (!is.null(stage1_prior$motif_map)) stage1_prior$motif_map else NULL
  if (isTRUE(include_T_direct) && is.null(motif_map)) {
    stop("include_T_direct=TRUE requires motif priors; use stage1_mode='motif' or 'hybrid'.")
  }

  # ---- extract shared T/A/X ----
  if (tf_mode %in% c("RNA", "both") && is.null(tf_genes)) {
    rna_genes <- rownames(Seurat::GetAssayData(obj, assay = rna_assay, slot = "counts"))
    tf_from_motif <- if (!is.null(motif_map)) intersect(unique(motif_map$tf_symbols), rna_genes) else character()
    tf_genes <- unique(c(tf_from_motif, genes_use))
    warning("tf_genes is NULL; auto-derived RNA regulators from available motif TFs + genes_use")
  }

  TAX <- extract_TAX_metacell(
    obj,
    metacell_ids = mc,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    genes_use = genes_use,
    peaks_use = peaks_use,
    tf_mode = tf_mode,
    tf_genes = tf_genes,
    covariates = covariates
  )

  if (!(ko_regulator %in% rownames(TAX$T))) {
    stop("ko_regulator not found in extracted regulators. Available examples: ",
         paste(head(rownames(TAX$T), 10), collapse=", "))
  }


  # ---- lambda estimation (optional) ----
  if (is.null(lambda_A)) {
    if (!is.null(motif_map)) {
      lambda_A <- estimate_lambda_stage1(TAX, motif_map, peaks = peaks_use, seed = seed, n_cores = n_cores)
    } else {
      lambda_A <- 1
      warning("lambda_A fallback to 1 in chip-only stage1 mode")
    }
  }
  if (is.null(lambda_X)) {
    lambda_X <- estimate_lambda_stage2(TAX, p2g_df = p2g, genes = genes_use, seed = seed, n_cores = n_cores)
  }

  # ---- fit models ----
  fit1 <- fit_stage1_T_to_A(TAX, stage1_prior, peaks = peaks_use, lambda = lambda_A, n_cores = n_cores)
  fit2 <- fit_stage2_A_to_X(
    TAX, p2g_df = p2g, genes = genes_use, lambda = lambda_X,
    obj = obj, atac_assay = atac_assay,
    peak_weights = peak_w,
    include_T_direct = include_T_direct,
    motif_map = motif_map,
    n_cores = n_cores
  )

  # ---- alpha calibration for CRISPRi-strength alignment ----
  alpha_used <- 1
  alpha_scan <- NULL
  if (ko_mode == "scale") {
    if (is.numeric(alpha) && length(alpha) == 1 && is.finite(alpha)) {
      alpha_used <- alpha
      if (alpha_used < 0 || alpha_used > 1) stop("alpha must be within [0,1]")
      if (!is.null(alpha_grid)) {
        a <- unique(as.numeric(alpha_grid))
        a <- a[is.finite(a) & a >= 0 & a <= 1]
        if (length(a) > 0) alpha_scan <- data.frame(alpha = a, loss = abs(a - alpha_used))
      }
    } else if (is.character(alpha) && length(alpha) == 1 && alpha == "auto") {
      if (is.null(alpha_external) || !is.numeric(alpha_external) || length(alpha_external) != 1 || !is.finite(alpha_external)) {
        stop("alpha='auto' requires alpha_external (externally estimated knockdown fraction in [0,1])")
      }
      if (alpha_external < 0 || alpha_external > 1) stop("alpha_external must be within [0,1]")
      alpha_used <- alpha_external
      if (!is.null(alpha_grid)) {
        a <- unique(as.numeric(alpha_grid))
        a <- a[is.finite(a) & a >= 0 & a <= 1]
        if (length(a) > 0) {
          alpha_scan <- data.frame(alpha = a, loss = abs(a - alpha_external))
          alpha_used <- a[which.min(abs(a - alpha_external))]
        }
      }
    } else {
      stop("alpha must be numeric scalar in [0,1], or use alpha='auto' with alpha_external")
    }
  }

  # ---- predict KO ----
  pred <- predict_virtual_ko(
    TAX, fit1, fit2,
    ko_regulators = ko_regulator,
    ko_value = if (ko_mode == "scale") alpha_used else ko_value,
    ko_mode = ko_mode,
    allow_motif_ko = allow_motif_ko,
    couple_rna_motif = couple_rna_motif
  )

  rank_genes <- rank_by_effect(pred$dX, top_n = 100)
  rank_peaks <- rank_by_effect(pred$dA, top_n = 100)

  # ---- fit diagnostics warning ----
  d1 <- fit1$diagnostics
  d2 <- fit2$diagnostics
  if (!is.null(d1) && d1$total_peaks > 0) {
    frac1 <- d1$fitted_peaks / d1$total_peaks
    if (is.finite(frac1) && frac1 < 0.5) warning(sprintf("Stage1 fitted peak fraction is low: %.2f", frac1))
  }
  if (!is.null(d2) && d2$total_genes > 0) {
    frac2 <- d2$fitted_genes / d2$total_genes
    if (is.finite(frac2) && frac2 < 0.5) warning(sprintf("Stage2 fitted gene fraction is low: %.2f", frac2))
  }

  list(
    setup = list(
      group.by = group.by,
      n_cells = n_cells,
      genes_use = genes_use,
      peaks_use = peaks_use,
      tf_mode = tf_mode,
      covariates = covariates,
      stage1_mode = stage1_mode,
      lambda_A = lambda_A,
      lambda_X = lambda_X,
      n_cores = n_cores,
      ko_regulator = ko_regulator,
      ko_value = ko_value,
      ko_mode = ko_mode,
      allow_motif_ko = allow_motif_ko,
      couple_rna_motif = couple_rna_motif,
      alpha = alpha_used
    ),
    gwas = list(snp_peak = snp_peak, peak_weights = peak_w),
    fits = list(stage1 = fit1, stage2 = fit2),
    diagnostics = list(stage1 = fit1$diagnostics, stage2 = fit2$diagnostics, alpha_scan = alpha_scan),
    pred = pred,
    rank_genes = rank_genes,
    rank_peaks = rank_peaks
  )
}
