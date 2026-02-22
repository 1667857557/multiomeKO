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
#' @param ko_regulator single re (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/OPTIMIZATION_PLAN.md b/OPTIMIZATION_PLAN.md
new file mode 100644
index 0000000000000000000000000000000000000000..0d3e92172bd400cb196f78a98e9f5bcd18ff05d3
--- /dev/null
+++ b/OPTIMIZATION_PLAN.md
@@ -0,0 +1,42 @@
+# MultiomeKO 优化方案（可执行）
+
+## 目标
+在保留“两阶段 Regulator→ATAC→RNA”框架下，提升：
+1) 稳健性；2) 生物学可解释性；3) 结果可复现性。
+
+## 已完成的首批优化
+- **CV 折数安全化**：新增 `.safe_nfolds()`，在样本过少时跳过 CV，避免 `cv.glmnet` 边界失败。
+- **默认 TF 回退策略**：当 `tf_mode` 包含 RNA 且 `tf_genes=NULL` 时，自动回退到 `genes_use` 并发出 warning。
+- **KO 语义可配置**：`predict_virtual_ko()` 新增 `ko_mode`：
+  - `set`：设为常数（原逻辑）
+  - `scale`：按比例缩放（更接近“部分抑制”）
+- **拟合诊断输出**：stage1/stage2 返回 `diagnostics`，记录拟合与跳过数量，便于评估模型退化。
+
+## 下一步（建议按优先级）
+### P0（必须）
+1. **增加稳定性评估函数**
+   - 对 metacell 重采样（不同 seed）重复 KO，输出 rank 一致性（Spearman/Kendall）。
+2. **负对照机制**
+   - 增加随机 regulator KO 与低表达 regulator KO，对比真实 KO 的分布分离度。
+3. **显式报告退化比例**
+   - 若 `fitted_peaks/total_peaks` 或 `fitted_genes/total_genes` 过低，主函数给出 warning。
+
+### P1（重要）
+4. **不确定性区间**
+   - 对 `dX` 与 `rank_genes` 增加 bootstrap CI。
+5. **先验鲁棒映射层**
+   - 统一处理 Links/Motif 不同版本字段名，加入 schema 校验。
+6. **lambda 策略增强**
+   - 支持按基因分层 lambda（而非全局单值）并对比稳定性。
+
+### P2（增强）
+7. **生物知识约束增强**
+   - stage1 增加 TF family/group 稀疏结构惩罚（可选）。
+8. **可视化报告**
+   - 自动输出 KO 效应火山图、通路富集与关键 peak-gene-regulator 子网。
+
+## 验证清单
+- 至少 2 套公开 multiome 数据复现；
+- 正负对照分离显著；
+- 不同 `n_cells`、`lambda`、`ko_mode` 下结果一致性可接受；
+- 与已知 TF 靶基因方向一致。
 
EOF
)gulator id, e.g. 'RNA:FOXA1' or 'MOTIF:MA0148.1_FOXA1'
#' @param group.by metadata columns for metacell stratification
#' @param n_cells metacell size
#' @param genes_use genes to model (must include target genes of interest)
#' @param tf_mode RNA/chromvar/both
#' @param tf_genes TF genes if needed
#' @param include_T_direct include direct regulator terms in stage2 (default FALSE)
#' @param ko_value KO value (used by predict_virtual_ko)
#' @param ko_mode KO mode: 'set' or 'scale'
#' @param finemap_snps NULL, data.frame(chr,pos,pip), GRanges, or path
#' @param pip_floor minimum PIP
#' @param lambda_A,lambda_X optional; if NULL estimated on subset
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
  ko_value = 0,
  ko_mode = c("set", "scale"),
  finemap_snps = NULL,
  pip_floor = 0,
  lambda_A = NULL,
  lambda_X = NULL,
  seed = 1,
  rna_assay = "RNA",
  atac_assay = "ATAC"
) {
  tf_mode <- match.arg(tf_mode)
  ko_mode <- match.arg(ko_mode)

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

  # ---- motif map for stage1 ----
  motif_map <- get_peak_motif_map(obj, atac_assay = atac_assay)

  # ---- extract shared T/A/X ----
  if (tf_mode %in% c("RNA", "both") && is.null(tf_genes)) {
    rna_genes <- rownames(Seurat::GetAssayData(obj, assay = rna_assay, slot = "counts"))
    tf_from_motif <- intersect(unique(motif_map$tf_symbols), rna_genes)
    tf_genes <- unique(c(tf_from_motif, genes_use))
    warning("tf_genes is NULL; auto-derived RNA regulators from motifs + genes_use")
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
    covariates = c("nCount_RNA", "nCount_ATAC")
  )

  if (!(ko_regulator %in% rownames(TAX$T))) {
    stop("ko_regulator not found in extracted regulators. Available examples: ",
         paste(head(rownames(TAX$T), 10), collapse=", "))
  }


  # ---- lambda estimation (optional) ----
  if (is.null(lambda_A)) {
    lambda_A <- estimate_lambda_stage1(TAX, motif_map, peaks = peaks_use, seed = seed)
  }
  if (is.null(lambda_X)) {
    lambda_X <- estimate_lambda_stage2(TAX, p2g_df = p2g, genes = genes_use, seed = seed)
  }

  # ---- fit models ----
  fit1 <- fit_stage1_T_to_A(TAX, motif_map, peaks = peaks_use, lambda = lambda_A)
  fit2 <- fit_stage2_A_to_X(
    TAX, p2g_df = p2g, genes = genes_use, lambda = lambda_X,
    obj = obj, atac_assay = atac_assay,
    peak_weights = peak_w,
    include_T_direct = include_T_direct,
    motif_map = motif_map
  )

  # ---- predict KO ----
  pred <- predict_virtual_ko(
    TAX, fit1, fit2,
    ko_regulators = ko_regulator,
    ko_value = ko_value,
    ko_mode = ko_mode
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
      lambda_A = lambda_A,
      lambda_X = lambda_X,
      ko_regulator = ko_regulator,
      ko_value = ko_value,
      ko_mode = ko_mode
    ),
    gwas = list(snp_peak = snp_peak, peak_weights = peak_w),
    fits = list(stage1 = fit1, stage2 = fit2),
    diagnostics = list(stage1 = fit1$diagnostics, stage2 = fit2$diagnostics),
    pred = pred,
    rank_genes = rank_genes,
    rank_peaks = rank_peaks
  )
}
