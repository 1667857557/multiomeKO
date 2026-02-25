# multiomeKO

`multiomeKO` is an R package for **interpretable in-silico knockout (virtual KO)** analysis on paired single-cell multiome data (scRNA + scATAC), designed for Seurat/Signac objects.

It implements a biologically constrained two-stage model:

1. **Regulator → ATAC peaks** (stage 1)
2. **ATAC peaks → gene expression** (stage 2)

Then it applies counterfactual perturbation (KO) to predict downstream changes in accessibility and expression.

---

## Key Features

- **End-to-end virtual KO runner** via `run_virtual_ko_optimized()`.
- **Shared metacell aggregation** for RNA and ATAC to preserve paired structure and improve robustness.
- **Poisson ridge models with library-size offsets** in both stages.
- **Biological priors (stage1 可插拔)**:
  - `stage1_mode = "motif"`: motif-constrained TF→peak (默认，最可解释);
  - `stage1_mode = "chip"`: ChIP/CUT&RUN peak set prior for non-TF regulators;
  - `stage1_mode = "hybrid"`: motif + ChIP prior联合。
- **Signac `Links()`-based peak→gene constraints** for stage2.
- **Optional GWAS fine-mapping integration** (SNP PIP → peak weights).
- **Adaptive stage-2 penalties** combining GWAS, link score, and distance-to-TSS information.
- **Counterfactual prediction modes**:
  - `ko_mode = "set"` (set regulator value to a constant);
  - `ko_mode = "scale"` (scale regulator by a factor).
- **CRISPRi 效应强度 α**:
  - `ko_mode = "scale"` 时可用 `alpha`/`alpha_grid`;
  - `alpha = "auto"` + `alpha_grid` 会返回 `diagnostics$alpha_scan`（当前为启发式校准）。
- **Diagnostics in model outputs** to track fitted/skipped peaks/genes.
- **Cross-platform parallel acceleration** (`n_cores`) using PSOCK clusters (Linux/Windows).

---

## Typical Workflow

### Inputs

A Seurat object with:

- RNA assay (default: `RNA`)
- ATAC assay (default: `ATAC`)
- peak-to-gene links available in `Signac::Links(obj[["ATAC"]])`
- motif data available for ATAC assay (`Signac::GetMotifData`)
- optionally chromVAR assay (if `tf_mode` includes `chromvar`)

### One-call execution

```r
library(multiomeKO)

res <- run_virtual_ko_optimized(
  obj = multiome_obj,
  ko_regulator = "RNA:FOXA1",              # or e.g. "MOTIF:MA0148.1_FOXA1"
  group.by = c("celltype", "SampleID"),
  n_cells = 50,
  genes_use = c("GATA3", "ESR1", "XBP1", "FOXA1"),
  tf_mode = "both",                        # "RNA" | "chromvar" | "both"
  tf_genes = c("FOXA1", "GATA3", "ESR1"),
  include_T_direct = FALSE,
  stage1_mode = "motif",
  chip_peak_map = NULL,
  ko_value = 0,
  alpha = "auto",
  alpha_grid = c(0.25, 0.5, 0.75, 1.0),
  ko_mode = "scale",                       # alpha/alpha_grid are used in scale mode
  finemap_snps = NULL,
  n_cores = 4,
  seed = 1
)

# Top affected genes/peaks
head(res$rank_genes, 20)
head(res$rank_peaks, 20)

# Counterfactual deltas on log1p scale
# rows: features, columns: metacells
dx <- res$pred$dX
# stage diagnostics
res$diagnostics
```

---

## Returned Object (from `run_virtual_ko_optimized`)

A list with major components:

- `setup`: run settings and effective parameters.
- `gwas`: SNP→peak mapping and peak weights (if provided).
- `fits`: fitted stage-1 and stage-2 models.
- `diagnostics`: fit coverage summaries for both stages.
- `pred`: wild-type and counterfactual expected values (`muA`, `muX`) and deltas (`dA`, `dX`).
- `rank_genes`, `rank_peaks`: top effects by mean absolute delta.

---

## Step-by-step API Example

If you prefer explicit control over each step:

```r
library(multiomeKO)

# 1) priors
p2g <- get_peak2gene_links(multiome_obj, atac_assay = "ATAC")
stage1_prior <- build_stage1_priors(multiome_obj, atac_assay = "ATAC", stage1_mode = "motif")
motif_map <- stage1_prior$motif_map

# 2) metacells
mc <- make_metacell_ids(
  multiome_obj,
  group.by = c("celltype", "SampleID"),
  n_cells = 50,
  seed = 1
)

# 3) extract shared T/A/X
TAX <- extract_TAX_metacell(
  obj = multiome_obj,
  metacell_ids = mc,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  genes_use = unique(p2g$gene),
  peaks_use = unique(p2g$peak),
  tf_mode = "both",
  tf_genes = c("FOXA1", "GATA3", "ESR1"),
  covariates = c("nCount_RNA", "nCount_ATAC")
)

# 4) regularization
lambda_A <- estimate_lambda_stage1(TAX, motif_map, peaks = unique(p2g$peak), seed = 1, n_cores = 4)
lambda_X <- estimate_lambda_stage2(TAX, p2g_df = p2g, genes = unique(p2g$gene), seed = 1, n_cores = 4)

# 5) fit models
fit1 <- fit_stage1_T_to_A(TAX, stage1_prior, peaks = unique(p2g$peak), lambda = lambda_A, n_cores = 4)
fit2 <- fit_stage2_A_to_X(
  TAX, p2g_df = p2g, genes = unique(p2g$gene), lambda = lambda_X,
  obj = multiome_obj, atac_assay = "ATAC", n_cores = 4
)

# 6) virtual KO prediction
pred <- predict_virtual_ko(
  TAX, fit1, fit2,
  ko_regulators = "RNA:FOXA1",
  ko_value = 0,
  ko_mode = "set"
)

top_genes <- rank_by_effect(pred$dX, top_n = 50)
```

---

## Notes and Best Practices

- Ensure `Links()` has been computed (`Signac::LinkPeaks`) before running.
- Prefer biologically meaningful `genes_use` for targeted modeling.
- Check `res$diagnostics` to confirm sufficient fitted coverage.
- Compare results across different metacell seeds/sizes for stability.
- 对于非 TF KO，可切换 `stage1_mode = "chip"` 或 `"hybrid"` 并提供 `chip_peak_map`。

---

## Benchmark / Leakage Best Practices

为避免 benchmark 泄漏、提升论文级可复现性，建议：

- `Links()` 必须在 **NTC** 子集上计算，再固定用于 KO/条件评估。
- HVG 与 peaks 筛选尽量在 **NTC** 上完成，避免使用 KO 信息。
- 药物/批次/条件必须 **分层训练与评估**，禁止跨层信息泄漏。
- 报告 `diagnostics`（含 stage1/stage2 覆盖率和 `alpha_scan`）以展示模型退化风险。

---

## Benchmark / Leakage Best Practices

为避免 benchmark 泄漏、提升论文级可复现性，建议：

- `Links()` 必须在 **NTC** 子集上计算，再固定用于 KO/条件评估。
- HVG 与 peaks 筛选尽量在 **NTC** 上完成，避免使用 KO 信息。
- 药物/批次/条件必须 **分层训练与评估**，禁止跨层信息泄漏。
- 报告 `diagnostics`（含 stage1/stage2 覆盖率和 `alpha_scan`）以展示模型退化风险。

---

## Dependencies

Declared imports include `Matrix`, `glmnet`, `Seurat`, `Signac`, `GenomicRanges`, `IRanges`, `GenomeInfoDb`.
