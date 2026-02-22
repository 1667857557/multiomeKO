# multiomeKO

`multiomeKO` is an R package for **interpretable two-stage virtual knock-out (KO)** analysis on paired single-cell multiome data (scRNA + scATAC), using Seurat/Signac objects.

## What this package does

The package models regulatory effects in two linked stages:

1. **Stage 1: Regulator -> ATAC peak accessibility**
   - Fits Poisson ridge models from regulator features (`T`) to peak counts (`A`).
   - Uses motif-based candidate restriction for biological interpretability.

2. **Stage 2: ATAC peaks -> Gene expression**
   - Fits Poisson ridge models from peak features to RNA counts (`X`).
   - Uses peak-to-gene links (Signac `Links`) and optional adaptive penalties from GWAS fine-mapping, link score, and distance-to-TSS priors.

Then the package performs **virtual KO counterfactual prediction**:
- Perturb one regulator (e.g. set to zero or scale down),
- Propagate effect through Stage 1 and Stage 2,
- Return predicted changes in accessibility (`dA`) and expression (`dX`),
- Rank peaks/genes by effect size.

## Main features

- Shared metacell construction across RNA and ATAC (paired structure preserved)
- Library-size offsets for count-based models
- Biologically constrained sparsity (motif prior, peak-to-gene prior)
- Optional GWAS SNP -> peak weighting
- End-to-end runner: `run_virtual_ko_optimized()`
- Feature effect ranking for downstream interpretation

## Core workflow

Typical high-level pipeline:

1. Build metacells with `make_metacell_ids()`
2. Extract matrices with `extract_TAX_metacell()`
3. Get priors with `get_peak2gene_links()` and `get_peak_motif_map()`
4. Estimate penalties (`estimate_lambda_stage1/2()`)
5. Fit two-stage models (`fit_stage1_T_to_A()`, `fit_stage2_A_to_X()`)
6. Run KO prediction (`predict_virtual_ko()`)
7. Rank impacts (`rank_by_effect()`)

Or use the one-shot interface:
- `run_virtual_ko_optimized()`

## Minimal example

> Note: this assumes your Seurat object already contains:
> - RNA assay (`RNA`),
> - ATAC assay (`ATAC`),
> - Signac peak-to-gene links (`Links`),
> - optionally chromVAR assay if using `tf_mode = "chromvar"` or `"both"`.

```r
library(multiomeKO)

# obj: Seurat multiome object (RNA + ATAC)
# genes_use: genes of interest for modeling
# ko_regulator: regulator name in TAX$T row naming convention
#   examples: "RNA:FOXA1" or "MOTIF:MA0148.1_FOXA1"

res <- run_virtual_ko_optimized(
  obj = obj,
  ko_regulator = "RNA:FOXA1",
  group.by = c("celltype", "SampleID"),
  n_cells = 50,
  genes_use = c("FOXA1", "GATA3", "ESR1", "KRT8", "KRT18"),
  tf_mode = "both",
  tf_genes = c("FOXA1", "GATA3", "ESR1"),
  include_T_direct = FALSE,
  ko_value = 0,
  ko_mode = "set",   # or "scale"
  seed = 1,
  rna_assay = "RNA",
  atac_assay = "ATAC"
)

# Top predicted affected genes/peaks
head(res$rank_genes, 20)
head(res$rank_peaks, 20)

# Counterfactual delta matrices (features x metacells)
# Gene expression changes:
dX <- res$pred$dX
# Accessibility changes:
dA <- res$pred$dA
```

## Optional GWAS prior example

```r
# finemap_snps must provide chr, pos, pip
# can be path, data.frame, or GRanges
res_gwas <- run_virtual_ko_optimized(
  obj = obj,
  ko_regulator = "RNA:FOXA1",
  group.by = c("celltype", "SampleID"),
  n_cells = 50,
  genes_use = genes_use,
  tf_mode = "both",
  tf_genes = tf_genes,
  finemap_snps = snp_df,
  pip_floor = 0.01
)
```

## Output structure (from `run_virtual_ko_optimized`)

- `setup`: run configuration
- `gwas`: SNP-to-peak mapping and peak weights (if provided)
- `fits`: stage1/stage2 model objects
- `diagnostics`: fit coverage diagnostics for stage1/stage2
- `pred`: wild-type and KO counterfactual predictions (`muA`, `muX`, `dA`, `dX`)
- `rank_genes`, `rank_peaks`: effect-based rankings

## Function index

- Data prep: `make_metacell_ids()`, `extract_TAX_metacell()`
- Priors: `get_peak2gene_links()`, `get_peak_motif_map()`
- GWAS: `read_finemap_snps()`, `snps_to_granges()`, `map_snps_to_peaks()`, `peak_weights_from_snps()`
- Penalty: `build_stage2_peak_penalty()`, `penalty_from_peak_weights()`
- Model fit: `estimate_lambda_stage1()`, `estimate_lambda_stage2()`, `fit_stage1_T_to_A()`, `fit_stage2_A_to_X()`
- Counterfactual: `predict_virtual_ko()`, `rank_by_effect()`, `run_virtual_ko_optimized()`
