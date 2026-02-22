# multiomeKO 代码逻辑逐功能审查（生物学 + 数学）

## 总体结论
- 该包**确实实现了两阶段虚拟 KO 主流程**：Regulator→ATAC（stage1）与 ATAC→RNA（stage2），再做反事实传播与效应排序。
- 主体思路在生物学上是合理的：共享 metacell、先验约束（motif/LinkPeaks）、深度 offset、可选 GWAS 权重。
- 但从数学与工程稳健性看，仍有若干“可运行但需谨慎解释”的点，尤其是并行随机性、特征跳过比例、训练/预测特征定义差异。

## 逐模块检查

### 1) 主流程 `run_virtual_ko_optimized`
**逻辑正确点**
- 流程顺序完整：Links→GWAS→metacell→TAX→lambda→fit→predict→rank。
- `ko_value/ko_mode/n_cores` 已在入口暴露并向下传递，接口可用性好。

**风险点**
- `genes_use <- intersect(genes_use, unique(p2g$gene))` 会直接过滤掉未在 Links 中出现的基因；若用户关注基因不在 Links，将被静默剔除后仅在长度为 0 时报错。
- `tf_genes` 自动推断依赖 motif 名称解析 TF symbol（`sub(".*_", "", motif)`），在 motif 命名不规范时可能引入错误映射。

### 2) `extract_TAX_metacell`
**逻辑正确点**
- RNA/ATAC 共用 metacell id，符合 paired multiome 假设。
- 使用 count 聚合并提供 `offset = log(libsize+1)`，适合 Poisson GLM。

**风险点**
- 同时使用 `offset` 与 `covariates = nCount_*` 可能带来“深度双重校正”风险（并非必错，但需做灵敏度分析）。

### 3) `estimate_lambda_stage1/2`
**逻辑正确点**
- 用 sampled feature 的 CV 估计全局 lambda，计算量可控。
- `.safe_nfolds()` 处理小样本 CV 边界问题，避免 `cv.glmnet` 崩溃。

**风险点（并行后新增）**
- 并行 CV 的随机折分可能导致跨平台/跨次运行不完全可复现；当前未显式设置并行 RNG stream。
- 使用全局单个 lambda（而非按 peak/gene 局部）可能在异质特征上过/欠收缩。

### 4) `fit_stage1_T_to_A`
**逻辑正确点**
- 候选 regulator 由 motif 命中约束，生物可解释性较强。
- `diagnostics` 记录 `skipped_no_motif`、`skipped_few_regs`，有助于识别模型退化。

**风险点**
- `length(regs) < 3` 直接跳过可能让不少 peak 仅保留截距，若比例过高会降低 KO 传播真实性。

### 5) `fit_stage2_A_to_X`
**逻辑正确点**
- 以 p2g 限定峰集，且 penalty 合并 GWAS+link+distance，符合先验融合方向。
- 支持 `include_T_direct`，可比较“完全中介”与“中介+直接”两种假设。

**风险点**
- `length(pks) < 5` 直接跳过同样会导致部分基因退化为截距主导。
- 并行时每个任务捕获对象较大（`TAX/p2g/obj` 等），PSOCK 在 Windows 下通信开销显著，核心数过高可能反而变慢。

### 6) `predict_virtual_ko`
**逻辑正确点**
- 先 stage1 生成 ATAC 反事实，再进 stage2，符合设定的因果方向。
- `ko_mode = set/scale` 使 KO 语义更灵活。

**关键数学注意点**
- stage2 训练使用的是观测 `A_log1p = log1p(A_count)`，预测时使用 `log1p(muA)`；
  严格来说 `E[log(1+A)] != log(1+E[A])`（Jensen gap），这会引入系统偏差，但在工程上常作为近似。

### 7) 并行实现 `.parallel_lapply`
**逻辑正确点**
- 采用 PSOCK，Linux/Windows 都可运行。
- `n_cores<=1` 自动顺序退化，兼容性好。

**风险点**
- 当前未做 worker 级错误聚合与重试策略；单任务报错会中断整轮。
- 未显式设置并行随机种子流，影响严格复现实验。

## 是否“逻辑正确”与“达成目标”的判断
- **逻辑正确（工程可运行）**：是。
- **生物学结论可直接高置信发布**：建议谨慎，先补齐稳定性验证。

## 建议优先改进（P0）
1. 在并行入口加入可复现 RNG（如 `clusterSetRNGStream`）与 deterministic fold 方案。
2. 对 stage1/stage2 的低拟合覆盖率提供可配置阈值（warning→error）。
3. 增加重复种子/重采样稳定性报告（rank correlation + overlap）。
4. 报告“被过滤基因/峰”的清单，减少 silent drop。

