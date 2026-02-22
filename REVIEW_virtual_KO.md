# Multiome 虚拟 KO 两阶段流程代码审阅

## 结论（摘要）
- **实现层面**：代码确实实现了“Regulator -> ATAC -> RNA”的两阶段虚拟 KO 主流程，并且在工程上串联完整（先验提取、拟合、反事实预测、排序输出）。
- **算法逻辑层面**：核心建模思路总体合理（共享 metacell、Poisson+offset、生物先验稀疏化、可选 GWAS 加权）。
- **关键风险**：目前实现中存在若干可能影响生物学解释有效性的点（参数估计稳健性、先验映射鲁棒性、KO 值定义、部分边界条件），因此更准确的结论是：
  - 已达成“可运行的两阶段虚拟 KO 原型”；
  - 但距离“稳健、可复现实验结论”的生物学计算目标仍有改进空间。

## 代码证据要点
1. **端到端两阶段流程已实现**：`run_virtual_ko_optimized()` 依次执行 peak2gene、可选 GWAS 权重、metacell 构建、T/A/X 抽取、stage1 拟合、stage2 拟合、KO 预测与排名。
2. **Stage1（Regulator->ATAC）**：使用 Poisson ridge（`glmnet`，`family='poisson'`, `alpha=0`）并加入 ATAC library-size offset。
3. **Stage2（ATAC->RNA）**：使用 Poisson ridge，输入 `log1p(ATAC)` 特征并加入 RNA offset；支持将 GWAS/link/distance 融合为 penalty.factor。
4. **反事实预测**：将 KO 调控因子置为 0，先预测 ATAC 反事实，再传导到 RNA，输出 `dA/dX`（`log1p(mu_cf)-log1p(mu_wt)`）。

## 主要优点
- 共享 metacell 构建 RNA/ATAC，符合 paired multiome 的结构约束。
- offset 机制可部分校正测序深度差异。
- stage1 通过 motif 约束候选调控因子，stage2 通过 LinkPeaks 约束峰-基因连接，具备生物可解释性。
- penalty 端整合 GWAS PIP + link score + 距离先验，方向正确。

## 主要问题与改进建议
1. **`tf_mode='RNA'` 默认却要求 `tf_genes`，可用性风险高**：默认调用可能直接报错，建议自动检测 TF 列表或把默认改为 `chromvar`/`both` 并提供更友好 fallback。
2. **CV 折数边界条件未充分保护**：`nfolds=min(5,nrow(X))` 在 metacell 很少时可能低于 `cv.glmnet` 可接受下限，建议 `max(3, min(5,nrow(X)))` 并在样本过少时跳过 CV。
3. **stage1/2 特征筛选阈值较硬（如 `<3` 或 `<5` 直接跳过）**：会导致部分 peak/gene 退化为仅 offset 预测，建议记录跳过比例并在结果中显式告警。
4. **motif/links 列名启发式推断较脆弱**：不同版本 Signac 输出字段差异较大，建议更严格校验与标准化适配层。
5. **KO 设为 0 的生物意义依赖尺度**：对 chromvar z-score 可能合理，但对 RNA log1pCPM 并不等价于“完全敲除”，建议支持相对抑制（如乘系数）或 baseline shift。
6. **缺乏不确定性评估**：当前输出是点估计排序，建议增加 bootstrap/metacell 重采样与置信区间。

## 是否达成生物学计算目标
- **部分达成**：
  - 若目标是“构建可解释、可运行的多组学虚拟 KO 推理框架”，目前代码已经达成。
  - 若目标是“用于高可信生物学结论与优先级决策（例如实验级靶点排序）”，目前仍建议补齐稳健性验证（重采样、外部验证、负对照/正对照、敏感性分析）。

## 建议的最小验证清单
- 在至少 2 套公开 multiome 数据上重复：检查 KO 对已知 TF 靶基因方向是否一致。
- 比较不同 metacell 粒度（`n_cells`）与不同 lambda 选择策略下排序稳定性。
- 用随机 KO/非表达 TF 作为负对照，确认假阳性率。
- 在 `include_T_direct=TRUE/FALSE` 下比较结果是否过拟合。
