# MultiomeKO 优化方案（可执行）

## 目标
在保留“两阶段 Regulator→ATAC→RNA”框架下，提升：
1) 稳健性；2) 生物学可解释性；3) 结果可复现性。

## 已完成的首批优化
- **CV 折数安全化**：新增 `.safe_nfolds()`，在样本过少时跳过 CV，避免 `cv.glmnet` 边界失败。
- **默认 TF 回退策略**：当 `tf_mode` 包含 RNA 且 `tf_genes=NULL` 时，自动回退到 `genes_use` 并发出 warning。
- **KO 语义可配置**：`predict_virtual_ko()` 新增 `ko_mode`：
  - `set`：设为常数（原逻辑）
  - `scale`：按比例缩放（更接近“部分抑制”）
- **拟合诊断输出**：stage1/stage2 返回 `diagnostics`，记录拟合与跳过数量，便于评估模型退化。

## 下一步（建议按优先级）
### P0（必须）
1. **增加稳定性评估函数**
   - 对 metacell 重采样（不同 seed）重复 KO，输出 rank 一致性（Spearman/Kendall）。
2. **负对照机制**
   - 增加随机 regulator KO 与低表达 regulator KO，对比真实 KO 的分布分离度。
3. **显式报告退化比例**
   - 若 `fitted_peaks/total_peaks` 或 `fitted_genes/total_genes` 过低，主函数给出 warning。

### P1（重要）
4. **不确定性区间**
   - 对 `dX` 与 `rank_genes` 增加 bootstrap CI。
5. **先验鲁棒映射层**
   - 统一处理 Links/Motif 不同版本字段名，加入 schema 校验。
6. **lambda 策略增强**
   - 支持按基因分层 lambda（而非全局单值）并对比稳定性。

### P2（增强）
7. **生物知识约束增强**
   - stage1 增加 TF family/group 稀疏结构惩罚（可选）。
8. **可视化报告**
   - 自动输出 KO 效应火山图、通路富集与关键 peak-gene-regulator 子网。

## 验证清单
- 至少 2 套公开 multiome 数据复现；
- 正负对照分离显著；
- 不同 `n_cells`、`lambda`、`ko_mode` 下结果一致性可接受；
- 与已知 TF 靶基因方向一致。
