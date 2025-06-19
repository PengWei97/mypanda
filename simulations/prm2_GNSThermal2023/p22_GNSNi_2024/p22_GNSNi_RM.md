# GNSNi 相场模拟备忘录（README）

本备忘文档用于记录 “梯度纳米结构纯镍（GNSNi）” 相场模拟过程中的关键思路、建模步骤、问题总结、Benchmark 验证、定量化指标提取策略，以及后续开展大尺度模拟的方案。  

---

## 一、模拟构思

### 1.3 模拟目标
- 验证实验结果；
- 剥离关键因素；

---

## 二、建模与实现

### 2.1 模型类型
- 相场模型类型：多晶晶粒生长模型（Allen-Cahn / multi-order parameter）；
- 考虑因素：
  - 晶界各向异性（GBAnisotropy）；
  - 晶粒尺寸梯度；
  - 晶界能与迁移率的空间分布；
  - 可选：热力学驱动力 / 形核 / 位错密度耦合。

### 2.2 初始结构生成
- 使用 Voronoi 方法生成初始晶粒；
- 结合 Python/Neper 脚本，设定尺寸梯度（如Z方向粗到细）；
- 相应 order parameters 分配策略。

### 2.3 参数设置
- 网格尺寸 / 时间步长；
- 相场参数：mobility、界面能、界面宽度等；
- 不同区域材料参数分布（使用 `bnds` or `mask` 设定）。

---

## 三、问题与解决方案

| 问题描述 | 出现场景 | 初步原因 | 解决方案 |
|----------|-----------|----------|-----------|
| EBSD数据计算缓慢 | 计算时 | 网格设定不合理 | MTEX & 设定 `pre_refine` 来粗化, Ref1|

---

## 四、Benchmark 验证

### 4.1 与经典模型对比
- 与等距晶粒、无梯度结构演化对比；
- 与已有文献（Verma & Mukherjee 2021）结果复现。

### 4.2 数值稳定性验证
- 网格数变化测试；
- 时间步长敏感性分析；
- 不同 GB 模型对比测试。

---

## 五、定量分析方法

### 5.1 晶粒尺寸统计
- 使用 Python / ImageJ / MTEX 分析晶粒尺寸分布；
- 分区域提取（纳米区 vs 微米区）；
- 演化曲线提取（平均尺寸 vs 时间）。

### 5.2 晶界速度分析
- 基于追踪晶界位置变化提取速度场；
- 结合粒子追踪算法或边界提取 + 差分。

### 5.3 区域演化对比
- 横截面提取固定区域尺寸；
- 不同梯度区的演化对比图像或曲线。

---

## 六、大尺度模拟计划

### 6.1 目标与挑战
- 目标：近似真实尺寸的模拟（>5 µm，百万节点）；
- 挑战：计算资源 / 结构构建复杂度 / 多尺度参数设定。

### 6.2 实施方案
- 使用 MPI + 并行框架（如 MOOSE）；
- 建立 coarse-grained 模型或降低 order parameter 数量；
- 区域细化策略（仅纳米区使用细网格）。

---

## 七、参考文献

1. 
2. Verma, M., & Mukherjee, R. (2021). *Grain growth stagnation in solid state thin films: A phase-field study.* J. Appl. Phys., **130(2)**.
3. Chen, L.Q. (2002). *Phase-field models for microstructure evolution*. Annual Review of Materials Research.
4. Moelans, N., et al. (2008). *An introduction to phase-field modeling of microstructure evolution*. CALPHAD.
5. MOOSE framework documentation: https://mooseframework.inl.gov/
6. Neper: https://neper.info

---
