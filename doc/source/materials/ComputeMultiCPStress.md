- **创建时间**: 2024年10月24日
- **创建者**: Peng Wei
- **文件目的**: 本文件旨在总结和记录基于 MOOSE 框架的晶体塑性有限元法（CPFEM）模型的学习过程，涵盖输入文件结构、材料类定义以及 CPFEM 理论的相关知识。
---

# 1. 学习资料
- [Compute Multiple Crystal Plasticity Stress](https://mooseframework.inl.gov/source/materials/crystal_plasticity/ComputeMultipleCrystalPlasticityStress.html)

# 2. Compute Multiple Crystal Plasticity Stress - 笔记
> 处理应力残差的[牛顿迭代](../../theory/Newton-Raphson_method.md) ，并结合 [CrystalPlasticityStressUpdateBase](./CPKalidiniUpdate.md) 继承的多个派生类的本构定律来计算**雅可比矩阵**

# 3. 概述
1. ComputeMultipleCrystalPlasticityStress：旨在实现不同的晶体塑性本构定律，并在模型之间具有一定的模块化程度；
2. 与 CrystalPlasticityStressUpdateBase 的派生类：对于特定的本构，可以将变量（滑移率、滑移阻力 & 内部状态变量）自包含在一个材料类中；

# 4. 晶体塑性控制方程
1. 以 an updated Lagrangian incremental form 来实现晶体滑移和应变增量；
2. 在每个高斯积分点（quadrature point），采用 PK2 stess 来确定是否局部收敛。一旦收敛，等效柯西应力值被计算。之后将所计算的柯西应力用于后面 traditional FEM residual 计算中。
   
## 4.1. Constitutive Equation
1. 参考 [R J Asaro. Crystal Plasticity. Journal of Applied Mechanics, 50(4b):921–934, December 1983.](https://doi.org/10.1115/1.3167205), 其形变梯度 $F$ 可以通过极分解为，
$$
\boldsymbol{F} = \boldsymbol{F}^e \boldsymbol{F}^p
$$

$F^p$ 定义在 intermediate configuration configure, $F^e F^p$ 被定义在 final deformed configuration；

The total plastic velocity gradient,
$$
\boldsymbol{L}^{\boldsymbol{p}}=\dot{\boldsymbol{F}}^p \boldsymbol{F}^{p-1}
$$

采用 backward time integration 来计算跟新的塑性形变梯度，
$$
\boldsymbol{F}_{n+1}^{p-1}=\boldsymbol{F}_n^{p-1}\left(\boldsymbol{I}-\boldsymbol{L}^p \Delta t\right)
$$
塑性形变梯度被用于计算弹性形变梯度，其中弹性拉格朗日应变 $\boldsymbol{E}^e$,
$$
\boldsymbol{E}^e=\frac{1}{2}\left(\boldsymbol{F}^{e \top} \boldsymbol{F}^e-\boldsymbol{I}\right)
$$
之后，PK2 stress被用于计算当前构型下的柯西应力，
$$
\boldsymbol{S}=\operatorname{det}\left(\boldsymbol{F}^e\right) \boldsymbol{F}^{e-1} \boldsymbol{\sigma}\left(\boldsymbol{F}^{e-1}\right)^{\top}
$$

为了考虑多个形变机制，通过将每个晶体塑性本构方程的塑性形变梯度求和的方式获取总的塑性速度梯度，
$$
\boldsymbol{L}^{\boldsymbol{p}}=\sum_{i=1}^m \mathcal{L}_i^p
$$
每个塑性本构中的速率梯度定义如下，是所有滑移系滑移增量的求和，
$$
\mathcal{L}_i^p=\sum_{\alpha=1}^n \dot{\gamma}_i^\alpha \boldsymbol{s}_{i, o}^\alpha \otimes \boldsymbol{m}_{i, o}^\alpha
$$
其中上面所列的 $\boldsymbol{s}_{i, o}^\alpha  \boldsymbol{m}_{i, o}^\alpha$ 是位于 the reference configuratio。而且，塑性滑移的演化必须指定确定的滑移系，
$$
\dot{\gamma}_i^\alpha=\hat{\dot{\gamma}}_i^\alpha\left(\tau_i^\alpha, s_i^\alpha\right)
$$
其中，不同的本构中采用不同的函数形式来构建。分解剪切应力表示为，
$$
\tau^\alpha=\operatorname{det}\left(\boldsymbol{F}^\theta\right)\left(\boldsymbol{F}^{\theta^{\top}} \boldsymbol{S} \boldsymbol{F}^{\theta^{-\top}}\right): \boldsymbol{s}_{i, o}^\alpha \otimes \boldsymbol{m}_{i, o}^\alpha
$$

# 继承关系
ComputeMultipleCrystalPlasticityStress > ComputeFiniteStrainElasticStress > ComputeStressBase, GuaranteeConsumer > ComputeGeneralStressBase > DerivativeMaterialInterface<Material>

1. ComputeGeneralStressBase: 定义材料性质 `_stress`, `_elastic_strain`, `_Jacobian_mult`, 获取材料性质 `_extra_stress`, `_mechanical_strain`;
   1. 函数 `initQpStatefulProperties()` ~ 初始化 `_elastic_strain[_qp]`, `_stress[_qp]`
   2. 函数 `computeQpProperties()`: 调用 `computeQpStress()`来计算 `_stress` & `_Jacobian_mult`
2. ComputeStressBase: 添加是否采用 `use_displaced_mesh` bool 选项
3. ComputeFiniteStrainElasticStress: declareProperty `_rotation_total`, getMaterialProperty `_elasticity_tensor`, `_rotation_total_old`, `_strain_increment`, `_rotation_increment`, `_stress_old`, `_elastic_strain_old`
   1. 函数 `initQpStatefulProperties`: 初始化 `_elastic_strain[_qp]`, `_stress[_qp]` && `_rotation_total`
   2. 函数 `computeQpStress`: ...
4. ComputeMultipleCrystalPlasticityStress: 
   1. 函数 `initialSetup`: 初始设定 `_models` 存储晶体塑性本构，`_eigenstrains` 存储本征应变模型；
   2. 函数 `initQpStatefulProperties`: 初始化 `_plastic_deformation_gradient`, `_pk2`, `_total_lagrangian_strain`, `_updated_rotation`，以及 `_models[i]->initQpStatefulProperties()` - 见[CPKalidiniUpdate.md](./CPKalidiniUpdate.md) & `_eigenstrains[i]->initQpStatefulProperties();`
   3. 函数 `computeQpStress`: 
      1. 初始化晶体塑性本构中的内部变量，包括 `_tau, _flow_direction, _slip_resistance, _slip_increment`
      2. 调用函数 `updateStress(_stress[_qp], _Jacobian_mult[_qp])` ~ 

# 计算逻辑
> ComputeMultipleCrystalPlasticityStress 输入 `_elasticity_tensor, _deformation_gradient, _deformation_gradient_old, _crysrot` 来计算 `_stress, _Jacobian_mult` 等等
1. 函数 `computeQpStress`: 
   1. 初始化晶体塑性本构中的内部变量，包括 `_tau, _flow_direction, _slip_resistance, _slip_increment`
   2. 调用函数 `updateStress(_stress[_qp], _Jacobian_mult[_qp])` ~ 更新柯西应力、雅可比矩阵
2. 函数 `updateStress(cauchy_stress, jacobian_mult)`： 
   1. 初始化设定子步次数和总数、形变梯度、形变梯度增量： `num_substep`  
   2. 调用本构-计算施密特因子，需要输入旋转矩阵 `_crysrot[_qp]`
   3. `do {} while(_convergence_failed)` ~ 如果 `_convergence_failed = false` 结束循环
      1. **go-1**
   4. `postSolveQp()`: 计算柯西应力、雅可比矩阵、总的拉格朗日应变、更新旋转矩阵；
3. -> **go-1** `do {} while(_convergence_failed)`
   1. preSolveQp(): 遍历本构-`setInitialConstitutiveVariableValues`, `_pk2[_qp] = _pk2_old[_qp];` old塑性形变梯度的逆；设定子步的dt，赋予 _substep_dt 给每个晶体塑性本构
      1. `setInitialConstitutiveVariableValues`: `s_{n+1} = s_{n}, _previous_substep_s_{n+1} = s_{n}`
   2. for (unsigned int istep = 0; istep < num_substep; ++istep)
      1. 为每个子步加权 `_delta_deformation_gradient ~ _temporary_deformation_gradient`
      2. 调用 `solveQp();` ~ go-2
      3. `if (_convergence_failed)`, 若 `true`，`substep_iter++; num_substep *= 2; break;` ~ 结束循环，之后继续 do {} while，在这里有跟大的 `num_substep`
   3.  若`false`，即收敛，跳出循环
4. -> **go-2** `solveQp()`:
   1. 循环所有塑性本构：
      1. `setSubstepConstitutiveVariableValues()` ~ `_slip_resistance[_qp] = _previous_substep_slip_resistance;`
      2. `calculateSlipResistance() {}`
   2. _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;
   3. 调用 `solveStateVariables()` ~ 基于本构中的滑移函数来求解内部状态变量应力 **go-3**
   4. 判定是否收敛，`_convergence_failed = true, return`
   5. 遍历所有塑性本构 ~ `updateSubstepConstitutiveVariableValues()` 在substep下更新本构变量
      1. `_previous_substep_slip_resistance = _slip_resistance[_qp];`
   6. `_inverse_plastic_deformation_grad_old = _inverse_plastic_deformation_grad;`
5. -> **go-3** `solveStateVariables()`:
   1. `unsigned int iteration = 0; bool iter_flag = true;`
   2. `do {} while(iter_flag && iteration < _maxiterg)` ~ 当 `iter_flag = false` 且 iteration 次数小于 `_maxiterg`-Maximum number of iterations for internal variable update，之后结束循环 - **go-4**
   3. `if (iteration == _maxiterg)`， 设定 `_convergence_failed = true`;
6. -> **go-4** `do {} while(iter_flag && iteration < _maxiterg)`
   1. 调用函数 `solveStress();` ~ solves for stress, updates plastic deformation gradient. - **go-5**
   2. `if (_convergence_failed)` 若为 `true` 则 return
   3. 遍历所有本构模型
      1. `_slip_resistance_before_update = _slip_resistance[_qp];`
      2. 基于函数 $\Delta g = \left| \Delta \gamma \cdot q^{\alpha \beta} \cdot h^{\beta} \right|$, $g^\alpha=g_o+\Delta \gamma^\alpha q^{\alpha \beta} h_o\left|1-\frac{g^\alpha}{g_{s a t}}\right|^a \operatorname{sign}\left(1-\frac{g^\alpha}{g_{s a t}}\right)$ 来更新 `_slip_resistance_increment`
      3. `if (!_models[i]->updateStateVariables())` ~ 当收敛时，确定当前时间步的状态变量和滑移阻力的值
         1. 遍历滑移系，确定当前的滑移阻力的值，基于滑移阻力增量；如果计算所得的滑移阻力小于0，那么返回false，其他情况返回true ~ `_slip_resistance`
      4. `_models[i]->calculateSlipResistance()` ~ calculateSlipResistance() {} 虚函数，可删除
      5. `if (!_models[i]->areConstitutiveStateVariablesConverged())` ~ 当 `asb(_slip_resistance - _slip_resistance_before_update) / _previous_substep_slip_resistance < _resistance_tol` 返回 true，否则返回 false；
7. -> **go-5**: `solveStress()`
   1. 创建变量 `rnorm, rnorm0, rnorm_prev` ~ 当前残差的L2范数、初始残差的L2范数、上一次迭代中的残差
   2. `calculateResidualAndJacobian()`;
   3. `while (rnorm > _rtol * rnorm0 && rnorm > _abs_tol && iteration < _maxiter)` 多次迭代，牛顿拉弗森迭代算法
      1. 计算应力增量 `dpk2`，获取 `_pk2`
      2. `if (_use_line_search && rnorm > rnorm_prev && !lineSearchUpdate(rnorm_prev, dpk2))`
# 问题

# 想法
