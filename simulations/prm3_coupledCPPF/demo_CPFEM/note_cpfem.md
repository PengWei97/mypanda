- **创建时间**: 2024年10月16日
- **创建者**: Peng Wei
- **文件目的**: 本文件旨在总结和记录基于 MOOSE 框架的晶体塑性有限元法（CPFEM）模型的学习过程，涵盖输入文件结构、材料类定义以及 CPFEM 理论的相关知识。
- TODO
---

## 内容概览
1. **Input files**: 讨论 MOOSE 中输入文件的构成与格式。
2. **Material class**: 介绍材料类的定义与实现。
3. **Theory of CPFEM**: 简述晶体塑性有限元法的理论基础。

# Input files
## 学习模板
- exception.i

### Material class
```bash
[Materials]
  [elasticity_tensor]
    type = ComputeElasticTenosrCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  []
  [stress]
    type = ComputeMultiCPStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    maximum_substep_iteration = 1
  []
  [trial_xtalpl]
    type = CPKalidiniUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.inp
  []
[]
```

- **elasticity_tensor**: 调用 `ComputeElasticTenosrCP` 来为 CP 弹性模量，TODO-阅读以及简化继承
- **stress**: 调用 `ComputeElasticTenosrCP` 计算 CP, 输入：；输出：
- **trial_xtalpl**: 结合 `ComputeMultiCPStress` 来提供 CP 本构，该本构采用经典的幂率定律，可以替换为基于位错密度的 CP 本构；

### Kernels
```bash
[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  add_variables = true
  generate_output = stress_zz
[]
```
- 采用 `action` 的形式，如下：
```c++
registerSyntax("QuasiStaticSolidMechanicsPhysics", "Physics/SolidMechanics/QuasiStatic/*");
```
用来调用 `QuasiStaticSolidMechanicsPhysics` action ~ 准静态固体力学模块。[Solid Mechanics QuasiStatic Physics System](https://mooseframework.inl.gov/syntax/Physics/SolidMechanics/QuasiStatic/index.html) - TODO

# Theory of CPFEM

# Refs
- [Compute Multiple Crystal Plasticity Stress](https://mooseframework.inl.gov/source/materials/crystal_plasticity/ComputeMultipleCrystalPlasticityStress.html)
- [CrystalPlasticityKalidindiUpdate](https://mooseframework.inl.gov/source/materials/crystal_plasticity/CrystalPlasticityKalidindiUpdate.html#)
- 