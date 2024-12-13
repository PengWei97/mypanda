# 环晶：晶粒长大benchmark算例

my_filename = 'c1_gg_pfm_GBIso'

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 100
  ymax = 100
  nx = 50
  ny = 50
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC]
      radius = 40
      x = 50.0
      y = 50.0
      int_width = 2
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[Materials]
  [./GBEovlution]
    type = GBEvolution
    GBenergy = 0.97
    GBMobility = 0.6e-6
    T = 300
    wGB = 5
  [../]
[]

[Postprocessors]
  [./grain_area_mat]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
    execute_on = 'initial TIMESTEP_END'
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2 # 使用二阶后向差分公式

  solve_type = NEWTON # 使用牛顿法求解非线性方程
  petsc_options_iname = '-pc_type -ksp_type -snes_type' # PETSc 求解器选项名称
  petsc_options_value = 'bjacobi gmres vinewtonrsls' # PETSc 求解器配置

  # solve_type = PJFNK
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 31 0.7'

  # l_tol = 1.0e-4 # 线性求解器容许误差
  # l_max_its = 30
  # nl_max_its = 25
  # nl_rel_tol = 1.0e-7

  end_time = 1.3
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2 # 8 
    # cycles_per_step = 2 # The number of adaptivity cycles per step
    refine_fraction = 0.5 # The fraction of elements or error to refine.
    coarsen_fraction = 0.05
    max_h_level = 2
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/o_${my_filename}
  [./my_exodus]
    type = Exodus
    time_step_interval = 5
  [../]
  csv = true
[]