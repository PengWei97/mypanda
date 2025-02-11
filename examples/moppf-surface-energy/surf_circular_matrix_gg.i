my_filename = 'case1_circular_gg'
my_wGB = 2.0
# my_interval = 1.0

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 50
  ny = 50
  nz = 10
  xmin = 0
  xmax = 500
  ymin = 0
  ymax = 500
  zmin = 0
  zmax = 100
  elem_type = HEX8

  parallel_type = distributed
[]

[Variables]
  # [./gr0]
  # [../]
  # [./gr1]
  # [../]
  [./PolycrystalVariables]
    var_name_base = gr
    op_num = 2
  [../]
  # [./c_v]
  # [../]
[]

[ICs]
  [./gr0]
    type = SmoothCircleIC
    variable = gr0
    x1 = 256
    y1 = 256
    z1 = 56

    radius = 200
    invalue = 0.0
    outvalue = 1.0
    z_threshold = 48
    int_width = ${my_wGB}
    3D_spheres = false
    zero_gradient = false
  [../]
  [./gr1]
    type = SmoothCircleIC
    variable = gr1
    x1 = 256
    y1 = 256
    z1 = 56

    radius = 200
    invalue = 1.0
    outvalue = 0.0
    z_threshold = 48
    int_width = ${my_wGB}
    3D_spheres = false
    zero_gradient = false
  [../]
[]

[GlobalParams]
  var_name_base = gr
  op_num = 2
[]
[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = FIRST
    family = LAGRANGE
  [../]
  [./var_indices]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./gr0dot] # 
    type = TimeDerivative
    variable = gr0
  [../]
  [./gr0bulk]
    type = AllenCahn
    variable = gr0
    f_name = f_total
    coupled_variables = gr1
    mob_name = mob_AC
  [../]
  [./gr0int]
    type = ACInterface
    variable = gr0
    kappa_name = kappa_op
    mob_name = mob_AC
  [../]

  [./gr1dot]
    type = TimeDerivative
    variable = gr1
  [../]
  [./gr1bulk]
    type = AllenCahn
    variable = gr1
    f_name = f_total
    coupled_variables = gr0
    mob_name = mob_AC
  [../]
  [./gr1int]
    type = ACInterface
    variable = gr1

    kappa_name = kappa_op
    mob_name = mob_AC
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./free_energy_solid]
    type = DerivativeParsedMaterial
    property_name = f_solid
    coupled_variables = 'gr0 gr1'
    constant_names = 'mu gamma_asymm'
    constant_expressions = '1.0 1.5'
    expression = 'mu*(gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2)' #  + 1.0/4.0
    derivative_order = 2
    enable_jit = true
  [../]
  [./free_energy_total]
    type = ParsedMaterial
    property_name = f_total
    material_property_names = 'f_solid'
    expression = 'f_solid'
  [../]
  [./mob_AC]
    type = GenericConstantMaterial
    prop_names = 'mob_AC kappa_op' # L kappa_op
    prop_values = '2.0 2.0'
  [../]
[]

[Preconditioning] # 预处理器，用于加速求解线性系统的收敛
  [./SMP]
    type = SMP # 表示使用对称多处理（SMP）并行算法。该预处理器设计用于多核处理器上，以优化内存访问和并行计算，提高求解效率。
    full = true
  [../]
[]

[Executioner]
  type = Transient # 使用瞬态求解器
  # BDF2 (Backward Differentiation Formula) 时间积分方案
  scheme = bdf2 # 这是一个二阶精度的显式时间步进方法，通常用于解决瞬态问题，特别适用于硬度较高的系统。

  solve_type = NEWTON # 使用牛顿法求解非线性方程
  petsc_options_iname = '-pc_type -ksp_type -snes_type' # PETSc 求解器选项名称
  petsc_options_value = 'bjacobi gmres vinewtonrsls' # PETSc 求解器配置

  # # Preconditioned JFNK (default)
  # solve_type = 'PJFNK' # JFNK (Preconditioned Jacobian-Free Newton Krylov) 求解方法

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  # petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 20
  l_tol = 1e-4
  nl_max_its = 10
  nl_rel_tol = 1e-9

  # end_time = 1e5
  # num_steps = 10
  # dt = 0.05
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 2.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2 # 8 
    cycles_per_step = 2 # The number of adaptivity cycles per step
    refine_fraction = 0.5 # The fraction of elements or error to refine.
    coarsen_fraction = 0.05
    max_h_level = 2
  [../]
[]

[Outputs]
  [./my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    type = Nemesis
    # time_step_interval = ${my_interval} # The interval at which time steps are output
    # sync_times = '10 50 100 500 1000 5000 10000 50000 100000'
    # sync_only = true
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]    
[]