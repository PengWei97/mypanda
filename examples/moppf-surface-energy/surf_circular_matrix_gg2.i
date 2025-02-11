my_filename = 'case2_circular_gg'

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  zmax = 100
  elem_type = HEX8

  parallel_type = distributed
[]

[Variables]
  [./PolycrystalVariables]
    var_name_base = gr
    op_num = 2
  [../]
  [./cv]
  [../]
[]

[Bounds]
  [./c_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cv
    bound_type = upper
    bound_value = 1.0
  [../]
  [./c_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = cv
    bound_type = lower
    bound_value = 0.0
  [../]
[]


[ICs]
  [./gr0]
    type = SmoothCircleIC
    variable = gr0
    x1 = 50
    y1 = 50
    z1 = 50

    radius = 30
    invalue = 0.0
    outvalue = 1.0
    z_threshold = 80
    # int_width = ${my_wGB}
    3D_spheres = false
    zero_gradient = false
  [../]
  [./gr1]
    type = SmoothCircleIC
    variable = gr1
    x1 = 50
    y1 = 50
    z1 = 50

    radius = 30
    invalue = 1.0
    outvalue = 0.0
    z_threshold = 80
    # int_width = ${my_wGB}
    3D_spheres = false
    zero_gradient = false
  [../]
  [./cv]
    type = BoundingBoxIC
    variable = cv
    x1 = 0
    y1 = 0
    z1 = 80 
    x2 = 100
    y2 = 100
    z2 = 100

    inside = 1.0
    outside = 0.0
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
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./gr0dot] # dgr0/dt
    type = TimeDerivative
    variable = gr0
  [../]
  [./gr0bulk] # df_loc/dgr0
    type = AllenCahn
    variable = gr0
    f_name = f_total
    coupled_variables = 'gr1 cv'
    mob_name = L_phi
  [../]
  [./gr0int] # df_gr/dgr0
    type = ACInterface
    variable = gr0
    kappa_name = kappa_phi
    mob_name = L_phi
  [../]

  [./gr1dot]
    type = TimeDerivative
    variable = gr1
  [../]
  [./gr1bulk]
    type = AllenCahn
    variable = gr1
    f_name = f_total
    coupled_variables = 'gr0 cv'
    mob_name = L_phi
  [../]
  [./gr1int]
    type = ACInterface
    variable = gr1
    kappa_name = kappa_phi
    mob_name = L_phi
  [../]

  [./cvdot]
    type = TimeDerivative
    variable = cv
  [../]
  [./cv_bulk]
    type = CahnHilliard
    variable = cv
    f_name = f_total
    coupled_variables = 'gr0 gr1'
    mob_name = M
  [../]
  [./cv_int]
    type = CHInterface
    variable = cv
    kappa_name = kappa_c
    mob_name = M
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
    property_name = f_total
    coupled_variables = 'gr0 gr1 cv'
    constant_names = 'beta_s beta_gb'
    constant_expressions = '1.0 1.0'
    expression = '1.0*(cv^4/4.0 - cv^2/2.0 + gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + beta_s*cv^2*(gr0^2 + gr1^2) + beta_gb*gr0^2*gr1^2)'
    derivative_order = 2
    enable_jit = true

    output_properties = 'f_total df_total/dgr0 df_total/dgr1 df_total/dcv'
    outputs = my_exodus
  [../]
  [./surface_mobility]
    type = ParsedMaterial
    property_name = M
    constant_names = 'Mb Ms'
    constant_expressions = '10e-4 10'
    coupled_variables = 'cv'
    expression = 'Mb + 16*Ms*(1-cv)^2*cv^2'

    output_properties = 'M'
    outputs = my_exodus
  [../]
  [./const_materials]
    type = GenericConstantMaterial
    prop_names = 'kappa_c kappa_phi L_phi'
    prop_values = '2 2 2'
  [../]
[]

[Postprocessors]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
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
  # scheme = bdf2 # 这是一个二阶精度的显式时间步进方法，通常用于解决瞬态问题，特别适用于硬度较高的系统。

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK' # JFNK (Preconditioned Jacobian-Free Newton Krylov) 求解方法

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  # solve_type = NEWTON # 使用牛顿法求解非线性方程
  # petsc_options_iname = '-pc_type -ksp_type -snes_type' # PETSc 求解器选项名称
  # petsc_options_value = 'bjacobi gmres vinewtonrsls' # PETSc 求解器配置

  l_max_its = 15
  l_tol = 1e-4
  nl_max_its = 10
  nl_rel_tol = 1e-9

  end_time = 20
  # num_steps = 3
  # dt = 0.05
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 2.5e-3
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3 # 8 
    cycles_per_step = 2 # The number of adaptivity cycles per step
    refine_fraction = 0.5 # The fraction of elements or error to refine.
    coarsen_fraction = 0.05
    max_h_level = 4
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