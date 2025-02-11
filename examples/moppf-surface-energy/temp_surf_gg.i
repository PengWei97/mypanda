# Example: Simulation of grain growth in materials exhibiting GBAnisotropy material class;
# The initial condition includes a circular grain embedded in a matirx griains, with the aim to study the influence of GB anisotropy on the grain growth procress.
# TODO - 不采用内置的GBEvolution来处理，直接采用材料类+kernels的形式

my_filename = 'case1_circular_gg'
my_wGB = 5.0
my_number_adaptivity = 3
my_interval = 2

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 5
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  zmax = 50
  elem_type = HEX8 # QUAD4

  # uniform_refine = 1
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  wGB = ${my_wGB}
  length_scale = 1.0 # 1.0e-9
  time_scale = 1.0 # 1.0e-9
[]

[UserObjects]
  [./term]
    type = Terminator
    expression = 'gr0_area < 1000'
  [../]
[]

[Variables]
  [./PolycrystalVariables]
  [../]
  # [./c_v] # vapor phase (c_v = 1), and solid film (c_v = 0)
  # [../]
[]

[ICs]
  # [./PolycrystalICs]
  #   [./BicrystalCircleGrainIC] # 3D cylinders
  #     # BicrystalCircleGrainICAction -> SmoothCircleIC -> SmoothCircleBaseIC
  #     radius = 35
  #     x = 50
  #     y = 50
  #     3D_sphere = false
  #     int_width = ${my_wGB}
  #   [../]
  # [../]
  [./vapor_phase_condition]
    type = SmoothCircleIC
    variable = gr0
    x1 = 50
    y1 = 50
    z1 = 0.0
    radius = 35
    
    invalue = 1.0
    outvalue = 0.0
    z_threshold = 150
    int_width = ${my_wGB}
    3D_spheres = false
    zero_gradient = false
  [../]
  [./BicrystalCircleGrainIC_gr1]
    type = SmoothCircleIC
    variable = gr1
    x1 = 50
    y1 = 50
    z1 = 0.0
    radius = 35
    
    invalue = 0
    outvalue = 1.0
    z_threshold = 150
    int_width = ${my_wGB}
    3D_spheres = false
    zero_gradient = false
  [../]
  # [./c_v]
  #   type = FunctionIC
  #   variable = c_v
  #   function = fn
  # [../]
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
  # [./c_v_dot]
  #   type = TimeDerivative
  #   variable = c_v
  # [../]
  # [./SHVapor]
  #   type = CahnHilliard
  #   variable = c_v
  #   mob_name = mob_surf
  #   f_name = fvapor
  # [../]
  # [./SHInterface]
  #   type = CahnHilliard
  #   variable = c_v
  #   mob_name = mob_surf
  #   f_name = fintf
  # [../]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

# [BCs]
#   [./Periodic]
#     [./All]
#       auto_direction = 'x y z'
#       variable = 'gr0 gr1'
#     [../]
#   [../]
# []

[Materials]
  # [./CuGrGranisotropic]
  #   type = GBEvolutionBenchmark # gamma ~ beta_gb
    
  #   T = 600 # K

  #   GBenergy = 1.0
  #   GBMobility = 1.0;

  #   output_properties = 'L mu gamma_asymm kappa_op' 
  #   outputs = my_exodus
  # [../]
  [./free_energy_solid_film]
    type = DerivativeParsedMaterial
    property_name = f_solid
    material_property_names = 'A L '
  [../]
  # [./free_energy_of_vapor]
  #   type = DerivativeParsedMaterial
  #   property_name = fvapor
  #   material_property_names = 'mu'
  #   expression = 'mu*(1/4*c_v^4-1/2*c_v^2)'
  #   coupled_variables = c_v
  #   derivative_order = 2
  # [../]
  # [./free_energy_of_interface]
  #   type = DerivativeParsedMaterial
  #   property_name = fintf
  #   coupled_variables = 'c_v gr0 gr1'
  #   material_property_names = 'mu' # mu ~ A0
  #   constant_names = 'beta_s'
  #   constant_expressions = '1.0'
  #   expression = 'mu*(beta_s*c_v^2*gr0^2+beta_s*c_v^2*gr1^2)'
  # [../]
  # [./surface_mobility]
  #   type = DerivativeParsedMaterial
  #   property_name = mob_surf
  #   coupled_variables = 'c_v'
  #   constant_names = 'Mb Ms'
  #   constant_expressions = '1.0e-4 10e-2'
  #   expression = 'Mb+16*Ms*(1-c_v)^2*c_v^2'
  # [../]
[]

[Postprocessors]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]

  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
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

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK' # JFNK (Preconditioned Jacobian-Free Newton Krylov) 求解方法

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 20
  l_tol = 1e-4
  nl_max_its = 10
  nl_rel_tol = 1e-9

  # end_time = 50
  num_steps = 10

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 2.5e-2
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3 # 8 
    cycles_per_step = 2 # The number of adaptivity cycles per step
    refine_fraction = 0.5 # The fraction of elements or error to refine.
    coarsen_fraction = 0.05
    max_h_level = ${my_number_adaptivity}
  [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    type = Checkpoint
    num_files = 6
    time_step_interval = 2
  [../] 
  [./my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    type = Nemesis
    time_step_interval = ${my_interval} # The interval at which time steps are output
    # sync_times = '10 50 100 500 1000 5000 10000 50000 100000'
    # sync_only = true
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  # [./pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial timestep_end final'  # Default is "final"
  #   level = 2                     # Default is 1
  #   heaviest_branch = true        # Default is false
  #   heaviest_sections = 2         # Default is 0
  # [../]
  [./my_console]
    type = Console
    output_linear = false
    # output_screen = false
    # time_step_interval = 5
  [../]
[]
