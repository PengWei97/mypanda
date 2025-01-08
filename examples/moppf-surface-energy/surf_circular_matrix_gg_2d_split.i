# This model was created to replicate the results of reference-1;
# reference-1: Verma, M., & Mukherjee, R. (2021). Grain growth stagnation in solid state thin films: A phase-field study. Journal of Applied Physics, 130(2).
# Data: 2025-01-08

my_filename = 'case4_circular_gg_2d_spilt'

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 10
  xmin = 0
  xmax = 200
  ymin = 0
  ymax = 100
  elem_type = QUAD4 # HEX8

  parallel_type = distributed
[]

[Variables]
  [./PolycrystalVariables]
    var_name_base = gr
    op_num = 2
  [../]
  [./cv]
    # order = THIRD
    # family = HERMITE

    order = FIRST
    family = LAGRANGE
  [../]
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Functions]
  [./function_gr0]
    type = ParsedFunction
    expression = 'if((x<50 | x > 150) & (y < 80), 1.0, 0.0)'
  [../]
[]

[ICs]
  [./gr0]
    type = FunctionIC
    variable = gr0
    function = function_gr0
  [../]
  [./gr1]
    type = BoundingBoxIC
    variable = gr1
    x1 = 50
    y1 = 0
    x2 = 150
    y2 = 80
    inside = 1.0
    outside = 0.0
  [../]
  [./cv]
    type = BoundingBoxIC
    variable = cv
    x1 = 0
    y1 = 80
    x2 = 200
    y2 = 100

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
    type = CoupledTimeDerivative
    variable = w
    v = cv
  [../]
  [./cv_bulk]
    type = SplitCHParsed
    variable = cv
    f_name = f_total
    kappa_name = kappa_c
    w = w
    coupled_variables = 'gr0 gr1'
  [../]
  [./cv_int]
    type = SplitCHWRes
    variable = w
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
    constant_names = 'beta_s beta_gb A0'
    constant_expressions = '1.0 1.0 1.0'
    expression = 'A0*(cv^4/4.0 - cv^2/2.0 + gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + beta_s*cv^2*(gr0^2 + gr1^2) + beta_gb*gr0^2*gr1^2)'
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

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type  -snes_type -ksp_gmres_restart'
  petsc_options_value = 'bjacobi vinewtonrsls 31'

  # scheme = bdf2
  # solve_type = NEWTON
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  # petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 15
  l_tol = 1e-4
  nl_max_its = 10
  nl_rel_tol = 1e-9

  end_time = 1.0e4
  # num_steps = 3
  # dt = 0.05
  # dtmax = 0.2

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
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
    # time_step_interval = 10 # The interval at which time steps are output
    # sync_times = '10 50 100 500 1000 5000 10000 50000 100000'
    # sync_only = true
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  print_linear_residuals = false
[]