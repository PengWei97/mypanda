# ref-1: https://mooseframework.inl.gov/modules/phase_field/Linearized_Interface_Grain_Growth.html
# 优化算法，提高计算效率，文中考虑了GBEvolution, GBAnisotropy, GBAnisotropyMisori

# my_filename = 'c1_linIf_gg_pfm_GBIso'
# my_filename = 'c2_linIf_gg_pfm_GBAniso'
my_filename = 'c3_linIf_gg_pfm_GBAniso_misori'

[GlobalParams]
  op_num = 2
  var_name_base = aux_gr
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
  [./phi0] # 实际计算变量，替换序参数 aux_gr0
  [../]
  [./phi1]
  [../]
[]

[ICs]
  [./phi0_IC]
    type = SmoothCircleICLinearizedInterface
    variable = phi0
    invalue = 1.0
    outvalue = 0.0
    bound_value = 5.0
    radius = 40
    int_width = 2
    x1 = 50.0
    y1 = 50.0
    profile = TANH
  [../]
  [./phi1_IC]
    type = SmoothCircleICLinearizedInterface
    variable = phi1
    invalue = 0.0
    outvalue = 1.0
    bound_value = 5.0
    radius = 40
    int_width = 2
    x1 = 50.0
    y1 = 50.0
    profile = TANH
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = param2_grn_2_rand_2D_45.tex
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2

    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL
  [../]
[]

[AuxVariables]
  [./aux_gr0]
  [../]
  [./aux_gr1]
  [../]
  [./bounds_dummy]
  [../]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./gr0]
    type = LinearizedInterfaceAux # 基于 gr0_aux = (1.0 + std::tanh(phi0 / std::sqrt(2.0))) / 2.0;
    variable = aux_gr0
    nonlinear_variable = phi0
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./gr1]
    type = LinearizedInterfaceAux
    variable = aux_gr1
    nonlinear_variable = phi1
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    #  phi1, Phi, phi2
    execute_on = 'initial timestep_end'
  [../]
[]

[Kernels]
  #phi0 Kernels
  [./phi0_dot]
    type = ChangedVariableTimeDerivative
    variable = phi0 # 场变量
    order_parameter = gr0 # 序参数，材料参数
  [../]
  [./phi0_ACInt]
    type = ACInterfaceChangedVariable
    variable = phi0
    kappa_name = kappa_op
    mob_name = L
    order_parameter = gr0
  [../]
  [./gr0_AC]
    type = ACGrGrPolyLinearizedInterface
    variable = phi0
    mob_name = L
    this_op = gr0
    other_ops = gr1
    v = phi1
  []
  #phi1 Kernels
  [./phi1_dot]
    type = ChangedVariableTimeDerivative
    variable = phi1
    order_parameter = gr1
  [../]
  [./phi1_ACInt]
    type = ACInterfaceChangedVariable
    variable = phi1
    kappa_name = kappa_op
    mob_name = L
    order_parameter = gr1
  [../]
  [./gr1_AC]
    type = ACGrGrPolyLinearizedInterface
    variable = phi1
    mob_name = L
    this_op = gr1
    other_ops = gr0
    v = phi0
  [../]
[]

[Materials]
  [./gr0]
    type = LinearizedInterfaceFunction # 基于\phi获取序参数，\eta，材料参数
    property_name = gr0
    phi = phi0
    outputs = my_exodus
    output_properties = 'gr0'
  [../]
  [./gr1]
    type = LinearizedInterfaceFunction
    property_name = gr1
    phi = phi1
    outputs = my_exodus
    output_properties = 'gr1' 
  [../]
  # [./GBEovlution]
  #   type = GBEvolution
  #   T = 500 # K
  #   wGB = 5 # nm
  #   GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
  #   Q = 0.23 # Migration energy in eV
  #   GBenergy = 0.708 # GB energy in J/m^2
  #   outputs = my_exodus
  #   output_properties = 'kappa_op gamma_asymm L mu'
  # [../]
  # [./GBAnisotropy]
  #   type = GBAnisotropy
  #   T = 500 # K
  #   wGB = 5

  #   Anisotropic_GB_file_name = param1_anisotropy.tex
  #   inclination_anisotropy = false
  #   # delta_sigma = 0.1
  #   # delta_mob = 0.0

  #   outputs = my_exodus
  #   output_properties = 'kappa_op gamma_asymm L mu'
  # [../]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori
    T = 500 # K
    wGB = 5.0
  
    GBsigma_HAGB = 0.708 # GB energy in J/m^2
    GBmob_HAGB = 2.5e-6

    TT1_sigma = 0.276 # 0.1019 0.3109 0.276
    CT1_sigma = 0.291 # 0.0616 0.1848 0.291
    # TT1_mob = ${my_tt1_mob}
    # CT1_mob = ${my_ct1_mob}

    grain_tracker = grain_tracker
    euler_angle_provider = euler_angle_file

    gb_energy_anisotropy = true
    gb_mobility_anisotropy = true

    output_properties = 'L mu misori_angle'
    outputs = my_exodus
  [../]
[]

[Bounds]
  [./phi0_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi0
    bound_type = upper
    bound_value = 5.0
  [../]
  [./phi0_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi0
    bound_type = lower
    bound_value = -5.0
  [../]
  [./phi1_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi1
    bound_type = upper
    bound_value = 5.0
  []
  [./phi1_lower_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = phi1
    bound_type = lower
    bound_value = -5.0
  [../]
[]

[Postprocessors]
  [./grain_area_mat]
    type = ElementIntegralMaterialProperty
    mat_prop = gr0
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
  scheme = bdf2

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -ksp_type -snes_type'
  petsc_options_value = 'bjacobi gmres vinewtonrsls'

  # solve_type = PJFNK
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_type -snes_type'
  # petsc_options_value = 'hypre boomeramg 31 0.7 boomeramg vinewtonrsls'

  end_time = 10
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