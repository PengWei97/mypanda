my_filename = "s1_fullCoupledCPPF"

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 5
  nz = 0
  xmax = 300
  ymax = 50
  zmax = 0
  elem_type = QUAD4

  parallel_type = distributed
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr

  length_scale = 1.0e-9
  time_scale = 1.0e-9

  displacements = 'disp_x disp_y'
[]

[Variables]
  [./PolycrystalVariables]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Physics/SolidMechanics/QuasiStatic/all] # ComputeFiniteStrain + StressDivergenceTensors
  strain = FINITE 

  volumetric_locking_correction = true  
  add_variables = true
  generate_output = stress_yy
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_2_rand_2D.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerMateProp # adding
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
    fill_method = symmetric9
    euler_angle_provider = euler_angle_file
  [../]
[]

[ICs]
  [./PolycrystalICs]
    # [./BicrystalCircleGrainIC]
    #   radius = 150
    #   x = 200
    #   y = 200
    #   int_width = 30
    # [../]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = 150
      y2 = 50
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./pk2]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./gss]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./slip_increment]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./slip_shear_strain_gr0]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./AC_ElasticDrivingForce_gr0]
    type = ACGrGrDeformationEnergyDrivingFC
    variable = gr0
    D_deformation_energy_name = 'delastic_energy/dgr0'
  [../]
  [./AC_ElasticDrivingForce_gr1]
    type = ACGrGrDeformationEnergyDrivingFC
    variable = gr1
    D_deformation_energy_name = 'delastic_energy/dgr1'
  [../]
  [./AC_PlasticDrivingForce_gr0]
    type = ACGrGrDeformationEnergyDrivingFC
    variable = gr0
    D_deformation_energy_name = 'dplastic_energy/dgr0'
  [../]
  [./AC_PlasticDrivingForce_gr1]
    type = ACGrGrDeformationEnergyDrivingFC
    variable = gr1
    D_deformation_energy_name = 'dplastic_energy/dgr1'
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./vonmises_stress]
    type = RankTwoScalarAux
    variable = vonmises_stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
  # [./pk2]
  #  type = RankTwoAux
  #  variable = pk2
  #  rank_two_tensor = second_piola_kirchhoff_stress
  #  index_j = 1
  #  index_i = 1
  #  execute_on = timestep_end
  # [../]
  # [./gss]
  #   type = MaterialStdVectorAux
  #   variable = gss
  #   property = slip_resistance
  #   index = 0
  #   execute_on = timestep_end
  # [../]
  # [./slip_inc]
  #   type = MaterialStdVectorAux
  #   variable = slip_increment
  #   property = slip_increment_gr0
  #   index = 0
  #   execute_on = timestep_end
  # [../]
  # [./slip_shear_strain_gr0]
  #   type = MaterialStdVectorAux
  #   variable = slip_shear_strain_gr0
  #   property = slip_shear_strain_gr0
  #   index = 0
  #   execute_on = timestep_end
  # [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = total_lagrangian_strain # elastic_strain # total_lagrangian_strain
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = -50*0.01*t
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 30 # 5 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2 
  [../]
  [./ElasticityTensor]
    type = CompPolyElastTensorCPCoupled # adding
    block = 0
    grain_tracker = grain_tracker
  [../]    
  [./stress]
    type = CompMultiCPStressFullCoupled 
    crystal_plasticity_models = 'trial_xtalpls'
    tan_mod_type = exact
    maximum_substep_iteration = 100
    min_line_search_step_size = 0.01

    grain_tracker = grain_tracker

    output_properties = 'elastic_energy_gr0 elastic_energy_gr1'
    outputs = my_exodus
  [../]
  [./trial_xtalpls]
    type = CPKalidindiFullCoupledUpdate # pubilic CPStressFullCoupledUpdateBase
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
  [../] 

  [./deformation_energy]
    type = CompDefEnergy # compDefEnergy
    grain_tracker = grain_tracker

    output_properties = 'elastic_energy delastic_energy/dgr0  delastic_energy/dgr1 plastic_energy dplastic_energy/dgr0  dplastic_energy/dgr1'
    outputs = my_exodus
  [../]
  [./free_energy] # caluate chemcal free energy density and derivative 
    type = DerivativeParsedMaterial
    coupled_variables = 'gr0 gr1'
    material_property_names = 'mu gamma_asymm'
    expression = 'mu*( gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2) + 1.0/4.0'
    derivative_order = 1
    enable_jit = true # Enable just-in-time compilation of function expressions for faster evaluation

    outputs = my_exodus
  [../]  
[]

[Postprocessors]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
  [../]
  [./vonmises_stress]
    type = ElementAverageValue
    variable = vonmises_stress
  [../]
  # [./slip_increment]
  #   type = ElementAverageValue
  #   variable = slip_increment
  # [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 15
  nl_max_its = 7
  nl_rel_tol = 1.0e-7
  start_time = 0.0
  
  # num_steps = 3
  end_time = 30.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  [./my_checkpoint]
    file_base = ./cp_${my_filename}/out_${my_filename}
    # interval = 10
    type = Checkpoint
    additional_execute_on = 'INITIAL FINAL'
  [../] 
  [my_exodus]
    type = Nemesis
    file_base = ./ex_${my_filename}/out_${my_filename} 
    time_step_interval = 10
    additional_execute_on = 'INITIAL FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  [pgraph]
    type = PerfGraphOutput
    # interval = 10
    level = 2                     # Default is 1
    heaviest_branch = false        # Default is false
    heaviest_sections = 7         # Default is 0
    additional_execute_on = 'FINAL'  # Default is "final"
  []
  print_linear_residuals = false
[]
    