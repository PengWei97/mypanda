my_filename = "case4_recovery_v7"
my_filename2 = "case4_recovery_v72s"

my_tt1_mob = 0.792e-14
my_ct1_mob = 1.048e-14

[GlobalParams]
  op_num = 12
  var_name_base = gr

  length_scale = 1.0e-6
  time_scale = 1.0

  grain_tracker = grain_tracker
  concurrent_recovery = true

  rho_end_l2 = 2.3e13
  rho_end_l3 = 8.5e13
  rho_critical = 4.0e13
  a_rho_l2 = 1.0e-4
  a_rho_l3 = 1.0e-4
[]

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    # EBSD Data with GNDs during isothermal annealing with GNDs at 700℃ 
    filename = Ti550du_120min_level2a3_rho.inl
  []
  parallel_type = distributed
[]

[UserObjects]
  [./ebsd_reader]
    # Get Euler angles, coordinates, grain ID, phase ID, symmetry, GNDs from EBSD file
    type = EBSDReader
    custom_columns = 1
  [../]
  [./ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = jp
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    compute_var_to_feature_map = false
  [../]
  [./grain_tracker]
    type = GrainTrackerMerge
    threshold = 0.5
    connecting_threshold = 1.0e-2
    halo_level = 3
    flood_entity_type = ELEMENTAL
    polycrystal_ic_uo = ebsd
    compute_var_to_feature_map = true

    execute_on = 'INITIAL TIMESTEP_BEGIN'

    merge_grains_based_misorientaion = true
    euler_angle_provider = ebsd_reader
  [../]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 10'
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    [../]
  [../]
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  [phi1]
    order = CONSTANT
    family = MONOMIAL
  []
  [Phi]
    order = CONSTANT
    family = MONOMIAL
  []
  [phi2]
    order = CONSTANT
    family = MONOMIAL
  []
  [./ebsd_rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolyStoredEnergyEBSD] # PolyStoredEnergyEBSDAction
    # ACSEDGPolyEBSD
    GNDs_provider = ebsd_reader
  [../]
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  []
  [phi1]
    type = OutputEulerAngles
    variable = phi1
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'phi1'
  []
  [Phi]
    type = OutputEulerAngles
    variable = Phi
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'Phi'
  []
  [phi2]
    type = OutputEulerAngles
    variable = phi2
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'phi2'
  []
  [grain_rho]
    type = EBSDReaderPointDataAux
    variable = ebsd_rho # GNDs
    ebsd_reader = ebsd_reader
    data_name = 'CUSTOM0'
  []
[]

[Modules]
  [PhaseField]
    [EulerAngles2RGB]
      crystal_structure = hexagonal # hexagonal cubic 
      euler_angle_provider = ebsd_reader
    []
  []
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori
    T = 973.15
    wGB = 1.0
  
    GBsigma_HAGB = 0.9565
    GBmob_HAGB = 2.7e-14 

    TT1_sigma = 0.276 # 0.1019 0.3109 0.276
    CT1_sigma = 0.291 # 0.0616 0.1848 0.291
    TT1_mob = ${my_tt1_mob}
    CT1_mob = ${my_ct1_mob}

    euler_angle_provider = ebsd_reader

    gb_energy_anisotropy = true
    gb_mobility_anisotropy = true

    output_properties = 'L mu misori_angle'
    outputs = my_exodus
  [../]
  [./deformed]
    type = DeformedGrainEBSDMaterial
    stored_factor = 5.0 # 5.0
    GNDs_provider = ebsd_reader

    output_properties = 'rho_eff' # test_loc
    outputs = my_exodus
  [../]
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./bnd_length]
    type = GrainBoundaryArea
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor
    flood_counter = grain_tracker
    output_centroids = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold' #  -snes_type
  petsc_options_value = 'hypre boomeramg 31 0.7' # vinewtonrsls

  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves
  l_max_its = 10 # Max number of linear iterations
  nl_max_its = 8 # Max number of nonlinear iterations
  dtmin = 1.0e-4

  start_time = 0.0
  end_time = 50 # 10200 # 21600
  # num_steps = 12

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 3
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 2
  [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    time_step_interval = 10
    type = Checkpoint
    additional_execute_on = 'INITIAL FINAL'
  [../]  
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    type = Nemesis
    time_step_interval = 10
    additional_execute_on = 'FINAL'
    # sync_times = '10 50 100 200 600	1200	1800	2400	3000	3600	4200	4800	5400	6000	6600	7200	7800	8400	9000	9600	10200	10800	11400	12000	12600	13200	13800	14400	15000	15600	16200	16800	17400	18000	18600	19200	19800	20400	21000	21600'
    # sync_only = true
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  [pgraph]
    type = PerfGraphOutput
    time_step_interval = 10
    level = 2                     # Default is 1
    heaviest_branch = false        # Default is false
    heaviest_sections = 7         # Default is 0
    additional_execute_on = 'FINAL'  # Default is "final"
  []
  print_linear_residuals = false
[]