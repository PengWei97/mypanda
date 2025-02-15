# p22_bm_test1: 验证基于MTEX+DREAM3D创建的xxx_rho.inl数据是否能在moose中模拟成功。
# p22_bm_case1_Iso: Level1a2_local1 + Isotropy

my_filename = "p22_bm_case1_Iso"
my_filename2 = "p22_bm_case1_Iso"

my_s3_mob = 5.0e-13
my_s9_mob = 5.0e-13

[GlobalParams]
  op_num = 15
  var_name_base = gr

  length_scale = 1.0e-6
  time_scale = 1.0

  grain_tracker = grain_tracker
  # concurrent_recovery = false

  # rho_end_l2 = 2.3e13
  # rho_end_l3 = 8.5e13
  # rho_critical = 4.0e13
  # a_rho_l2 = 1.0e-4
  # a_rho_l3 = 1.0e-4
[]

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = GNSNi_10min_Level1a2_local1_denoising_with_rho.inl
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

    merge_grains_based_misorientaion = false
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
  [unique_grains] # Grain ID
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices] # index of the etaX
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
  # [./PolyStoredEnergyEBSD] # PolyStoredEnergyEBSD->PolyStoredEnergyEBSDAction & PolycrystalElasticDrivingForce->PolycrystalElasticDrivingForceAction
  #   # ACSEDGPolyEBSD
  # [../]
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
      crystal_structure = cubic # hexagonal cubic 
      euler_angle_provider = ebsd_reader
    []
  []
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori # GBAnisotropyMisori GBAnisotropyWtStoredEnergy
    T = 863.15
    wGB = 1.0
  
    GBsigma_HAGB = 1.09 # Ref1. https://www.sciencedirect.com/science/article/pii/S1359645421003165?via%3Dihub
    GBmob_HAGB = 5.0e-12

    Sigma3_sigma = 0.36 # Ref1
    Sigma9_sigma = 0.52 # Ref2. https://www.sciencedirect.com/science/article/pii/S1359645409003425
    Sigma3_mob = ${my_s3_mob}
    Sigma9_mob = ${my_s9_mob}

    euler_angle_provider = ebsd_reader
    crystal_structure = FCC

    gb_energy_anisotropy = false
    gb_mobility_anisotropy = false

    # ebsd_reader = ebsd_reader
    # stored_energy_mobility = true

    output_properties = 'L mu misori_angle twinning_type'
    outputs = my_exodus
  [../]
  # [./deformed]
  #   type = DeformedGrainEBSDMaterial
  #   stored_factor = 5.0
  #   GNDs_provider = ebsd_reader

  #   output_properties = 'rho_eff dstored_energy/dgr0'
  #   outputs = my_exodus

  #   Elas_Mod = 7.30e10 # J/m^3
  #   Burg_vec = 1.24e-10 # m
  # [../]
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
  end_time = 1.0e6
  # num_steps = 3

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 2
  [../]
[]

[Outputs]
  [./my_checkpoint]
    file_base = ./${my_filename}/out_${my_filename}
    time_step_interval = 5
    type = Checkpoint
    additional_execute_on = 'INITIAL FINAL'
  [../]  
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    type = Nemesis
    # time_step_interval = 5
    additional_execute_on = 'FINAL'
    sync_times = '0.1 0.25 0.50 0.75 1.00 2.5 5 7.5 10 12.5 15 17.5 20 30 50 75 100 150 200 300 400 500 600	800 1000 1200 1600 1800 2000 2400 3000	3600 3720 3840 3960 4080 4200 4320 4440 4560 4680 4800 4920 5040 5160 5280 5400 5520 5640 5760 5880 6000 6120 6240 6360 6480 6600 6720        6840 6960 7080 7200'
    # sync_times = '1 2 5 10 30 50 75 100 150 200 300 400 500 600	800 1000 1200 1600 1800 2000	2400 3000	3600'
    sync_only = true
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  [pgraph]
    type = PerfGraphOutput
    time_step_interval = 5
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 5         # Default is 0
    execute_on = 'TIMESTEP_END FINAL'
  []
  print_linear_residuals = false
[]

# [Problem]
#   solve = false
#   restart_file_base = ./test_isotropy_v1/out_test_isotropy_v1_cp/0000
# []