# C1_L13_iso_noStored: 算例1，模拟区域Layer13，只考虑gb isotropy，不考虑存储能驱动，my_high_mob = 5.0e-11, slurm-c1.log

my_filename = "C1_L13_iso_noStored"
my_filename2 = "C1_L13_iso_noStored"

my_s3_mob = 5.0e-14
my_s9_mob = 5.0e-14
my_high_mob = 5.0e-11

[GlobalParams]
  op_num = 15
  var_name_base = gr

  length_scale = 1.0e-6
  time_scale = 1.0

  grain_tracker = grain_tracker
  concurrent_recovery = false

  # rho_end_l2 = 2.3e13
  # rho_end_l3 = 8.5e13
  # rho_critical = 4.0e13
  # a_rho_l2 = 1.0e-4
  # a_rho_l3 = 1.0e-4
[]

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = Ni_10min_layer13_rho.inl
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
  [./PolyStoredEnergyEBSD] # PolyStoredEnergyEBSD->PolyStoredEnergyEBSDAction & PolycrystalElasticDrivingForce->PolycrystalElasticDrivingForceAction
    # ACSEDGPolyEBSD
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
      crystal_structure = cubic # hexagonal cubic 
      euler_angle_provider = ebsd_reader
    []
  []
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyWtStoredEnergy # GBAnisotropyMisori GBAnisotropyWtStoredEnergy
    T = 863.15
    wGB = 1.0
  
    GBsigma_HAGB = 1.09 # Ref1. https://www.sciencedirect.com/science/article/pii/S1359645421003165?via%3Dihub
    GBmob_HAGB = ${my_high_mob}

    Sigma3_sigma = 0.36 # Ref1
    Sigma9_sigma = 0.52 # Ref2. https://www.sciencedirect.com/science/article/pii/S1359645409003425
    Sigma3_mob = ${my_s3_mob}
    Sigma9_mob = ${my_s9_mob}

    euler_angle_provider = ebsd_reader
    crystal_structure = FCC

    gb_energy_anisotropy = false
    gb_mobility_anisotropy = false

    ebsd_reader = ebsd_reader
    stored_energy_mobility = false

    output_properties = 'L misori_angle twinning_type'
    outputs = my_exodus
  [../]
  [./deformed]
    type = DeformedGrainEBSDMaterial
    stored_factor = 5.0 # 5.0
    GNDs_provider = ebsd_reader

    output_properties = 'rho_eff' #  dstored_energy/dgr0
    outputs = my_exodus

    Elas_Mod = 7.30e10 # J/m^3
    Burg_vec = 1.24e-10 # m
  [../]
  # [./free_energy] # caluate chemcal free energy density and derivative 
  #   type = DerivativeParsedMaterial
  #   coupled_variables = 'gr0 gr1'
  #   material_property_names = 'mu gamma_asymm'
  #   expression = 'mu*( gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2) + 1.0/4.0'
  #   derivative_order = 1
  #   enable_jit = true # Enable just-in-time compilation of function expressions for faster evaluation

  #   outputs = my_exodus
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
  end_time = 600.0
  # num_steps = 10

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
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
  # [./my_checkpoint]
  #   file_base = ./${my_filename}/out_${my_filename}
  #   time_step_interval = 5
  #   type = Checkpoint
  #   additional_execute_on = 'INITIAL FINAL'
  # [../]  
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    type = Nemesis
    # time_step_interval = 5
    additional_execute_on = 'FINAL'
    sync_times = '1 2 3 4 5 6 7 8 9 10 15	17 20	23 25	30	35	40	45	50	55	60	65	70	75	80	85	90	95	100	105	110	115	120	125	130	135	140	145	150	155	160	165	170	175	180	185	190	195	200	205	210	215	220	225	230	235	240	245	250	255	260	265	270	275	280	285	290	295	300	305	310	315	320	325	330	335	340	345	350	355	360	365	370	375	380	385	390	395	400	405	410	415	420	425	430	435	440	445	450	455	460	465	470	475	480	485	490	495	500	505	510	515	520	525	530	535	540	545	550	555	560	565	570	575	580	585	590	595	600'
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
    heaviest_branch = false        # Default is false
    heaviest_sections = 7         # Default is 0
    additional_execute_on = 'FINAL'  # Default is "final"
  []
  print_linear_residuals = false
[]

# [Problem]
#   restart_file_base = c4_ly1_noStored_aniso/out_c4_ly1_noStored_aniso_cp/LATEST
#   force_restart = true
# []