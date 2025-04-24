# demo1_c1_IsoElaticTensor

my_filename = "demo1_c0_f6"

my_s3_mob = 1.0e-11
my_s9_mob = 1.0e-11
my_Hmob = 1.0e-11

[GlobalParams]
  op_num = 15
  var_name_base = gr

  length_scale = 1.0e-6
  time_scale = 1.0

  grain_tracker = grain_tracker
[]

[Mesh]
  [./ebsd_mesh]
    type = EBSDMeshGenerator
    filename = GNSNi_17min_R3_layer1a2_local2.inl
    pre_refine = 6 # Mesh can go two levels coarser than the EBSD grid
  [../]
  parallel_type = distributed
[]

[UserObjects]
  [./ebsd_reader]
    # Get Euler angles, coordinates, grain ID, phase ID, symmetry, GNDs from EBSD file
    type = EBSDReader
    custom_columns = 0
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
    expression = 'grain_tracker < 20'
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
  [./bnds]
  [../]
  [./unique_grains] # Grain ID
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices] # index of the etaX
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Phi]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./phi2]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./phi1]
    type = OutputEulerAngles
    variable = phi1
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'phi1'
  [../]
  [./Phi]
    type = OutputEulerAngles
    variable = Phi
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'Phi'
  [../]
  [./phi2]
    type = OutputEulerAngles
    variable = phi2
    euler_angle_provider = ebsd_reader
    
    output_euler_angle = 'phi2'
  [../]
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
    GBmob_HAGB = ${my_Hmob}

    Sigma3_sigma = 0.36 # Ref1
    Sigma9_sigma = 0.52 # Ref2. https://www.sciencedirect.com/science/article/pii/S1359645409003425
    Sigma3_mob = ${my_s3_mob}
    Sigma9_mob = ${my_s9_mob}

    euler_angle_provider = ebsd_reader
    crystal_structure = FCC

    gb_energy_anisotropy = false
    gb_mobility_anisotropy = false

    output_properties = 'L mu misori_angle twinning_type'
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

  start_time = 0.0
  num_steps = 3

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-2
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 1
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 1
  [../]
[]

[Outputs]
  [my_exodus]
    file_base = ./ex_${my_filename}/out_${my_filename} 
    type = Nemesis
    additional_execute_on = 'FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
  [./pgraph]
    type = PerfGraphOutput
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 5         # Default is 0
    execute_on = 'TIMESTEP_END FINAL'
  [../]
  print_linear_residuals = false
[]