# Example: Simulation of grain growth in materials exhibiting GBAnisotropy material class;
# The initial condition includes a circular grain embedded in a matirx griains, with the aim to study the influence of GB anisotropy on the grain growth procress.

my_filename = 'case4_polycrystal_gg_old'
my_wGB = 10.0
my_number_adaptivity = 3
my_interval = 2

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
  zmin = 0
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 6
  var_name_base = gr
  grain_num = 6

  wGB = ${my_wGB}
  length_scale = 1.0e-9
  time_scale = 1.0e-9
[]

[UserObjects]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 4'
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = jp
    int_width = ${my_wGB}
    rand_seed = 12
  [../]
  [./grain_tracker]
    type = GrainTracker
    flood_entity_type = elemental
    threshold = 0.1
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
  [../]
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
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
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./bnds_aux]
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
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy
    T = 600 # K

    # molar_volume_value = 7.11e-6 #Units:m^3/mol
    Anisotropic_GB_file_name = a_poly_grain_anisotropy.tex
    inclination_anisotropy = false # true

    output_properties = 'L mu gamma_asymm kappa_op' 
    outputs = my_exodus
  [../]
[]

[Postprocessors]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.1
  [../]
  [./avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
    execute_on = 'initial timestep_end'
  [../]
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
[]

[Executioner]
  type = Transient
  scheme = bdf2

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 10
  l_tol = 1e-4
  nl_max_its = 10
  nl_rel_tol = 1e-9

  end_time = 50
  # num_steps = 10

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 4
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
