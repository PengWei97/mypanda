# Example: Simulation of grain growth in materials exhibiting GBAnisotropy material class;
# The initial condition includes a triple junction, with the aim to study the influence of GB anisotropy on the grain growth procress.

my_filename = 'case3_triple_junction_gg_old'
my_wGB = 1.6
my_number_adaptivity = 3
# my_interval = 2

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 3
  var_name_base = gr
  wGB = ${my_wGB}
  length_scale = 1.0e-9
  time_scale = 1.0e-9

  v = 'gr0 gr1 gr2' # Names of the grains
  theta1 = 90 # Angle the first grain makes at the triple junction
  theta2 = 90 # Angle the second grain makes at the triple junction
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

# [UserObjects]
#   [./term]
#     type = Terminator
#     expression = 'gr2_area < 2000'
#   [../]
# []

[ICs]
  [gr0_IC]
    type = TricrystalTripleJunctionIC
    variable = gr0
    op_index = 1
  []
  [gr1_IC]
    type = TricrystalTripleJunctionIC
    variable = gr1
    op_index = 2
  []
  [gr2_IC]
    type = TricrystalTripleJunctionIC
    variable = gr2
    op_index = 3
  []
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

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy
    T = 600 # K

    # molar_volume_value = 7.11e-6 #Units:m^3/mol
    Anisotropic_GB_file_name = triple_junction_anisotropy.tex
    inclination_anisotropy = false # true

    output_properties = 'L mu gamma_asymm kappa_op' 
    outputs = my_exodus
  [../]
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
  [./gr2_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr2
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

[Executioner]
  type = Transient
  scheme = bdf2

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 20
  l_tol = 1e-4
  nl_max_its = 10
  nl_rel_tol = 1e-9

  # end_time = 20
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
    # time_step_interval = ${my_interval} # The interval at which time steps are output
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
