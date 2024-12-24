# Example: Simulation of grain growth in materials exhibiting GBAnisotropy material class;
# The initial condition includes two circular grain embedded in a matirx griains, with the aim to study the influence of GB anisotropy on the grain growth procress.

my_filename = 'case2_two_circular_gg_old'
my_wGB = 3.0
my_number_adaptivity = 3
# my_interval = 2

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 10
  xmin = 0
  xmax = 200
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
[]

[UserObjects]
  [./term]
    type = Terminator
    expression = 'gr2_area < 1000'
  [../]
[]
[Variables]
  [./gr0] # for matrix grains
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SpecifiedSmoothCircleIC # creates multiple SmoothCircles
      x_positions = '50.0  150.0'
      y_positions = '50.0  50.0'
      z_positions = '0.0    0.0'
      radii       = '40.0  40.0'
      invalue = 0.0
      outvalue = 1.0
      int_width = 3.0
    [../]
  [../]
  [./gr1] # the left circular grains
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 50.0
      y1 = 50.0
      radius = 40.0
      invalue = 1.0
      outvalue = 0.0
      int_width = 3.0
    [../]
  [../]
  [./gr2] # for the right circular grain
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 150.0
      y1 = 50.0
      radius = 40.0
      invalue = 1.0
      outvalue = 0.0
      int_width = 3.0
    [../]
  [../]
[]

# [ICs]
#   [./PolycrystalICs]
#     [./BicrystalCircleGrainIC] # SmoothCircleIC
#       radius = 40
#       x = 50
#       y = 50
#       int_width = ${my_wGB}
#     [../]
#   [../]
# []

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

# [BCs]
#   [./Periodic]
#     [./All]
#       auto_direction = 'x y'
#       variable = 'gr0 gr1'
#     [../]
#   [../]
# []

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy
    T = 600 # K

    # molar_volume_value = 7.11e-6 #Units:m^3/mol
    Anisotropic_GB_file_name = two_circular_grain_anisotropy.tex
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

  end_time = 20
  # num_steps = 10

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
