# Begin
1. case1_GBIsotropy
2. case2_GBAnisotropy
3. case3_GBIsotropy-StoredEnergy
4. case4_GBAnisotropy-StoredEnergy

# EX_14min
1. 2024.11.15: 基于大面积离位EBSD表征-14min，采用 `GNSNi2_14min_Stage3_local4a2` 执行 `case1_GBIsotropy_stage3_local4a2` 模拟；


# 错误
## recover

- 当采用 `PolycrystalColoringIC` 时后，recover 模拟时显示如上错误：
```bash
*** ERROR ***
/home/pw-moose/projects/mypanda/simulations/prm2_GNSThermal2023/p23_GNSNi_2024/ex_14min/inp_GNSNi_level1_590du.i:313.1:
The following error occurred in the Problem 'MOOSE Problem' of type FEProblem.

Initial conditions have been specified during a checkpoint restart, by IC object 'PolycrystalColoringIC' for variable 'PolycrystalColoringIC0'.
This is only allowed if you specify 'allow_initial_conditions_with_restart' to the [Problem], as initial conditions can override restarted fields

Abort(1) on node 23 (rank 23 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 23
```

- 当采用在input file 中设定如下时，
```bash
[Problem]
  solve = false
  restart_file_base = Case1_GBIsotropy_Stage3_Local4a2/out_Case1_GBIsotropy_Stage3_Local4a2_cp/0000
  allow_initial_conditions_with_restart = true
[]
```
出现下错误：
```bash
*** ERROR ***
RestartableEquationSystems::load(): Previously stored elements/nodes do not match the current element/nodes

Abort(1) on node 32 (rank 32 in comm 0): application called MPI_Abort(MPI_COMM_WORLD, 1) - process 32
```