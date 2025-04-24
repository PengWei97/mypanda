mpiexec -n 32 ~/projects/mypanda/mypanda-opt -i inp_elastic_energy_and_EBSD.i > demo1_c1_e5p.log & 

# --recover  c1_layer1_initial_mesh1/out_c1_layer1_initial_mesh1_cp/LATEST > 02_recover_m1.log & 

# mpiexec -n 30 ~/projects/mypanda/mypanda-opt -i inp_p23_bm1_layer1.i --recover  c1_layer1_initial/out_c1_layer1_initial_cp/LATEST > 02_recover.log & 

# tar -cvf - ex_demo1_c1_e5p/*.e-s0022.* | pigz -9 -p 20 > ex_demo1_c1_e5p.tgz
