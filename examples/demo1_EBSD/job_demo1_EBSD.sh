mpiexec -n 32 ~/projects/mypanda/mypanda-opt -i inp2_elastic_energy_and_EBSD.i > demo1_c1_e5p.log & 

# --recover  c1_layer1_initial_mesh1/out_c1_layer1_initial_mesh1_cp/LATEST > 02_recover_m1.log & 

# mpiexec -n 30 ~/projects/mypanda/mypanda-opt -i inp_p23_bm1_layer1.i --recover  c1_layer1_initial/out_c1_layer1_initial_cp/LATEST > 02_recover.log & 

# tar -cvf - ex_3__17min_R3_l1a2_fe/*.e-s000[37].* | pigz -9 -p 20 > ex_3__17min_R3_l1a2_fe.tgz
