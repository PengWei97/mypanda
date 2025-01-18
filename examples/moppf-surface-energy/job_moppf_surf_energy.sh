rm -rf ex_* csv_* case?_* *_cp *.tgz *.log

mpiexec -np 30 ~/projects/mypanda/mypanda-opt -i surf_circular_matrix_gg_3d_split.i > 01_3d_split_NoAdaptivity_Newton.log & 
# /home/pw-moose/projects/mypanda/examples/moppf-surface-energy/surf_circular_matrix_gg_3d_split.i

# tar -cvf - ex_c4_circ_gg_3d/* | pigz -9 -p 20 > ex_c4_circ_gg_3d.tgz 
# tar -cvf - ex_case4_circular_gg/*.e-s000?.* ex_case4_circular_gg/*.e-s001[012345].* | pigz -9 -p 20 > ex_case4_circular_gg.tgz 
# ll ex_case4_circular_gg_noBounds/*.e-s*.02