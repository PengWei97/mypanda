# rm -rf ex_* csv_* test_*
mpiexec -n 35 ~/projects/mypanda/mypanda-opt -i inp_GNSNi_level1_590du.i > 01.log & 

# tar -cvf - ex_case1_GBIsotropy_stage3_local4a2_v3/*.e-s0021.* | pigz -9 -p 20 > ex_case1_GBIsotropy_stage3_local4a2_v3.tgz

# scp pw-moose@pw-moose-PowerEdge-T640:./*.tgz C:\Users\12481\Downloads