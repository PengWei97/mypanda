# rm -rf ex_* csv_* test_*
mpiexec --oversubscribe -n 35 ~/projects/mypanda/mypanda-opt -i inp_GNSNi_level1a2_600du.i > 01.log & 

# tar -cvf - ex_test7_anisotropy/*.e-s0013.* | pigz -9 -p 20 > ex_test7_anisotropy.tgz 

# scp pw-moose@pw-moose-PowerEdge-T640:./*.tgz C:\Users\12481\Downloads