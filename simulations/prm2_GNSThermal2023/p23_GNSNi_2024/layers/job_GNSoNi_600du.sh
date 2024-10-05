#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 2
#SBATCH -n 256

# rm -rf ex_* csv_* test_*
mpiexec --oversubscribe -n 35 ~/projects/mypanda/mypanda-opt -i inp_GNSNi_layer1_600du.i > slurm-c1.log & 

# tar -cvf - ex_c7_ly1_Stored_aniso/*.e-s0055.* ex_c7_ly1_Stored_aniso/*.e-s0075.* | pigz -9 -p 20 > ex_c7_ly1_Stored_aniso.tgz 

# tar -cvf - ex_c4_ly1_noStored_aniso/*.e-s*.* | pigz -9 -p 20 > ex_c4_ly1_noStored_aniso.tgz 
# tar -cvf - ex_c6_ly1_Stored_anisoRho/*.e-s0060.* | pigz -9 -p 20 > ex_c6_ly1_Stored_anisoRho.tgz 

# tar -cvf - ex_c8_ly1_Stored_iso/*.e-s0012.* | pigz -9 -p 20 > ex_c8_ly1_Stored_iso.tgz c