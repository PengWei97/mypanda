# mkdir -p /home/pengwei/projects/mypanda/src/kernels/grain_growth/
# mkdir -p /home/pengwei/projects/mypanda/include/kernels/grain_growth/

# cp /home/pengwei/projects/panda/src/kernels/grain_growth/ACGrGrDeformationEnergyDrivingFC.C /home/pengwei/projects/mypanda/src/kernels/grain_growth/ACGrGrDeformationEnergyDrivingFC.C

# cp /home/pengwei/projects/panda/include/kernels/grain_growth/ACGrGrDeformationEnergyDrivingFC.h /home/pengwei/projects/mypanda/include/kernels/grain_growth/ACGrGrDeformationEnergyDrivingFC.h

# rm -rf cp_s1_fullCoupledCPPF csv_s1_fullCoupledCPPF ex_s1_fullCoupledCPPF
# mpiexec -np 10 ~/projects/mypanda/mypanda-opt -i m1_coupledCPPF_fullyCoupled.i > 01.log & 
