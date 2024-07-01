rm -rf csv* cp* ex*

mpiexec -np 10 ~/projects/mypanda/mypanda-opt -i m1_coupledCPPF_fullyCoupled.i > 01.log & 

# tar -cvf - ex_s3_fullCoupledCPPF/*.e* | pigz -9 -p 20 > ex_s3_fullCoupledCPPF.tgz