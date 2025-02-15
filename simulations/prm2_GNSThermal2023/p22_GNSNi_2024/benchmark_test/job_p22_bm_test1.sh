# rm -rf ex_* csv_* p22_*
mpiexec -n 35 ~/projects/mypanda/mypanda-opt -i inp_p22_bm_test1.i > log_p22_case1_Iso.log & 

# --oversubscribe 
# tar -cvf - ex_p22_bm_case1_Iso/*.e-s0017.35.* | pigz -9 -p 20 > out_p22_bm_case1_Isos.tgz 