mpiexec -np 10 ~/projects/mypanda/mypanda-opt -i inp_GNSTi_AGG_level2a3_550du.i > 01.log &

# tar -cvf - csv_case3_noStored_v1/* | pigz -9 -p 20 > csv_case3_noStored_v1.tgz 
