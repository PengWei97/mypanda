rm -rf ex_* csv_* case?_*

mpiexec -np 35 ~/projects/mypanda/mypanda-opt -i old_model_triple_junction_gg_anisotropy.i
# tail -f 01.log  > 01.log & 

# tar -cvf - ex_case3_triple_junction_gg_old/* | pigz -9 -p 20 > ex_case3_triple_junction_gg_old.tgz 
