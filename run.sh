~/UGP/allocator/src.allocator.out 64 8
make
touch output.csv


for p in 16 36 49 64
do
mpicc -o f1 src.c -lm


for N in 256 1024 4096 16384 65536 262144 
do

for i in `seq 1 5` 
do
mpirun -np $p -hostfile hosts ./f1 $N 50 1>>output.csv
echo "Method_1 P=$p N=$N Completed"
done


for i in `seq 1 5` 
do
mpirun -np $p -hostfile hosts ./f1 $N 50 2>>output.csv
echo "Method_2 P=$p N=$N Completed"
done


for i in `seq 1 5` 
do
mpirun -np $p -hostfile hosts ./f1 $N 50 3>>output.csv
echo "Method_3 P=$p N=$N Completed"
done



done
done

python plot.py
