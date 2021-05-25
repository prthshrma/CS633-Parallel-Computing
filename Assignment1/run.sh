#!/bin/sh

make

rm data.txt

for k in 16 36 49 64;do
	~/UGP/allocator/src/allocator.out $k 8
	for i in 256 1024 4096 16384 65536 262144 1048576;do
        	for j in `seq 1 5`; do
                	mpiexec -n $k -f hosts ./src $i 50
        	done
	done
done

python3 plot.py
