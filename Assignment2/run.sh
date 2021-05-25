#!/bin/sh

make

rm data.txt


for j in `seq 1 10`; do
	for nodes in 4 16;do
		for ppn in 1 8;do
			n=`expr $nodes / 4`
			python3 script.py 4 $n $ppn
			for D in 16 256 2048;do
				processes=`expr $nodes \* $ppn`
				mpirun -np $processes -f hostfile ./src $D $ppn
			done
		done
	done
done


python3 plot.py

