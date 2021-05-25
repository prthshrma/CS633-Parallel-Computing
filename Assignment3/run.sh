#!/bin/sh

make

rm output.txt

for nodes in 1 2;do
	for ppn in 1 2 4;do
		python3 script.py 2 $nodes $ppn
		processes=`expr $nodes \* $ppn`
		mpirun -np $processes -f hostfile ./src tdata.csv
	done
done
