#!/bin/bash
var=$1
for i in `seq 1 $var`;
do
	sbatch ./run_batch_devSimNull.sh $i
done
