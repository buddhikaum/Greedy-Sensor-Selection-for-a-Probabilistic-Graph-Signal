#!/bin/bash 

#N=20
obs_frac=.5
ct=0
for obs_frac in .5
do
	ct=$((ct+1))
	for N in 100 200 
	do
		echo running for $N
		rm greedy_${N}.txt
		for i in `seq 1 10`
		do
			echo $i
			python GenerateCov.py $i $N
			python Algo1.py $i $N $obs_frac >> greedy_${N}.txt

		done
	done
done
