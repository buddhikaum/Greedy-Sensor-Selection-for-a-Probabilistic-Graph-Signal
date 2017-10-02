# Greedy-Sensor-Selection-for-a-Probabilistic-Graph-Signal
Implementation of the algorithms in
"Efficient Sensor Selection with Application to Time Varying Graphs",
in 2017 IEEE 7th International Workshop on Computational Advances in
Multi-Sensor Adaptive Processing (CAMSAP) (IEEE CAMSAP 2017)

First algorithm imporoves the time complexity of the 
existing greedy sensor selection algorithm in [1].

Next two algorithms for efficient selection of sensors for 
estimating a probabilistic graph signal. 
First algorithm (SelectOriginalOnly in SelectTimeVarying.py) provides 
an efficent algorithm for computing trace of the posterior covariance matrix
in ordere to select sensors using a greedy approach such that the 
mean squared error will be minimized. 
The second algorithm (main method in SelectTimeVarying.py) provides 
an algorithm to update the sensors when there is a change 
in a single edge in the original graph. Based on provious computations 
this algorithm selects sensors efficiently compared to selecting sensors 
from scratch. 

Laplacian.py generates the unnormalized laplacian from a random graph 
on which the covariance of the graph signal depeneds on.

The proposed algorithm has a time complexity of O(n^3) compared to 
previously published greedy algorithm by Liu et. al. [1].

Usage : 
# For algorithm 1 run the script file 
bash algo1_script.txt 
# For algorithms 2 and 3 in the paper 
python Laplacian.py [random seed] [number of vertices] [b]

#b is 1,-1 indicating whether we are removing an edge or adding an edge.

To compute sensors
python SelectTimeVarying.py [random seed] [number of vertices] [fraction observed] [b]

#fraction observed is a value less than 1. 

Eg. 

python Laplacian.py 1 400 1

python SelectTimeVarying.py 1 400 .5 1






[1] S. Liu, S. P. Chepuri, M. Fardad, E. Maşazade, G. Leus and P. K. Varshney, 
"Sensor Selection for Estimation with Correlated Measurement Noise," 
in IEEE Transactions on Signal Processing, 
vol. 64, no. 13, pp. 3509-3522, July1, 1 2016.
