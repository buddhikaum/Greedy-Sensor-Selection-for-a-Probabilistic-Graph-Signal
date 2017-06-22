# Greedy-Sensor-Selection-for-a-Probabilistic-Graph-Signal
Provides two algorithms for efficient selection of sensors for 
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

Laplacian.py generates the unnormalized laplacian from random graph 
on which the covariance of the graph signal depeneds on 

Usage : 
python Laplacian.py <random seed> <number of vertices> <b>

#b is 1,-1 indicating whether we are removing an edge or adding an edge.

To compute sensors
python SelectTimeVarying.py <random seed> <number of vertices> <fraction observed> <b>

#fraction observed is a value less than 1. 

Eg. 

python Laplacian.py 1 400 1

python SelectTimeVarying.py 1 400 .5 1
