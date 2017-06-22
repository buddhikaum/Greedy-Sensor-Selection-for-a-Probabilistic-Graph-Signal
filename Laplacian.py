import networkx as nx
import sys
import numpy as np
from numpy import linalg as LA

#This script is for generating a covariance 
#matrix based on graph laplacian
rseed = int(sys.argv[1])
N = int(sys.argv[2])
b = int(sys.argv[3])
G = nx.binomial_graph(N,.5,seed=rseed)
A = nx.adjacency_matrix(G)
L = np.zeros((N,N))
L = nx.laplacian_matrix(G).todense()
edge1=0
edge2=1
if b==1:
	checkval = 0
else:
	checkval = -1
while L[edge1,edge2]!=checkval:
	
	rval = np.random.uniform(0,1,1)
	temp1 = int(rval*N)
	rval = np.random.uniform(0,1,1)
	temp2 = int(rval*N)
	#print temp1,temp2
	if temp1==temp2:
		continue
	edge1 = temp1
	edge2 = temp2
	#print L[edge1,edge2]

if b==1:
	print 'Edge to be added ',edge1,edge2
else:
	print 'Edge to be removed ',edge1,edge2
edge = np.array([edge1,edge2])

np.savetxt('edge.txt',edge,fmt='%d')

np.savetxt('cov.txt',LA.inv(L+.01*np.identity(N)))
rdiag = np.ones(N)
full_R = np.diag(rdiag)
np.savetxt('meas_cov.txt',full_R)
