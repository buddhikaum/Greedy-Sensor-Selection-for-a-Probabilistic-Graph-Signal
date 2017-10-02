import numpy as np 
from numpy import linalg as LA
import sys
import time 



rseed = int(sys.argv[1])
N = int(sys.argv[2])
np.random.seed(rseed)
X = np.zeros((2,N))
for i in range(0,N):
	X[:,i] = np.random.multivariate_normal(np.zeros(2),
									np.identity(2))
Sigma = np.zeros((N,N))
for i in range(0,N):
	for j in range(0,N):
		norm_2 = LA.norm(X[:,i]-X[:,j])
		#print norm_2
		Sigma[i,j] = np.exp(-.2*norm_2)

X = np.zeros((2,N))
for i in range(0,N):
	X[:,i] = np.random.multivariate_normal(np.zeros(2),
									np.identity(2))
	X[:,i] = np.random.uniform(0,1,2)
full_R = np.zeros((N,N))
for i in range(0,N):
	for j in range(0,N):
		norm_2 = LA.norm(X[:,i]-X[:,j])
		full_R[i,j] = np.exp(-.1*norm_2)

np.savetxt('cov.txt',Sigma)
np.savetxt('meas_cov.txt',full_R)

