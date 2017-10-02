
import numpy as np 
from numpy import linalg as LA
import sys
import time 



def GetMeasSet(N,M,R,Sigma_pre):

	Sr = Sigma_pre+R
	D,V = LA.eigh(Sr)
	a = D[0]/2
	Dnew = D-a
	S = np.dot(V,np.dot(np.diag(Dnew),V.T))
	Si = S #this is S_0^{-1},i=0 case

	avail_list = range(0,N)
	selected_list = []
	Ps = np.dot(Sigma_pre,LA.inv(S))
	Q = np.dot(Ps,Si)

	for i in range(0,M):

		SigmaDiag = np.zeros(N)
		crnt_max = 0
		max_idx = -1
		for j in avail_list:
			SigmaDiag[j] = np.dot(Q[:,j],Q[:,j].T)
		for j in avail_list:
			crnt_sensor = SigmaDiag[j]/(a+Si[j,j])
			if crnt_sensor>crnt_max:
				crnt_max = crnt_sensor
				max_idx = j
	
		#print i,max_idx,crnt_max
		avail_list.remove(max_idx)
		selected_list.append(max_idx)
	
		Q = Q - np.outer(Q[:,max_idx],Si[:,max_idx])/(a+Si[max_idx,max_idx])
		Si = Si - np.outer(Si[max_idx,:],Si[max_idx,:])\
											/(a+Si[max_idx,max_idx])

	return selected_list

rseed = int(sys.argv[1])
N = int(sys.argv[2])
obs_frac = float(sys.argv[3])
Sigma = np.genfromtxt('cov.txt')
full_R = np.genfromtxt('meas_cov.txt')
M=int(N*obs_frac)
t0 = time.time()
GetMeasSet(N,M,full_R,Sigma)
print 'time ',time.time()-t0
