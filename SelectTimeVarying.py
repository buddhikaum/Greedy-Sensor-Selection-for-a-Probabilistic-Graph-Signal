
import numpy as np 
import math
import scipy.optimize as opt 
from scipy import linalg as LA
import sys
import time 

def SelectStartingJ(Sigma,R,N,M,avail_list,selected_list):

	#This is the subrouting for Algo. 3
	#Sensor selectin starting from \tau
	rsigmasqrd = R[0,0]
	for i in range(len(selected_list),M):
		SigmaDiag = np.zeros(N)
		#for j in range(0,N):
		for j in avail_list:
			SigmaDiag[j] = np.dot(Sigma[j,:],Sigma[j,:].T)
		max_er = 0
		max_idx=0

		for j in avail_list:
			crntval = SigmaDiag[j]/(rsigmasqrd+Sigma[j,j])
			if crntval>max_er:
				max_er = crntval 
				max_idx = j
		avail_list.remove(max_idx)
		selected_list.append(max_idx)
		Sigma = Sigma - np.outer(Sigma[:,max_idx],Sigma[:,max_idx])\
									/(rsigmasqrd+Sigma[max_idx,max_idx])
	return selected_list,avail_list
def SelectOriginalOnly(Sigma,R,N,M):

	#Selecting Sensors From Scratch
	selected_list = []
	avail_list = range(0,N)
	rsigmasqrd = R[0,0]
	T=0
	for i in range(0,M):
		t1 = time.time()
		SigmaDiag = np.zeros(N)
		for j in avail_list:
			SigmaDiag[j] = np.dot(Sigma[j,:],Sigma[j,:].T)
		max_er = 0
		max_idx=0
		for j in avail_list:
			crntval = SigmaDiag[j]/(rsigmasqrd+Sigma[j,j])
			#print j,crntval	
			if crntval>max_er:
				max_er = crntval 
				max_idx = j
		avail_list.remove(max_idx)
		selected_list.append(max_idx)
		Sigma = Sigma - np.outer(Sigma[:,max_idx],Sigma[:,max_idx])\
									/(rsigmasqrd+Sigma[max_idx,max_idx])
	
	return selected_list

def SelectOriginal(Sigma,R,N,M):

	#Algo. 2 in paper. 
	#Selection sensor and saving computations 
	selected_list = []
	avail_list = range(0,N)
	second_trace = np.zeros(M)
	min_trace = np.zeros(M)
	second_diag = np.zeros(M)
	sigma_eig = np.zeros(M)
	Omega = np.zeros((M,N))
	rsigmasqrd = R[0,0]
	Pi = np.zeros((M,N))
	Gamma = np.zeros((M,N))
	Sdiag = np.zeros((M,N))
	for i in range(0,M):
		#np.savetxt('sigma_'+str(i)+'.txt',Sigma)
		SigmaDiag = np.zeros(N)
		for j in range(0,N):
			SigmaDiag[j] = np.dot(Sigma[j,:],Sigma[j,:].T)
		crntQ = np.zeros(N)
		max_er = 0
		max_idx=0
		for j in avail_list:
			crntval = SigmaDiag[j]/(rsigmasqrd+Sigma[j,j])
			crntQ[j] = SigmaDiag[j]
			Sdiag[i,j] = Sigma[j,j]
			if crntval>max_er:
				max_er = crntval
				max_idx = j
		for kk in range(0,N):
			Omega[i,kk] = crntQ[kk]
		Gamma[i,:] = Sigma[max_idx,:]
		avail_list.remove(max_idx)
		selected_list.append(max_idx)
		Sigma_pre = Sigma
		Pi[i,:] = np.dot(Sigma[max_idx,:],Sigma)
		Sigma = Sigma - np.outer(Sigma[:,max_idx],Sigma[:,max_idx])\
									/(rsigmasqrd+Sigma[max_idx,max_idx])
		
	return selected_list,Omega,Pi,Gamma,Sdiag

rseed = int(sys.argv[1])
np.random.seed(rseed)
N = int(sys.argv[2])
obs_frac = float(sys.argv[3])
b = int(sys.argv[4])
M=int(N*obs_frac)
Sigma = np.genfromtxt('cov.txt')
full_R = np.genfromtxt('meas_cov.txt')

#raw_input()
V = np.zeros((N,1))
rval = np.random.uniform(0,.2,1)
rval=1
edge = np.genfromtxt('edge.txt',dtype=np.int32)
#print edge
#print rval
V[edge[0]] = rval 
V[edge[1]] =-rval
Sigmahat = LA.inv(LA.inv(Sigma)+b*np.dot(V,V.T));
SigmaOld = Sigma
#Selecting sensors from scratch
t0 = time.time()
selected_list2 = SelectOriginalOnly(Sigmahat,full_R,N,M)
print 'Time for Selecting Sensors From Scratch',time.time()-t0
kk=0
[selected_list,Omega,Pi,Gamma,Sdiag] = SelectOriginal(Sigma,full_R,N,M)
while kk<M and selected_list[kk]==selected_list2[kk]:
	kk = kk+1
print 'Number of valid sensors from Original graph',kk
#Selecting sensors using Algo. 3
t0 = time.time();
#w_0 and v_0
v_t = np.dot(V.T,np.dot(Sigma,V))[0,0]
w_t = np.dot(V.T,np.dot(np.dot(Sigma,Sigma),V))[0,0]
#This is \phi_0
psi_t = np.dot(Sigma,V) 
psi_t = psi_t[:,0]
rsigmasqrd = full_R[0,0]
eta_t = np.dot(np.dot(Sigma,Sigma),V)[:,0] # This is \eta_0
updated_selection = []
updated_avail = range(0,N)
T=0
allVqdot  = np.zeros(N)
t0 = time.time()
for i in range(0,M):
	
	crnt_sensor = selected_list[i]
	crnt_diag = Gamma[i,crnt_sensor]
	allVqdotSqrd = np.multiply(psi_t,psi_t) # This is phi_t^2 elementwise
	qRqnew = Omega[i,:] -2*b*psi_t*eta_t/(1+b*v_t)+w_t*allVqdotSqrd/(1+b*v_t)**2
	new_sensr_val_vec = qRqnew/((rsigmasqrd+Sdiag[i,:])*(1+b*v_t)-b*allVqdotSqrd)
	crnt_sensor_val = new_sensr_val_vec[crnt_sensor]
	maxval = crnt_sensor_val
	maxidx = crnt_sensor
	for j in updated_avail :
		if j!=crnt_sensor:
			if maxval < new_sensr_val_vec[j]:
				maxval = new_sensr_val_vec[j]
				maxidx = j

	if maxidx!=crnt_sensor:
		print 'changed at ',i
		diagMat = np.zeros((N,N))
		for kk in range(0,i):
			diagMat[selected_list[kk],selected_list[kk]]=1
		Qupdated = LA.inv(LA.inv(SigmaOld)+diagMat+b*np.outer(V,V))
		updated_selection,updated_avail = SelectStartingJ(Qupdated,full_R,N,M,\
										updated_avail,updated_selection)
		break
	updated_avail.remove(crnt_sensor)
	updated_selection.append(crnt_sensor)
	eta_t = eta_t - Pi[i,:]*psi_t[crnt_sensor]/(rsigmasqrd+crnt_diag)\
			+(Omega[i,crnt_sensor]*psi_t[crnt_sensor]/(rsigmasqrd+crnt_diag)**2\
						-eta_t[crnt_sensor]/(rsigmasqrd+crnt_diag))*Gamma[i,:]
	psi_t = psi_t - psi_t[crnt_sensor]*Gamma[i,:]/(rsigmasqrd+crnt_diag)
	w_t = (eta_t[edge[0]]-eta_t[edge[1]])*rval
	v_t = (psi_t[edge[0]]-psi_t[edge[1]])*rval

print 'Time for Algo. 3 ',time.time()-t0
for i in range(0,M):
	if updated_selection[i]!=selected_list2[i]:
		print 'error with updated set'
		break
